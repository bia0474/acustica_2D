import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

base = "/home/processamento/acustica_2D"

# ----------------------------------------------------------
# Constantes do directWaveMute (C++) -- manter IDENTICAS ao C++
# (no C++ elas estao fixas dentro da funcao, nao vem do parameters.txt)
# ----------------------------------------------------------
V_DIRECT = np.float32(1500.0)
SHIFT = np.float32(0.3)
WINDOW = np.float32(0.11)

# ----------------------------------------------------------
# Parametros: lidos do mesmo parameters.txt do C++
# ----------------------------------------------------------
def read_parameters(fname):
    p = {}
    with open(fname) as f:
        for line in f:
            if "=" in line:
                k, v = line.split("=", 1)
                p[k.strip()] = v.strip()
    return p

par = read_parameters(f"{base}/inputs/parameters.txt")

nt = int(par["nt"])
dt = np.float32(par["dt"])
dx = np.float32(par["dx"])
dz = np.float32(par["dz"])
Nboundary = int(par["Nboudary"])   # a chave no parameters.txt tem a grafia do C++
Nsource = int(par["Nsource"])
nrec = int(par["nrec"])

# ----------------------------------------------------------
# Geometria (mesma logica de readReceivers / readSources)
# ----------------------------------------------------------
def read_xz(fname, n, nb):
    df = pd.read_csv(fname)                       # pula o cabecalho, como no C++
    x = df.iloc[:n, 1].to_numpy().astype(int) + nb
    z = df.iloc[:n, 2].to_numpy().astype(int) + nb
    return x, z

rec_x, rec_z = read_xz(f"{base}/inputs/{os.path.basename(par['receivers_file'])}", nrec, Nboundary)
sx, sz = read_xz(f"{base}/inputs/{os.path.basename(par['sources_file'])}", Nsource, Nboundary)

assert nrec >= Nsource, "o C++ usa receivers[shot]; precisa de nrec >= Nsource"

# ----------------------------------------------------------
# t1 e t2 por tiro, em float32 como no C++
# ----------------------------------------------------------
dz_rec = (rec_z[:Nsource] - sz).astype(np.float32) * dz
dx_rec = (rec_x[:Nsource] - sx).astype(np.float32) * dx
dist = np.sqrt(dz_rec * dz_rec + dx_rec * dx_rec).astype(np.float32)

t1 = (dist / V_DIRECT + SHIFT).astype(np.float32)
t2 = (t1 + WINDOW).astype(np.float32)

# ----------------------------------------------------------
# Mute: replica directWaveMute (float32, t = it*dt)
# ----------------------------------------------------------
def direct_wave_mute(trace, t1_s, t2_s, dt):
    out = trace.astype(np.float32).copy()
    t = np.arange(out.size, dtype=np.float32) * dt
    out[t < t1_s] = 0.0
    ramp = (t >= t1_s) & (t < t2_s)
    out[ramp] *= (t[ramp] - t1_s) / (t2_s - t1_s)
    return out

# ----------------------------------------------------------
# Leitura dos gathers gravados pelo C++
# ----------------------------------------------------------
def read_cmp(fname, nt):
    if not os.path.exists(fname) or os.path.getsize(fname) == 0:
        raise RuntimeError(f"Arquivo vazio ou inexistente: {fname}. Rode a modelagem direta de novo.")
    n_floats = os.path.getsize(fname) // 4
    if n_floats % nt != 0:
        raise ValueError(f"{n_floats} amostras nao sao divisiveis por nt={nt}.")
    n_pairs = n_floats // nt
    return np.fromfile(fname, dtype=np.float32).reshape(n_pairs, nt)   # (shot, tempo)

gather = read_cmp(f"{base}/outputs/cmp_gather.bin", nt)               # sem mute (C++)
gather_mute_cpp = read_cmp(f"{base}/outputs/cmp_gather_mute.bin", nt)  # com mute (C++)
n_pairs = gather.shape[0]
assert n_pairs == Nsource, f"gather tem {n_pairs} tiros, parameters.txt diz Nsource={Nsource}"

# mute em Python aplicado ao gather SEM mute
gather_mute_py = np.zeros_like(gather)
for s in range(Nsource):
    gather_mute_py[s] = direct_wave_mute(gather[s], t1[s], t2[s], dt)

# ----------------------------------------------------------
# Comparacao numerica: Python x C++
# ----------------------------------------------------------
diff = gather_mute_cpp - gather_mute_py
print(f"max |diff| (C++ - Python): {np.abs(diff).max():.3e}")
print(f"amostras diferentes: {np.count_nonzero(diff)} de {diff.size}")
bad_shots = np.where(np.abs(diff).max(axis=1) > 0)[0]
print("tiros com diferenca:", bad_shots) #número de amostras diferentes, os tiros com diferença e a amostra onde ela ocorre.
if bad_shots.size:
    s = bad_shots[0]
    idx = np.where(diff[s] != 0)[0]
    print(f"tiro {s}: amostras {idx.min()}..{idx.max()} (t1 = {t1[s] / dt:.2f} amostras, t2 = {t2[s] / dt:.2f})")

# ----------------------------------------------------------
# Offsets e ordenacao (mesmo esquema da sua plotagem)
# ----------------------------------------------------------
pairs = pd.read_csv(f"{base}/inputs/pares_geometry.csv")
offsets = pairs["offset_m"].values[:n_pairs]      # assume mesma ordem dos tiros (linha i = tiro i)

print(f"max ||offset_m| - dist|: {np.abs(np.abs(offsets) - dist).max():.3f} m") #Confere se |offset_m| do pares_geometry.csv bate com a distância calculada da geometria.

order = np.argsort(offsets)
offsets_sorted = offsets[order]
t = np.arange(nt) * float(dt)

# meia-largura de traco, para as curvas t1/t2 caírem no centro de cada traco
half = (offsets_sorted.max() - offsets_sorted.min()) / (n_pairs - 1) / 2.0
extent = [offsets_sorted.min() - half, offsets_sorted.max() + half, t.max(), t.min()]

vclip = np.percentile(np.abs(gather), 98)         # mesmo clip em todas as figuras

def plot_cmp(g, title, curves=False):
    data_sorted = g.T[:, order]                    # (tempo, shot) ordenado por offset
    plt.figure(figsize=(9, 6))
    plt.imshow(data_sorted, aspect="auto", cmap="gray", extent=extent, vmin=-vclip, vmax=vclip)
    if curves:
        plt.plot(offsets_sorted, t1[order], "r--", lw=1.5, label="t1 (inicio da rampa)")
        plt.plot(offsets_sorted, t2[order], "c--", lw=1.5, label="t2 (fim da rampa)")
        plt.legend(loc="lower right")
    plt.xlabel("Offset (m)")
    plt.ylabel("Tempo (s)")
    plt.title(title)
    plt.colorbar(label="Amplitude")
    plt.tight_layout()

# ----------------------------------------------------------
# Plot 1: uma traca (antes x depois), eixo em segundos
# ----------------------------------------------------------
shot_idx = 12
plt.figure(figsize=(9, 4))
plt.plot(t, gather[shot_idx], label="sem mute (C++)")
plt.plot(t, gather_mute_cpp[shot_idx], label="mute C++")
plt.plot(t, gather_mute_py[shot_idx], "--", label="mute Python")
plt.axvline(t1[shot_idx], color="r", ls="--", label="t1")
plt.axvline(t2[shot_idx], color="c", ls="--", label="t2")
plt.xlabel("Tempo (s)")
plt.ylabel("Amplitude")
plt.title(f"Tiro {shot_idx} (offset = {offsets[shot_idx]:.0f} m)")
plt.legend()
plt.tight_layout()

# ----------------------------------------------------------
# Plot 2-5: gathers completos do CMP
# ----------------------------------------------------------
plot_cmp(gather, "Sismograma completo (CMP) sem MUTE", curves=True)
plot_cmp(gather_mute_cpp, "Sismograma completo (CMP) com MUTE - C++")
plot_cmp(gather_mute_py, "Sismograma completo (CMP) com MUTE - Python")
plot_cmp(diff, "Diferenca (C++ - Python)")

plt.show()