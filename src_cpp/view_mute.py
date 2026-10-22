import numpy as np
import matplotlib.pyplot as plt
from numba import njit, prange

# ----------------------------------------------------------
# Leitura dos arquivos (equivalente a readReceivers / readSources)
# ----------------------------------------------------------
def read_receivers(receivers_file, nrec, nboundary):
    data = np.loadtxt(receivers_file, delimiter=",", skiprows=1, dtype=int)
    data = np.atleast_2d(data)[:nrec]
    return data[:, 1] + nboundary, data[:, 2] + nboundary

def read_sources(sources_file, nsource, nboundary):
    data = np.loadtxt(sources_file, delimiter=",", skiprows=1, dtype=int)
    data = np.atleast_2d(data)[:nsource]
    return data[:, 1] + nboundary, data[:, 2] + nboundary

# ----------------------------------------------------------
# Mute
# ----------------------------------------------------------
@njit(parallel=True)
def Mute(seismogram, shot, rec_x, rec_z, shot_x, shot_z, dt, shift, window, v0=1500.0):
    result = np.zeros_like(seismogram)
    Nt = seismogram.shape[0]
    Nrec = seismogram.shape[1]
    for rec in prange(Nrec):
        dist = np.sqrt((rec_z[rec] - shot_z[shot])**2 + (rec_x[rec] - shot_x[shot])**2)
        t1 = dist / v0 + shift
        t2 = t1 + window
        for i in range(Nt):
            t = i * dt
            if t < t1:
                result[i, rec] = 0.0
            elif t < t2:
                result[i, rec] = (t - t1) / (t2 - t1) * seismogram[i, rec]
            else:
                result[i, rec] = seismogram[i, rec]
    return result

# ----------------------------------------------------------
# Parâmetros (ajuste ao seu caso)
# ----------------------------------------------------------
base = "/home/processamento/acustica_2D"
receivers_file = f"{base}/inputs/receivers.csv"
sources_file   = f"{base}/inputs/sources.csv"

Nboundary = 100       # mesmo valor do C++
Nsource = 25         # 25 tiros, 1 receptor por tiro
Nt = 8001
dx = 10.0            # espaçamento da malha (m)
dt = 0.0005          # passo de tempo (s)
tlag = 0.0           # atraso da Ricker (s)

shift = 0.48
window = 0.1
v0 = 1500.0

# ----------------------------------------------------------
# Geometria (índices -> metros)
# ----------------------------------------------------------
rec_x, rec_z = read_receivers(receivers_file, Nsource, Nboundary)
shot_x, shot_z = read_sources(sources_file, Nsource, Nboundary)

rec_x_m, rec_z_m = rec_x * dx, rec_z * dx
shot_x_m, shot_z_m = shot_x * dx, shot_z * dx

# ----------------------------------------------------------
# Lê os 25 sismogramas (1 traça por tiro) e aplica o mute
# ----------------------------------------------------------
gather = np.zeros((Nt, Nsource), dtype=np.float32)
muted_gather = np.zeros_like(gather)

for s in range(Nsource):
    f = f"{base}/outputs/seismogram_shot{s}.bin"   # ajuste se o C++ começa em 1
    gather[:, s] = np.fromfile(f, dtype=np.float32)

    muted_gather[:, s:s+1] = Mute(
        gather[:, s:s+1], s,
        rec_x_m[s:s+1], rec_z_m[s:s+1],
        shot_x_m, shot_z_m,
        dt, shift + tlag, window, v0
    )

# tempo de corte (início da rampa) por par fonte-receptor
dist = np.sqrt((rec_z_m - shot_z_m)**2 + (rec_x_m - shot_x_m)**2)
travel_idx = (dist / v0 + tlag + shift) / dt

# ----------------------------------------------------------
# Plot 1: uma traça (antes x depois)
# ----------------------------------------------------------
shot_idx = 12
plt.figure()
plt.plot(gather[:, shot_idx], label="seismogram")
plt.plot(muted_gather[:, shot_idx], label="muted")
plt.axvline(travel_idx[shot_idx], color="r", ls="--", label="t1")
plt.legend()
plt.title(f"Tiro {shot_idx}")

# ----------------------------------------------------------
# Plot 2: CMP gather completo (original x muted)
# ----------------------------------------------------------
fig, ax = plt.subplots(1, 2, sharey=True, figsize=(10, 6))
for a, d, title in zip(ax, (gather, muted_gather), ("original", "muted")):
    a.imshow(d, aspect="auto", cmap="gray")
    a.plot(np.arange(Nsource), travel_idx, "r", lw=2, alpha=0.7)
    a.set_title(title)
plt.show()