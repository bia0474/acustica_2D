import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.stats import pearsonr
from matplotlib.colors import Normalize
from scipy.ndimage import laplace

#----------------------------------
# PARAMETERS
#----------------------------------

parameters = {}

with open("/home/processamento/acustica_2D/inputs/parameters.txt") as file:
    for line in file:
        key, value = line.strip().split(" = ")
        parameters[key] = value

nx = int(parameters["nx"])
nz = int(parameters["nz"])
nx_abc = int(parameters["nx_abc"])
nz_abc = int(parameters["nz_abc"])
nt = int(parameters["nt"])

dx = float(parameters["dx"])
dz = float(parameters["dz"])
dt = float(parameters["dt"])

nrec = int(parameters["nrec"])
Nboudary = int(parameters["Nboudary"])

'''
#----------------------------------
# plot one snap of the simulation
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_500.bin", dtype=np.float32)

wavefield = data.reshape((nx, nz))

plt.figure(figsize=(8,6))

plt.imshow(wavefield.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

plt.colorbar(label="Amplitude")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Wavefield Snapshot")

plt.show()
'''

#----------------------------------
# plot two snapshots side by side
#----------------------------------

data_fwd = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_fwd_3000.bin",dtype=np.float32)

data_back = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_back_3000.bin",dtype=np.float32)

wavefield_fwd = data_fwd.reshape((nx, nz))

wavefield_back = data_back.reshape((nx, nz))

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Snapshot forward
im1 = axes[0].imshow(wavefield_fwd.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

axes[0].set_xlabel("x (m)")

axes[0].set_ylabel("z (m)")

axes[0].set_title("Forward Wavefield")

fig.colorbar(im1, ax=axes[0], label="Amplitude")

# Snapshot backward
im2 = axes[1].imshow(wavefield_back.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

axes[1].set_xlabel("x (m)")

axes[1].set_ylabel("z (m)")

axes[1].set_title("Backward Wavefield")

fig.colorbar(im2, ax=axes[1], label="Amplitude")

plt.tight_layout()

plt.show()

'''
#--------------------------------------------------------
# plot the PVxz and PVOFxz to the snapshot corresponding
#--------------------------------------------------------

PVx_of = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx500.bin", dtype=np.float32).reshape((nx, nz))
PVz_of = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz500.bin", dtype=np.float32).reshape((nx, nz))

PVx_pv = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX500.bin", dtype=np.float32).reshape((nx, nz))
PVz_pv = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ500.bin", dtype=np.float32).reshape((nx, nz))

wavefield = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_500.bin", dtype=np.float32).reshape((nx, nz))

x = np.arange(nx) * dx
z = np.arange(nz) * dz
X, Z = np.meshgrid(x, z, indexing='ij')

step = 30  # menos denso

# função de plotagem com cor por direção

def plot_panel(ax, Vx, Vz, title):

    ax.imshow(wavefield.T, cmap="gray", origin="upper", extent=[0, nx*dx, nz*dz, 0], interpolation="bilinear", aspect="auto")

    Xs = X[::step, ::step]
    Zs = Z[::step, ::step]

    Vx_s = Vx[::step, ::step]
    Vz_s = Vz[::step, ::step]

    # normalização para vetor unitário
    norm = np.sqrt(Vx_s**2 + Vz_s**2)
    norm[norm == 0] = np.nan  # evita divisão por zero; NaN faz o quiver ignorar essas setas

    Vxs = Vx_s / norm
    Vzs = -Vz_s / norm  # mesmo sinal invertido que você já usava

    # separa por direção: vermelho = descendo (Vz>0 no eixo original), verde = subindo
    mask_down = Vz_s > 0
    mask_up   = ~mask_down

    ax.quiver(Xs[mask_down], Zs[mask_down], Vxs[mask_down], Vzs[mask_down], color='red', pivot='mid', scale=40, width=0.003)
    ax.quiver(Xs[mask_up], Zs[mask_up], Vxs[mask_up], Vzs[mask_up], color='lime', pivot='mid', scale=40, width=0.003)

    ax.set_title(title)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")

# figura com 2 painéis

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

plot_panel(axes[0], PVx_pv, PVz_pv, "(a) Poynting vector (Yoon & Marfurt)")

plot_panel(axes[1], PVx_of, PVz_of, "(b) Optical flow (Horn-Schunck)")

plt.tight_layout()

plt.show()


# ---------------------------------------------------------
# diferença angular (arccos do produto escalar normalizado)
# ---------------------------------------------------------

eps = 1e-12

mag_pv = np.sqrt(PVx_pv**2 + PVz_pv**2)
mag_of = np.sqrt(PVx_of**2 + PVz_of**2)

amp_thr = 0.01 * np.max(np.abs(wavefield))
mask = np.abs(wavefield) < amp_thr

mag_pv_safe = np.where(mag_pv > eps, mag_pv, np.nan)
mag_of_safe = np.where(mag_of > eps, mag_of, np.nan)

cos_theta = (PVx_pv/mag_pv_safe)*(PVx_of/mag_of_safe) + (PVz_pv/mag_pv_safe)*(PVz_of/mag_of_safe)
cos_theta = np.clip(cos_theta, -1.0, 1.0)

angle_diff = np.degrees(np.arccos(cos_theta))
angle_diff[mask] = np.nan

# --------------------------------------------------
# plot + estatística resumo
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(angle_diff.T, cmap="inferno", origin="upper", extent=[0, nx*dx, nz*dz, 0], vmin=0, vmax=90)

ax.contour(X, Z, wavefield, levels=[0.3*np.max(wavefield)], colors="cyan", linewidths=0.8)

fig.colorbar(im, ax=ax, label="Diferença angular (graus)")

ax.set_title("Diferença angular: Poynting (Y&M) vs Optical Flow (H-S)")

ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")

plt.tight_layout()

plt.show()

print(f"Erro angular mediano: {np.nanmedian(angle_diff):.2f} graus")
print(f"Erro angular médio:   {np.nanmean(angle_diff):.2f} graus")

#---------------------------------------------
# animation with Poynting vectors (Direction)
#---------------------------------------------

fig, ax = plt.subplots(figsize=(8,6))

step = 30

x = np.arange(nx) * dx
z = np.arange(nz) * dz

X, Z = np.meshgrid(x, z, indexing="ij")

def normalize_vectors(Vx, Vz, eps=1e-12):
    mag = np.sqrt(Vx**2 + Vz**2)
    mag[mag < eps] = eps
    return Vx / mag, Vz / mag

def get_direction_colors(V):
    # V > 0 -> seta aponta para cima na tela -> verde
    # V < 0 -> seta aponta para baixo na tela -> vermelho
    return np.where(V > 0, "black", "green").flatten()

# primeiro frame
wave = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_100.bin", dtype=np.float32).reshape(nx,nz)

PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX100.bin", dtype=np.float32).reshape(nx,nz)

PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ100.bin", dtype=np.float32).reshape(nx,nz)

PVx, PVz = normalize_vectors(PVx, PVz)

img = ax.imshow(wave.T, cmap="seismic", origin="upper", extent=[0,nx * dx,nz * dz,0], interpolation="bilinear", aspect="auto")

Vplot = -PVz[::step,::step]
quiv = ax.quiver(X[::step,::step], Z[::step,::step], PVx[::step,::step], Vplot, color=get_direction_colors(Vplot), pivot="mid", scale=30, width=0.003)

def update(n):

    wave = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVx = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVz = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVx, PVz = normalize_vectors(PVx, PVz)

    img.set_data(wave.T)

    Vplot = -PVz[::step,::step]
    quiv.set_UVC(PVx[::step,::step], Vplot)
    quiv.set_color(get_direction_colors(Vplot))

    return img, quiv

ani = animation.FuncAnimation(fig, update, frames=range(100, 4000, 100), interval=900, blit=True)

ani.save("/home/processamento/acustica_2D/outputs/poynting_direction.gif", writer=animation.PillowWriter(fps=0.8))

plt.colorbar(img, label="Amplitude")

plt.show()


#--------------------------------------------------
# animation with Poynting vectors using Optical Flow
#--------------------------------------------------

fig2, ax2 = plt.subplots(figsize=(8,6))

# primeiro frame
wave2 = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_100.bin", dtype=np.float32).reshape(nx,nz)

PVx_OF = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx100.bin", dtype=np.float32).reshape(nx,nz)

PVz_OF = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz100.bin", dtype=np.float32).reshape(nx,nz)

PVx_OF, PVz_OF = normalize_vectors(PVx_OF, PVz_OF)

img2 = ax2.imshow(wave2.T, cmap="seismic", origin="upper", extent=[0,nx * dx,nz * dz,0], interpolation="bilinear", aspect="auto")

Vplot2 = -PVz_OF[::step,::step]
quiv2 = ax2.quiver(X[::step,::step], Z[::step,::step], PVx_OF[::step,::step], Vplot2, color=get_direction_colors(Vplot2), pivot="mid", scale=30, width=0.003)

def update2(n):

    wave2 = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVx_OF = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorOFx{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVz_OF = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorOFz{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVx_OF, PVz_OF = normalize_vectors(PVx_OF, PVz_OF)

    img2.set_data(wave2.T)

    Vplot2 = -PVz_OF[::step,::step]
    quiv2.set_UVC(PVx_OF[::step,::step], Vplot2)
    quiv2.set_color(get_direction_colors(Vplot2))

    return img2, quiv2

ani2 = animation.FuncAnimation(fig2, update2, frames=range(100, 4000, 100), interval=900, blit=True)

ani2.save("/home/processamento/acustica_2D/outputs/poynting_of.gif", writer=animation.PillowWriter(fps=0.8))

plt.colorbar(img2, label="Amplitude")

plt.show()
'''
#------------------------------------
# animatiom 2D acustic wave - backward
#------------------------------------

fig, ax = plt.subplots(figsize=(8,6))

frames = range(6000, 10, -10)

clip = max(np.percentile(np.abs(np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_back_{f}.bin", dtype=np.float32)), 99.5)
    for f in frames[::40]
)

first = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_back_10.bin", dtype=np.float32)
first = first.reshape((nx, nz))

img = ax.imshow(first.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto", animated=True, vmin=-clip, vmax=clip)

plt.colorbar(img, label="Amplitude")

ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")

ax.set_title("2D Acoustic Wave - Backward")

def update(frame):

    data = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_back_{frame}.bin", dtype=np.float32)

    wavefield = data.reshape((nx, nz))

    img.set_array(wavefield.T)
    return [img]

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=True)

plt.show()

#------------------------------------
# animatiom 2D acustic wave - farward
#------------------------------------

fig, ax = plt.subplots(figsize=(8,6))

first = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_fwd_10.bin", dtype=np.float32)
first = first.reshape((nx, nz))

img = ax.imshow(first.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto", animated=True)

plt.colorbar(img, label="Amplitude")

ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")

ax.set_title("2D Acoustic Wave")

def update(frame):

    data = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_fwd_{frame}.bin", dtype=np.float32)

    wavefield = data.reshape((nx, nz))

    img.set_array(wavefield.T)
    return [img]

frames = range(10, 6000, 10)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=True)

plt.show()


#----------------------------------
# velocity model 
#----------------------------------

vel = np.fromfile("/home/processamento/acustica_2D/inputs/velocityModel.bin", dtype=np.float32)
vel = vel.reshape((nx, nz))  # mesma ordem row-major usada no C++ (c[i * nz_abc + j])

plt.figure(figsize=(8,6))

plt.imshow(vel.T, origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

plt.colorbar(label="Velocity (m/s)")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Velocity Model")

plt.show()

'''
#----------------------------------
# velocity model marmousi juan
#----------------------------------

vel = np.fromfile("/home/processamento/acustica_2D/src_cpp/vp_marmousi-ii_shape_(2801, 13601)_dh10m_Nz351_Nx851.bin", dtype=np.float32)

vel = vel.reshape((nz, nx))

plt.figure(figsize=(8,6))

plt.imshow(vel, origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

plt.colorbar(label="Velocity (m/s)")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Velocity Model")

plt.show()
'''
#----------------------------------
# Plot the sismogram
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/outputs/seismogram.bin", dtype=np.float32)

seismogram = data.reshape((nrec, nt), order="C")

plt.figure(figsize=(10,8))

vmax = np.percentile(np.abs(seismogram), 99)

plt.imshow(seismogram.T, cmap="gray", aspect="auto", vmin=-vmax, vmax=vmax)

plt.xlabel("Receiver")
plt.ylabel("Time sample")

plt.title("Seismogram")

plt.colorbar()

plt.show()

#----------------------------------
# Plot the sismogram WITH MUTE
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/outputs/seismogram_mute.bin", dtype=np.float32)

seismogram = data.reshape((nrec, nt), order="C")

plt.figure(figsize=(10,8))

vmax = np.percentile(np.abs(seismogram), 99)

plt.imshow(seismogram.T, cmap="gray", aspect="auto", vmin=-vmax, vmax=vmax)

plt.xlabel("Receiver")
plt.ylabel("Time sample")

plt.title("Seismogram - with mute")

plt.colorbar()

plt.show()

#----------------------------------
# PLOT THE MIGRATED IMAGE
#----------------------------------

image = np.fromfile("/home/processamento/acustica_2D/outputs/image.bin", dtype=np.float32)

image = image.reshape(nx, nz)

plt.figure(figsize=(6, 8))

vmax = np.percentile(np.abs(image), 99) if np.any(image) else 1.0

plt.imshow(image.T, cmap="gray", aspect="auto", vmin=-vmax, vmax=vmax)

plt.xlabel("x (m)")

plt.ylabel("z (m)")

plt.title("Migrated Image (RTM)")

plt.colorbar(label="Amplitude")

plt.show()

#----------------------------------
# Carrega os dados 
#----------------------------------

vel = np.fromfile("/home/processamento/acustica_2D/inputs/velocityModel.bin", dtype=np.float32)

vel = vel.reshape((nx, nz)) 

image = np.fromfile("/home/processamento/acustica_2D/outputs/image.bin", dtype=np.float32)

image = image.reshape((nx, nz))

#----------------------------------
# Filtro Laplaciano 
#----------------------------------

image_filt = laplace(image)

vmax_filt = np.percentile(np.abs(image_filt), 99)

plt.figure(figsize=(6, 8))

plt.imshow(image_filt.T, cmap="gray", aspect="auto", vmin=-vmax_filt, vmax=vmax_filt)

plt.xlabel("x (m)")

plt.ylabel("z (m)")

plt.title("Imagem migrada — após filtro Laplaciano")

plt.colorbar(label="Amplitude")

plt.show()

#----------------------------------
# Posições da fonte e dos receptores
#----------------------------------

src_ix, src_iz = 250, 100

src_x = src_ix * dx
src_z = src_iz * dz

rx_init = 60
rx_end = 440
nrec = 381

rec_x = np.arange(rx_init, rx_end + 1) * dx
rec_z = 60 * dz

#----------------------------------
# imagem sobreposta com modelo
#----------------------------------

vmax = np.percentile(np.abs(image_filt), 99) if np.any(image_filt) else 1.0

fig3, ax3 = plt.subplots(figsize=(8, 8))

ax3.imshow(vel.T, origin="upper", extent=[0, nx * dx, nz * dz, 0], cmap="PuOr", alpha=0.4, aspect="auto")

norm = Normalize(vmin=-vmax, vmax=vmax)

cmap_img = plt.get_cmap("Greys")

rgba_img = cmap_img(norm(image_filt.T))

alpha_img = np.clip((np.abs(image_filt.T) / vmax) ** 0.5, 0, 1)

rgba_img[..., 3] = alpha_img

ax3.imshow(rgba_img, origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

ax3.scatter(rec_x, np.full_like(rec_x, rec_z), marker="v", color="green", s=20, label="Receivers", zorder=5)

ax3.scatter([src_x], [src_z], marker="*", color="yellow", s=220, edgecolor="k", linewidth=0.6, label="Sources", zorder=6)

ax3.set_xlabel("Distance (m)")

ax3.set_ylabel("Depth (m)")

ax3.set_title("RTM - após filtro Laplaciano")

ax3.legend(loc="upper right", fontsize=9, framealpha=0.9)

sm = plt.cm.ScalarMappable(cmap=cmap_img, norm=norm)

sm.set_array([])

fig3.colorbar(sm, ax=ax3, fraction=0.046, pad=0.04)

plt.tight_layout()

plt.show()