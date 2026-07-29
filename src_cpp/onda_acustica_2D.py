
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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

data = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_1500.bin", dtype=np.float32)

wavefield = data.reshape((nx, nz))

plt.figure(figsize=(8,6))

plt.imshow(wavefield.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

plt.colorbar(label="Amplitude")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Wavefield Snapshot")

plt.show()

#----------------------------------
# plot the PVxz
#----------------------------------

PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx1500.bin", dtype=np.float32)
PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz1500.bin", dtype=np.float32)

#PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX1500.bin", dtype=np.float32)
#PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ1500.bin", dtype=np.float32)

PVx = PVx.reshape((nx, nz))
PVz = PVz.reshape((nx, nz))

x = np.arange(nx) * dx
z = np.arange(nz) * dz

X, Z = np.meshgrid(x, z, indexing='ij')

step = 20

plt.figure(figsize=(8,6))

plt.quiver(X[::step, ::step], Z[::step, ::step], PVx[::step, ::step], -PVz[::step, ::step], color='black', pivot='tail', scale=80)

#width: espessura do corpo da seta
#headwidth: largura da ponta da seta
#headlength: comprimento da ponta da seta

plt.show()
'''
#-------------------------------------------------------
# plot the PVxz and PVOFxz to the snapshot corresponding
#-------------------------------------------------------

PVx_of = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx1500.bin", dtype=np.float32).reshape((nx, nz))
PVz_of = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz1500.bin", dtype=np.float32).reshape((nx, nz))

PVx_pv = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX1500.bin", dtype=np.float32).reshape((nx, nz))
PVz_pv = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ1500.bin", dtype=np.float32).reshape((nx, nz))

wavefield = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_1500.bin", dtype=np.float32).reshape((nx, nz))

x = np.arange(nx) * dx
z = np.arange(nz) * dz
X, Z = np.meshgrid(x, z, indexing='ij')

step = 30  # menos denso

#-------------------------------------------
# função de plotagem com cor por direção
#-------------------------------------------

def plot_panel(ax, Vx, Vz, title):
    ax.imshow(wavefield.T, cmap="gray", origin="upper", extent=[0, nx*dx, nz*dz, 0], interpolation="bilinear", aspect="auto")

    Xs = X[::step, ::step]
    Zs = Z[::step, ::step]

    Vx_s = Vx[::step, ::step]
    Vz_s = Vz[::step, ::step]

    # --- normalização para vetor unitário ---
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

#-------------------------------------------
# figura com 2 painéis
#-------------------------------------------

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

plot_panel(axes[0], PVx_pv, PVz_pv, "(a) Poynting vector (Yoon & Marfurt)")
plot_panel(axes[1], PVx_of, PVz_of, "(b) Optical flow (Horn-Schunck)")

plt.tight_layout()

plt.show()

'''
#-------------------------------------------
# animation with Poynting vectors
#------------------------------------------
fig, ax = plt.subplots(figsize=(8,6))

step = 8

x = np.arange(nx) * dx
z = np.arange(nz) * dz

X, Z = np.meshgrid(x, z, indexing="ij")

# primeiro frame
wave = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_250.bin", dtype=np.float32).reshape(nx,nz)

#PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX250.bin", dtype=np.float32).reshape(nx,nz)

#PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ250.bin",dtype=np.float32).reshape(nx,nz)

PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx250.bin", dtype=np.float32)

PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz250.bin", dtype=np.float32)

img = ax.imshow(wave.T, cmap="seismic", origin="upper", extent=[0,nx * dx,nz * dz,0], interpolation="bilinear", aspect="auto")

quiv = ax.quiver(X[::step,::step], Z[::step,::step], PVx[::step,::step], -PVz[::step,::step], color="black", pivot="mid", scale=30, width=0.003)

def update(n):

    wave = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVx = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionX{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVz = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorDirectionZ{n}.bin", dtype=np.float32).reshape(nx,nz)

    PVx = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorOFx{n}.bin", dtype=np.float32).reshape(nx,nz)
    
    PVz = np.fromfile(f"/home/processamento/acustica_2D/outputs/PoyntingVectorOFz{n}.bin", dtype=np.float32).reshape(nx,nz)

    img.set_data(wave.T)

    quiv.set_UVC(PVx[::step,::step], -PVz[::step,::step])

    return img, quiv

ani = animation.FuncAnimation(fig, update, frames=range(250, 3900, 250), interval=300, blit=True)

plt.colorbar(img, label="Amplitude")

plt.show()


#----------------------------------
# animatiom 1D acustic wave
#----------------------------------

fig, ax = plt.subplots(figsize=(8,6))

first = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_1500.bin", dtype=np.float32)
first = first.reshape((nx, nz))

img = ax.imshow(first.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto", animated=True)

plt.colorbar(img, label="Amplitude")

ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")

ax.set_title("2D Acoustic Wave")

def update(frame):

    data = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_{frame}.bin", dtype=np.float32)

    wavefield = data.reshape((nx, nz))

    img.set_array(wavefield.T)

    return [img]

frames = range(250, 31500, 250)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=True)

plt.show()
'''

#----------------------------------
# velocity model marmousi Juan
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
# velocity model
#----------------------------------

vel = np.fromfile("/home/processamento/acustica_2D/src_cpp/vp_marmousi-ii_shape_(2801, 13601)_dh25m_Nz141_Nx681.bin", dtype=np.float32)

vel = vel.reshape((nx, nz))

plt.figure(figsize=(8,6))

plt.imshow(vel.T, origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

plt.colorbar(label="Velocity (m/s)")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Velocity Model")

plt.show()

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

'''