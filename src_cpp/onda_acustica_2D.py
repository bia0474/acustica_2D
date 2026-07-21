
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

#----------------------------------
# plot one snap of the simulation
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_1000.bin", dtype=np.float32)

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

PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx1000.bin", dtype=np.float32)
PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz1000.bin", dtype=np.float32)

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

#-------------------------------------------
# plot the PVxz e the snapshot corresponding
#------------------------------------------

PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx1000.bin", dtype=np.float32)
PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz1000.bin", dtype=np.float32)

data = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_1000.bin", dtype=np.float32)

wavefield = data.reshape((nx, nz))

PVx = PVx.reshape((nx, nz))
PVz = PVz.reshape((nx, nz))

x = np.arange(nx) * dx
z = np.arange(nz) * dz

X, Z = np.meshgrid(x, z, indexing='ij') #indexa a malha igual em C

step = 50

plt.figure(figsize=(8,6))


img = plt.imshow(wavefield.T, cmap="seismic", origin="upper", extent=[0, nx * dx, nz * dz, 0], aspect="auto")

plt.quiver(X[::step, ::step], Z[::step, ::step], PVx[::step, ::step], -PVz[::step, ::step], color='black', pivot='tail', scale=100)

#width: espessura do corpo da seta
#headwidth: largura da ponta da seta
#headlength: comprimento da ponta da seta

plt.colorbar(img, label="Amplitude")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title('Snapshot + Vetores de Poynting')

plt.show()

#-------------------------------------------
# animation with Poynting vectors
#------------------------------------------
fig, ax = plt.subplots(figsize=(8,6))

step = 20

x = np.arange(nx) * dx
z = np.arange(nz) * dz

X, Z = np.meshgrid(x, z, indexing="ij")

# primeiro frame
wave = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_250.bin", dtype=np.float32).reshape(nx,nz)

PVx = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFx250.bin", dtype=np.float32).reshape(nx,nz)

PVz = np.fromfile("/home/processamento/acustica_2D/outputs/PoyntingVectorOFz250.bin",dtype=np.float32).reshape(nx,nz)

img = ax.imshow(wave.T, cmap="seismic", origin="upper", extent=[0,nx * dx,nz * dz,0], aspect="auto")

quiv = ax.quiver(X[::step,::step], Z[::step,::step], PVx[::step,::step], -PVz[::step,::step], color="black",pivot="tail", scale=110)

def update(n):

    wave = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_{n}.bin", dtype=np.float32).reshape(nx,nz)

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

first = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_500.bin", dtype=np.float32)
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

frames = range(250, 3750, 250)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=True)

plt.show()

#----------------------------------
# velocity model
#----------------------------------

vel = np.fromfile("/home/processamento/acustica_2D/src_cpp/Vp_camadas_501x501.bin", dtype=np.float32)

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

