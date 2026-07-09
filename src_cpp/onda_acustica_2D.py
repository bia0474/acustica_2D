
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#----------------------------------
# PARAMETERS
#----------------------------------

L = T = dx = dz = dt = nrec = None

with open("/home/processamento/acustica_2D/inputs/parameters.txt", "r") as file:
    for linha in file:
        if linha.startswith("L ="):
            L = float(linha.split("=")[1])
        elif linha.startswith("T ="):
            T = int(linha.split("=")[1])
        elif linha.startswith("dx ="):
            dx = float(linha.split("=")[1])
        elif linha.startswith("dz ="):
            dz = float(linha.split("=")[1])
        elif linha.startswith("dt ="):
            dt = float(linha.split("=")[1])
        elif linha.startswith("nrec ="):
            nrec = int(linha.split("=")[1])

nx_complete_model = 501
nz_complete_model = 501
nx_snap = 381
nz_snap = 381
nt = int(T/dt) + 1

#----------------------------------
# plot one snap of the simulation
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_500.bin", dtype=np.float32)

wavefield = data.reshape((nx_snap, nz_snap))

plt.figure(figsize=(8,6))

plt.imshow(wavefield.T, cmap="seismic", origin="upper", extent=[0, nx_snap*dx, nz_snap*dz, 0], aspect="auto")

plt.colorbar(label="Amplitude")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Wavefield Snapshot")

plt.show()

#----------------------------------
# animatiom 1D acustic wave
#----------------------------------

fig, ax = plt.subplots(figsize=(8,6))

first = np.fromfile("/home/processamento/acustica_2D/outputs/snapshot_500.bin", dtype=np.float32)
first = first.reshape((nx_snap, nz_snap))

img = ax.imshow(first.T, cmap="seismic", origin="upper", extent=[0, nx_snap*dx, nz_snap*dz, 0], aspect="auto", animated=True)

plt.colorbar(img)

ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")

ax.set_title("2D Acoustic Wave")

def update(frame):

    data = np.fromfile(f"/home/processamento/acustica_2D/outputs/snapshot_{frame}.bin", dtype=np.float32)

    wavefield = data.reshape((nx_snap, nz_snap))

    img.set_array(wavefield.T)

    return [img]

frames = range(100, 3900, 100)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=True)

plt.show()

#----------------------------------
# velocity model
#----------------------------------

vel = np.loadtxt("/home/processamento/acustica_2D/inputs/velocityModel.csv", delimiter=",", skiprows=1)

vel = vel.reshape((nx_complete_model, nz_complete_model))

plt.figure(figsize=(8,6))

plt.imshow(vel.T, origin="upper", extent=[0, nx_complete_model*dx, nz_complete_model*dz, 0], aspect="auto")

plt.colorbar(label="Velocity (m/s)")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Velocity Model")

plt.show()

#----------------------------------
# Plot the sismogram
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/outputs/seismogram.bin", dtype=np.float32)

nt = data.size // nrec

seismogram = data.reshape((nrec, nt), order="C")

plt.figure(figsize=(10,8))

vmax = np.percentile(np.abs(seismogram), 80)

plt.imshow(seismogram.T, cmap="gray", aspect="auto", vmin=-vmax, vmax=vmax)

plt.xlabel("Receiver")
plt.ylabel("Time sample")

plt.title("Seismogram")

plt.colorbar()

plt.show()
