
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#----------------------------------
# PARAMETERS
#----------------------------------

L = 5000
dx = 10
dz = 10
T = 2.0 ##total simulation time

nx = 501
nz = 501

nrec = 381
"""
#----------------------------------
# plot one snap of the simulation
#----------------------------------

data = np.fromfile("snapshot_900.bin", dtype=np.float32)

wavefield = data.reshape((nx, nz))

plt.figure(figsize=(8,6))

plt.imshow(wavefield.T, cmap="seismic", origin="upper", extent=[0, nx*dx, nz*dz, 0], aspect="auto")

plt.colorbar(label="Amplitude")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Wavefield Snapshot")

plt.show()

#----------------------------------
# animatiom 1D acustic wave
#----------------------------------

fig, ax = plt.subplots(figsize=(8,6))

first = np.fromfile("snapshot_100.bin", dtype=np.float32)
first = first.reshape((nx, nz))

img = ax.imshow(first.T, cmap="seismic", origin="upper", extent=[0, nx*dx, nz*dz, 0], aspect="auto", animated=True)

plt.colorbar(img)

ax.set_xlabel("x (m)")
ax.set_ylabel("z (m)")

ax.set_title("2D Acoustic Wave")

def update(frame):

    data = np.fromfile(f"snapshot_{frame}.bin", dtype=np.float32)

    wavefield = data.reshape((nx, nz))

    img.set_array(wavefield.T)

    return [img]

frames = range(100, 3900, 100)

ani = animation.FuncAnimation(fig, update, frames=frames, interval=100, blit=True)

plt.show()
"""
#----------------------------------
# velocity model
#----------------------------------

vel = np.fromfile("velocity.bin", dtype=np.float32)

vel = vel.reshape((nx, nz))

plt.figure(figsize=(8,6))

plt.imshow(vel.T, origin="upper", extent=[0, nx*dx, nz*dz, 0], aspect="auto")

plt.colorbar(label="Velocity (m/s)")

plt.xlabel("x (m)")
plt.ylabel("z (m)")

plt.title("Velocity Model")

plt.show()
"""
#----------------------------------
# Plot the sismogram
#----------------------------------

data = np.fromfile("/home/processamento/acustica_2D/seismogram.bin", dtype=np.float32)

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
"""