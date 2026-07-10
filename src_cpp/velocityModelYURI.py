import numpy as np
import matplotlib.pyplot as plt

nz = 621
nx = 621

# Modelo de 3 camadas
camadas = np.zeros((nz, nx),dtype=np.float32)
camadas[:50, :] = 1500
camadas[50:75, :] = 2700
camadas[75:125, :] = 3000
camadas[125:250, :] = 4400
camadas[250:, :] = 5500

#gradiente de velocidades

ax = 0
b = 3000 

for i in range(75, 125):
    camadas[i, :] = ax + b
    ax += 10

# Parâmetros da falha 1
x0 = 50 # posição da falha
m0 = 2.2 # inclinação da falha
D0 = 25  # rejeito ao longo da falha

# Parâmetros da falha 2
x1 = 450 # posição da falha
m1 = -2.2 # inclinação da falha
D1 = -25  # rejeito ao longo da falha

# Vetores unitários paralelos às falhas
tx0 = 1 / np.sqrt(1 + m0**2)
tz0 = m0 / np.sqrt(1 + m0**2)

tx1 = 1 / np.sqrt(1 + m1**2)
tz1 = m1 / np.sqrt(1 + m1**2)

dx0 = int(round(D0 * tx0))
dz0 = int(round(D0 * tz0))

dx1 = int(round(D1 * tx1))
dz1 = int(round(D1 * tz1))

X, Z = np.meshgrid(np.arange(nx), np.arange(nz))

falha0 = Z - m0*(X - x0) < 0
falha1 = Z - m1*(X - x1) < 0

X[falha0] -= dx0
Z[falha0] -= dz0
X[falha1] -= dx1
Z[falha1] -= dz1

mask = ((X >= 0) & (X < nx) &
        (Z >= 0) & (Z < nz))

camadas[mask] = camadas[Z[mask], X[mask]]

camadas[50:100, :] = 2700

# Visualização
plt.figure(figsize=(12, 6))
img = plt.imshow(camadas,vmax=5500, aspect='auto', cmap='jet')
cbar = plt.colorbar(img)
cbar.ax.invert_yaxis()
cbar.ax.set_title('Vp (m/s)', pad=20)
cbar.set_ticks([1500, 2700, 3000, 3500, 4500, 5500])
cbar.set_ticklabels(['Camada 1 = 1500m/s','Camada 2 = 2700m/s','Camada 3 (inicio) = 3000m/s', 'Camada 3 (final) = 3500m/s','Camada 4 = 4500m/s','Camada 5 = 5500m/s'])
plt.title('Modelo com Falha Normal')
plt.xlabel('distância (m)')
plt.ylabel('profundidade (m)')
plt.show()

camadas.tofile(f"Vp_camadas_{nz}x{nx}.bin")