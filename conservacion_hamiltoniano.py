import numpy as np
import matplotlib.pyplot as plt

# Leer los datos del archivo
datos = np.loadtxt("hamiltoniano.txt")

# Extraer la columna "t" y la columna "H'"
t = datos[:, 0]
H = datos[:, 1]

# Crear el gráfico
fig, ax = plt.subplots()
ax.plot(t, H)

# Establecer las etiquetas de los ejes
ax.set_xlabel("$t$ (s)")
ax.set_ylabel("$H'$ (1/s$^2$)")

# Establecer los límites del eje y
ax.set_ylim(1E-15,1E-11)

# Mostrar y guardad el gráfico
plt.savefig('hamiltoniano.png', dpi=300, bbox_inches='tight')
plt.show()
