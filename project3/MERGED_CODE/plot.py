import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

# Función para leer los datos de las energías desde el archivo 'trajectories.xyz'
def read_trajectory(filename):
    steps = []
    KE = []
    PE = []
    TE = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if "Step" in line:
                parts = line.split(":")
                step = int(parts[0].split()[1])
                energies = parts[1].split(",")
                ke = float(energies[0].split('=')[1].strip())
                pe = float(energies[1].split('=')[1].strip())
                te = float(energies[2].split('=')[1].strip())

                steps.append(step)
                KE.append(ke)
                PE.append(pe)
                TE.append(te)

    return np.array(steps), np.array(KE), np.array(PE), np.array(TE)

# Función para suavizar los datos usando un filtro de promedio móvil
def smooth_data(data, window_size):
    return uniform_filter1d(data, size=window_size)

# Función para graficar las energías
def plot_energies(steps, KE, PE, TE, smooth=False, window_size=50):
    if smooth:
        KE = smooth_data(KE, window_size)
        PE = smooth_data(PE, window_size)
        TE = smooth_data(TE, window_size)

    plt.figure(figsize=(10, 8))

    # Gráfica de KE vs Step
    plt.subplot(311)
    plt.plot(steps, KE, label="Kinetic Energy (KE)", color='r')
    plt.xlabel('Step')
    plt.ylabel('KE (J)')
    plt.title('Kinetic Energy vs Step')
    plt.grid(True)
    plt.legend()

    # Gráfica de PE vs Step
    plt.subplot(312)
    plt.plot(steps, PE, label="Potential Energy (PE)", color='b')
    plt.xlabel('Step')
    plt.ylabel('PE (J)')
    plt.title('Potential Energy vs Step')
    plt.grid(True)
    plt.legend()

    # Gráfica de TE vs Step
    plt.subplot(313)
    plt.plot(steps, TE, label="Total Energy (TE)", color='g')
    plt.xlabel('Step')
    plt.ylabel('TE (J)')
    plt.title('Total Energy vs Step')
    plt.grid(True)
    plt.legend()

    plt.tight_layout()
    plt.show()

# Main
if __name__ == "__main__":
    trajectory_file = "trajectories.xyz"

    # Leer los datos del archivo
    steps, KE, PE, TE = read_trajectory(trajectory_file)

    # Graficar las energías suavizadas
    plot_energies(steps, KE, PE, TE, smooth=True, window_size=50)

