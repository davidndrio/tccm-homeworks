import numpy as np
import matplotlib.pyplot as plt

# Función para leer los datos de las energías desde el archivo 'trajectories.xyz'
def read_trajectory(filename):
    steps = []
    KE = []
    PE = []
    TE = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            try:
                # Buscar las líneas que contienen la información de las energías
                if "Step" in line:
                    parts = line.split(":")
                    step = int(parts[0].split()[1])
                    energies = parts[1].split(",")
                    ke = float(energies[0].split('=')[1].strip())
                    pe = float(energies[1].split('=')[1].strip().replace('-1', 'E-1'))  # Ajustar formato si es necesario
                    te = float(energies[2].split('=')[1].strip())

                    steps.append(step)
                    KE.append(ke)
                    PE.append(pe)
                    TE.append(te)
            except ValueError as e:
                print(f"Error parsing line: {line.strip()} -> {e}")
                continue

    return np.array(steps), np.array(KE), np.array(PE), np.array(TE)

# Función para graficar las energías
def plot_energies(steps, KE, PE, TE):
    plt.figure(figsize=(10, 6))

    # Gráfica de KE vs Step
    plt.subplot(311)
    plt.plot(steps, KE, label="Kinetic Energy (KE)", color='r')
    plt.xlabel('Step')
    plt.ylabel('KE (J)')
    plt.title('Kinetic Energy vs Step')
    plt.grid(True)

    # Gráfica de PE vs Step
    plt.subplot(312)
    plt.plot(steps, PE, label="Potential Energy (PE)", color='g')
    plt.xlabel('Step')
    plt.ylabel('PE (J)')
    plt.title('Potential Energy vs Step')
    plt.grid(True)

    # Gráfica de TE vs Step
    plt.subplot(313)
    plt.plot(steps, TE, label="Total Energy (TE)", color='b')
    plt.xlabel('Step')
    plt.ylabel('TE (J)')
    plt.title('Total Energy vs Step')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

# Ruta del archivo 'trajectories.xyz'
filename = 'trajectories.xyz'

# Leer los datos de las energías
steps, KE, PE, TE = read_trajectory(filename)

# Verificar si se leyeron datos correctamente
if len(steps) == 0:
    print("Error: No valid energy data found in the file.")
else:
    # Graficar las energías
    plot_energies(steps, KE, PE, TE)

