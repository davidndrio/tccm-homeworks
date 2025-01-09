import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d

# Function to read energy data from 'trajectories.xyz'
def read_trajectory(filename):
    steps = []
    KE = []
    PE = []
    TE = []

    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            # Check if the line contains the expected "Energies (in J/mol)" part
            if "Energies (in J/mol)" in line:
                try:
                    # Extract energy values
                    energy_part = line.split("Energies (in J/mol):")[1].strip()
                    energy_values = energy_part.split(", ")
                    ke = float(energy_values[0].split("=")[1].strip())
                    pe = float(energy_values[1].split("=")[1].strip())
                    te = float(energy_values[2].split("=")[1].strip())

                    # Extract step information (from the same line)
                    step_part = line.split("Step=")[1].split(",")[0].strip()
                    step = int(step_part)

                    # Append to respective lists
                    steps.append(step)
                    KE.append(ke)
                    PE.append(pe)
                    TE.append(te)

                except IndexError:
                    print(f"Skipping line due to incorrect format: {line}")

    # Convert lists to numpy arrays for consistency
    return np.array(steps), np.array(KE), np.array(PE), np.array(TE)

# Function to smooth the data using a moving average filter
def smooth_data(data, window_size):
    return uniform_filter1d(data, size=window_size)

# Function to plot the energies
def plot_energies(steps, KE, PE, TE, smooth=False, window_size=50):
    if smooth:
        KE = smooth_data(KE, window_size)
        PE = smooth_data(PE, window_size)
        TE = smooth_data(TE, window_size)

    plt.figure(figsize=(10, 8))

    # Plot KE vs Step
    plt.subplot(311)
    plt.plot(steps, KE, label="Kinetic Energy (KE)", color='r')
    plt.xlabel('Step')
    plt.ylabel('KE (J)')
    plt.title('Kinetic Energy vs Step')
    plt.grid(True)
    plt.legend()

    # Plot PE vs Step
    plt.subplot(312)
    plt.plot(steps, PE, label="Potential Energy (PE)", color='b')
    plt.xlabel('Step')
    plt.ylabel('PE (J)')
    plt.title('Potential Energy vs Step')
    plt.grid(True)
    plt.legend()

    # Plot TE vs Step
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

    # Read the data from the file
    steps, KE, PE, TE = read_trajectory(trajectory_file)

    # Plot the smoothed energies
    plot_energies(steps, KE, PE, TE, smooth=True, window_size=50)

