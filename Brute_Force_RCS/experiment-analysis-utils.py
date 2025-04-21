# IMPORTANT: the code in this folder is not currently set up to store results of the experiments correctly
# this file used to be a jupyter notebook in another file structure, so naming conventions are off

# TODO
# - rewrite code to either set up a folder structure that allows experimentation or the
# programmer provides a spot to generate the graphs

import json
import os
import numpy as np
import matplotlib.pyplot as plt

from evaluation_utils import tvd_truedist_empdist, total_variation_distance, xeb_truedist_empdist_ideal, xeb_truedist_empdist_noisy

import os
import json

def compute_avg_xeb_varyingqubits(
    shots: int, 
    min_qubits: int, 
    max_qubits: int, 
    noiseRate: float, 
    gigaShots: int, 
    isLogn: bool, 
    subLabel: str, 
    output_dir: str
) -> None:
    """
    Computes and saves average XEB values for a range of qubit counts.

    Parameters:
        shots (int): Number of shots per experiment.
        min_qubits (int): Minimum number of qubits.
        max_qubits (int): Maximum number of qubits.
        noiseRate (float): Noise rate for simulation.
        gigaShots (int): Number of unique sets of random unitaries to sample over.
        isLogn (bool): If True, uses log(n) depth; otherwise, uses n depth.
        subLabel (str): Subfolder for organizing output files (include '/' if needed).
        output_dir (str): Directory where the output file should be saved.
    """
    # Generate filename
    depth_type = "logndepth" if isLogn else "ndepth"
    filename = f"XEB_{depth_type}_shots{shots}_noiseRate{int(noiseRate * 100):02d}_gigaShots{gigaShots}.json"

    # Ensure subLabel ends with a '/'
    if subLabel and not subLabel.endswith('/'):
        subLabel += '/'

    # Construct full output path
    full_path = os.path.join(output_dir, f"data/Noisy/XEB/{subLabel}")
    
    # Ensure the directory exists
    os.makedirs(full_path, exist_ok=True)
    
    file_path = os.path.join(full_path, filename)

    # Initialize an empty dictionary if the file doesn't exist
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            data = json.load(f)
    else:
        data = {}

    # Gather XEBs for qubit counts
    for qubit_count in range(min_qubits, max_qubits + 1):
        XEBs = []
        for i in range(gigaShots):
            # Adjust depth based on whether depth is logN
            depth = None if isLogn else qubit_count
            # Compute XEB of true distribution and noisy empirical distribution
            XEB = xeb_truedist_empdist_noisy(num_qubits=qubit_count, noise_rate=noiseRate, shots=shots, depth=depth)
            XEBs.append(XEB)

        data[str(qubit_count)] = XEBs

        # Save data into file
        with open(file_path, "w") as f:
            json.dump(data, f, indent=4)





def plot_avg_xeb_varyingqubits(shots: int, min_qubits: int, max_qubits: int, noiseRate: float, gigaShots: int, isLogn: bool, subLabel: str) -> None:
    """
    Plots the average XEB scores for different qubit counts given a file.

    Parameters:
        shots (int): Number of shots per experiment.
        min_qubits (int): Minimum number of qubits.
        max_qubits (int): Maximum number of qubits.
        noiseRate (float): Noise rate for simulation.
        gigaShots (int): Number of unique sets of random unitaries to sample over.
        isLogn (bool): If True, uses log(n) depth; otherwise, uses n depth.
        subLabel (str): Subfolder for organizing output files (include '/' if needed).
    """
    # Generate filename
    depth_type = "logndepth" if isLogn else "ndepth"
    filename = f"XEB_{depth_type}_shots{shots}_noiseRate{int(noiseRate * 100):02d}_gigaShots{gigaShots}.json"
    output_path = f"data/Noisy/XEB/{subLabel}{filename}"

    # Load data from JSON
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"File {output_path} not found!")

    with open(output_path, "r") as f:
        data = json.load(f)

    qubit_counts = []
    avg_xebs = []

    # Calculate average XEBs
    for qubit_count, xeb_list in data.items():
        qubit_counts.append(int(qubit_count))  # Convert keys to int
        avg_xebs.append(np.mean(xeb_list))  # Compute average XEB

    # Plot
    plt.figure(figsize=(8, 5))
    plt.plot(qubit_counts, avg_xebs, marker='o', linestyle='-', color='b', label="Avg XEB")

    # Labels and Title
    plt.xlabel("Qubit Count")
    plt.ylabel("XEB Score")
    plt.title(filename)  # Use filename as title
    plt.legend()
    plt.grid(True)
    plt.show()


def compute_avg_xeb_by_depth(shots: int, num_qubits: int, min_depth: int, max_depth: int, noiseRate: float, gigaShots: int, subLabel: str) -> None:
    """
    Computes and saves average XEB values for a range of circuit depths.

    Parameters:
        shots (int): Number of shots per experiment.
        num_qubits (int): Fixed number of qubits.
        min_depth (int): Minimum circuit depth.
        max_depth (int): Maximum circuit depth.
        noiseRate (float): Noise rate for simulation.
        gigaShots (int): Number of unique sets of random unitaries to sample over.
        subLabel (str): Subfolder for organizing output files (include '/' if needed).
    """
    NoiseFolder = "Noisy"
    if noiseRate == 0.0:
        NoiseFolder = "Ideal"

    # Sets up filename and path based on parameters
    filename = f"XEB_varyingdepth_numqubits{num_qubits}_noiseRate{int(noiseRate * 100):02d}_mindepth{min_depth}_maxdepth{max_depth}.json"
    output_path = f"data/{NoiseFolder}/XEB/{subLabel}{filename}"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Initialize an empty dictionary if the file doesn't exist
    if os.path.exists(output_path):
        with open(output_path, "r") as f:
            data = json.load(f)
    else:
        data = {}

    # Gather XEBs for different circuit depths
    for depth in range(min_depth, max_depth + 1):
        XEBs = []
        for i in range(gigaShots):
            # Compute XEB of true distribution and noisy empirical distribution
            XEB = xeb_truedist_empdist_noisy(num_qubits, noiseRate, shots, depth=depth)
            XEBs.append(XEB)

        data[str(depth)] = XEBs

        # Save data into file
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)


def plot_avg_xeb_by_depth(num_qubits: int, shots: int, min_depth: int, max_depth: int, noiseRate: float, gigaShots: int, subLabel: str) -> None:
    """
    Plots and saves the average XEB scores for different depths.
    """
    NoiseFolder = "Noisy" if noiseRate != 0.0 else "Ideal"
    
    filename = f"XEB_varyingdepth_numqubits{num_qubits}_noiseRate{int(noiseRate * 100):02d}_mindepth{min_depth}_maxdepth{max_depth}.json"
    output_dir = f"data/{NoiseFolder}/XEB/{subLabel}"
    output_path = os.path.join(output_dir, filename)
    
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"File {output_path} not found!")
    
    with open(output_path, "r") as f:
        data = json.load(f)
    
    depths = list(map(int, data.keys()))
    avg_xebs = [np.mean(xeb_list) for xeb_list in data.values()]
    
    plt.figure(figsize=(8, 5))
    plt.plot(depths, avg_xebs, marker='o', linestyle='-', color='b', label=f"Avg XEB Varying Depth with {num_qubits}")
    
    plt.xlabel("Depth")
    plt.ylabel("XEB Score")
    plt.title(filename)
    plt.legend()
    plt.grid(True)
    
    # Save the plot in the same folder as the JSON file
    plot_filename = filename.replace(".json", ".png")
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path)
    plt.close()
    
    print(f"Plot saved at {plot_path}")


def compute_avg_tvd(shots: int, min_qubits: int, max_qubits: int, noiseRate: float, gigaShots: int, isLogn: bool, subLabel: str = "") -> None:
    """
    Computes and saves TVD values for a range of qubit counts.

    Parameters:
        shots (int): Number of shots per experiment.
        min_qubits (int): Minimum number of qubits.
        max_qubits (int): Maximum number of qubits.
        noiseRate (float): Noise rate for simulation.
        gigaShots (int): Number of unique sets of random unitaries to sample over.
        isLogn (bool): If True, uses log(n) depth; otherwise, uses n depth.
        subLabel (str): Optional subdirectory label.
    """
    NoiseFolder = "Noisy"
    if noiseRate == 0.0:
        NoiseFolder = "Ideal"

  # Set the depth_type based on the value of isLogn
    depth_type = "logn" if isLogn else "n"

    # Sets up filename and path based on parameters
    filename = f"TVD_{depth_type}_numqubits{min_qubits}-{max_qubits}_noiseRate{int(noiseRate * 100):02d}.json"
    output_path = f"data/{NoiseFolder}/TVD/{subLabel}{filename}"

    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Initialize an empty dictionary if the file doesn't exist
    if os.path.exists(output_path):
        with open(output_path, "r") as f:
            data = json.load(f)
    else:
        data = {}

    # Gather TVDs for qubit counts
    for qubit_count in range(min_qubits, max_qubits + 1):
        TVDs = []
        for i in range(gigaShots):
            # Adjusts depth based on whether depth is logN
            depth = None if isLogn else qubit_count
            # Compute TVD of true distribution and noisy empirical distribution
            TVD = tvd_truedist_empdist(num_qubits=qubit_count, noise_rate=noiseRate, shots=shots, depth=depth)
            TVDs.append(TVD)

        data[str(qubit_count)] = TVDs

        # Save data into file
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)

def plot_avg_tvd(min_qubits, max_qubits, noiseRate, isLogn, subLabel=""):
    """
    Loads and plots TVD values for a range of qubit counts.

    Parameters:
    min_qubits (int): Minimum number of qubits.
    max_qubits (int): Maximum number of qubits.
    noiseRate (float): Noise rate for simulation.
    isLogn (bool): If True, uses log(n) depth; otherwise, uses n depth.
    subLabel (str): Optional subdirectory label.
    """

    # Ensure subLabel has a trailing slash if it's not empty
    subLabel = subLabel.rstrip("/") + "/" if subLabel else ""

    # Determine the folder based on the noise rate (either "Noisy" or "Ideal")
    NoiseFolder = "Noisy" if noiseRate > 0.0 else "Ideal"

    # Generate the filename using parameters, determining depth type based on isLogn
    depth_type = "logndepth" if isLogn else "ndepth"
    filename = f"TVD_{depth_type}_numqubits{min_qubits}-{max_qubits}_noiseRate{int(noiseRate*100):02d}.json"
    output_path = f"data/{NoiseFolder}/TVD/{subLabel}{filename}"

    # Load the data from the JSON file
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"File {output_path} not found!")

    with open(output_path, "r") as f:
        data = json.load(f)

    qubit_counts = []  # List to store qubit counts
    avg_tvds = []  # List to store the average TVD values

    # Calculate the average TVD values for each qubit count
    for qubit_count, tvd_list in data.items():
        qubit_counts.append(int(qubit_count))  # Convert the qubit count to an integer
        avg_tvds.append(np.mean(tvd_list))  # Compute the average TVD for the current qubit count

    # Plotting the results
    plt.figure(figsize=(8, 5))
    plt.plot(qubit_counts, avg_tvds, marker='o', linestyle='-', color='b', label="Avg TVDs")

    # Set labels, title, and display options
    plt.xlabel("Qubit Count")
    plt.ylabel("TVD Score")
    plt.title(filename)  # Use the filename as the plot title
    plt.legend()
    plt.grid(True)

    # Save the plot in the same directory as the JSON file
    plot_filename = filename.replace(".json", ".png")
    plot_path = os.path.join(os.path.dirname(output_path), plot_filename)
    plt.savefig(plot_path)

    plt.show()


def compute_avg_tvd_by_depth(num_qubits, shots, min_depth, max_depth, noiseRate, gigaShots, subLabel):
    """
    Computes and saves average TVD values for a range of circuit depths.

    Parameters:
    shots (int): Number of shots per experiment.
    num_qubits (int): Fixed number of qubits.
    min_depth (int): Minimum circuit depth.
    max_depth (int): Maximum circuit depth.
    noiseRate (float): Noise rate for simulation.
    gigaShots (int): Number of unique sets of random unitaries to sample over.
    subLabel (str): Subfolder for organizing output files (include '/' if needed).
    """

    # Determine the noise folder (either "Noisy" or "Ideal" based on the noise rate)
    NoiseFolder = "Noisy" if noiseRate > 0.0 else "Ideal"

    # Set up the filename and output path using the parameters
    filename = f"TVD_varyingdepth_numqubits{num_qubits}_noiseRate{int(noiseRate*100):02d}_mindepth{min_depth}_maxdepth{max_depth}.json"
    output_path = f"data/{NoiseFolder}/TVD/{subLabel}{filename}"

    # Create the necessary directory if it doesn't exist
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    # Initialize data dictionary if the file doesn't exist
    if os.path.exists(output_path):
        with open(output_path, "r") as f:
            data = json.load(f)
    else:
        data = {}

    # Gather TVD values for different circuit depths
    for depth in range(min_depth, max_depth + 1):
        TVDs = []
        for i in range(gigaShots):
            # Compute the TVD of the true distribution and noisy empirical distribution for the current depth
            TVD = tvd_truedist_empdist(num_qubits, noise_rate=noiseRate, shots=shots, depth=None)
            TVDs.append(TVD)

        data[str(depth)] = TVDs

        # Save the data into the output JSON file
        with open(output_path, "w") as f:
            json.dump(data, f, indent=4)


def plot_avg_tvd_by_depth(num_qubits, shots, min_depth, max_depth, noiseRate, gigaShots, subLabel):
    """
    Plots and saves the average TVD scores for different depths.
    """

    # Determine the noise folder (either "Noisy" or "Ideal" based on the noise rate)
    NoiseFolder = "Noisy" if noiseRate != 0.0 else "Ideal"

    # Generate the filename for the TVD values based on the parameters
    filename = f"TVD_varyingdepth_numqubits{num_qubits}_noiseRate{int(noiseRate*100):02d}_mindepth{min_depth}_maxdepth{max_depth}.json"
    output_dir = f"data/{NoiseFolder}/TVD/{subLabel}"
    output_path = os.path.join(output_dir, filename)

    # Load the data from the JSON file
    if not os.path.exists(output_path):
        raise FileNotFoundError(f"File {output_path} not found!")

    with open(output_path, "r") as f:
        data = json.load(f)

    # Extract depths and compute the average TVD values for each depth
    depths = list(map(int, data.keys()))
    avg_tvds = [np.mean(tvd_list) for tvd_list in data.values()]

    # Plotting the results
    plt.figure(figsize=(8, 5))
    plt.plot(depths, avg_tvds, marker='o', linestyle='-', color='b', label=f"Avg TVD Varying Depth with {num_qubits}")

    # Set labels, title, and display options
    plt.xlabel("Depth")
    plt.ylabel("TVD Score")
    plt.title(filename)
    plt.legend()
    plt.grid(True)

    # Save the plot in the same folder as the JSON file
    plot_filename = filename.replace(".json", ".png")
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path)
    plt.close()

    print(f"Plot saved at {plot_path}")
