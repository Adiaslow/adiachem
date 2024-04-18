import argparse
from matplotlib.pyplot import figure, subplot, imshow, scatter, plot, axis, axvline, savefig, subplots
from numpy import array, convolve, ones, arange, linalg, argmax
from numpy import sum as npsum
from numpy import abs as npabs
import pywt
from scipy.signal import find_peaks
import sys
import time

start_time = time.time()

parser = argparse.ArgumentParser(
                    prog='Count Number of Secondary Structures',
                    description='Takes in a PDB file, calculates the Morlet Wavelet Transform'
                    'of the flattened distance matrix, and finds local minima to estimate'
                    'the number of secondary structures.')
parser.add_argument('-i', '--input', dest='pdb_file', required=True, help='Input PDB file.')
parser.add_argument('-a', '--atoms', dest='atoms', default="alpha-carbons",
                    choices=["all", "backbone", "backbone-no-carbonyl", "alpha-carbons"],
                    help='Atoms to consider.')
parser.add_argument('-c', '--chain', dest='chain', default="A",
                    help='Which chain to consider.')
parser.add_argument('-r', '--resolution', dest='resolution', type=int, default=32,
                    help='The resolution of the wavelet transform.')
parser.add_argument('-f', '--frequency-index', dest='frequency_index', type=int, default=-1,
                    help='The index of the frequency to use for finding local minima.')
parser.add_argument('-p', '--prominence', dest='prominence', default=0.5, type=float,
                    help='The required prominence of a minima to be considered for counting.')

args = parser.parse_args()

pdb_file = args.pdb_file
atoms = args.atoms
chain = args.chain
resolution = args.resolution
frequency_index = args.frequency_index
prominence = args.prominence

valid_atoms = {
    "all": {"CA", "C", "N", "O"},
    "backbone": {"CA", "C", "N", "O"},
    "backbone-no-carbonyl": {"CA", "C", "N"},
    "alpha-carbons": {"CA"}
}


def extract_vertices(pdb, atoms, chain):
    vertices = []
    with open(pdb, "r") as file:
        for line in file:
            if line.startswith('ATOM'):
                temp_atom = line.split()
                atom_name = temp_atom[2]
                if atom_name in valid_atoms[atoms] and temp_atom[4] == chain:
                    atom_vertex = [float(temp_atom[6]), float(temp_atom[7]), float(temp_atom[8])]
                    vertices.append(atom_vertex)
    return array(vertices)


def moving_average(data, window_size):
    window = ones(window_size) / window_size
    smoothed_data = convolve(data, window, mode='valid')
    return smoothed_data


def find_large_valleys(data, prominence):
    valleys, _ = find_peaks(data, width=vertices.shape[0], prominence=prominence)

    return valleys


def find_largest_amplitude(coefficients):
    max_index = argmax([npsum(npabs(coef)) for coef in coefficients])
    return max_index

vertices = extract_vertices(pdb_file, atoms, chain)

vertex_count = vertices.shape[0]

distance_matrix = linalg.norm(vertices[:, None] - vertices, axis=-1)

signal = distance_matrix.ravel()

wavelet = 'morl'
scales = arange(1, resolution)

coefficients, frequencies = pywt.cwt(signal, scales, wavelet)

index_of_highest_amplitude = find_largest_amplitude(coefficients)
frequency_used = frequencies[index_of_highest_amplitude]
# print(index_of_highest_amplitude)
data = npabs(coefficients[index_of_highest_amplitude])

data = moving_average(data, vertex_count * 10)

max_amplitude = max(data)
# print(f"Max Amplitude = {max_amplitude}, Left = {data[0]}, Right = {data[-1]}")
if data[0] < max_amplitude * 0.5 and data[-1] < max_amplitude * 0.5:
    # print("Flipped")
    data = max_amplitude-data



valleys = find_large_valleys(-data, prominence)

cmap = 'turbo'
fig, axs = subplots(2, 2, figsize=(8, 8), constrained_layout=True)

# Plot 1: Spatial Arrangement
axs[0, 0].scatter(vertices[:, 0], vertices[:, 2], alpha=1, s=5)
axs[0, 0].plot(vertices[:, 0], vertices[:, 2], alpha=1, linewidth=1)
axs[0, 0].axis('equal')
axs[0, 0].set_xlabel('X')
axs[0, 0].set_ylabel('Z')
axs[0, 0].set_title('Spatial Arrangement')

# Plot 2: Distance Matrix
axs[0, 1].imshow(distance_matrix, cmap=cmap, extent=[0, vertex_count, 0, vertex_count])
axs[0, 1].set_xlabel('Index')
axs[0, 1].set_ylabel('Index')
axs[0, 1].set_title('Distance Matrix')

# Plot 3: Wavelet Transform Spectrogram
im1 = axs[1, 0].imshow(abs(coefficients), aspect='auto', cmap=cmap, extent=[0, len(signal), frequencies[index_of_highest_amplitude], frequencies[0]])
axs[1, 0].set_xlabel('Index')
axs[1, 0].set_ylabel('Scale/Frequency (ω)')
axs[1, 0].set_title('Wavelet Transform Spectrogram')

# Plot 4: Convolved Data
axs[1, 1].plot(data)
axs[1, 1].margins(0)
axs[1, 1].set_xlabel('Convolved Index')
axs[1, 1].set_ylabel('Magnitude')
axs[1, 1].set_title(f"ω = {round(frequency_used, 4)}, n = {len(valleys)}")

# Highlight valleys
for valley in valleys:
    axs[1, 1].axvline(x=valley, color='orange', linestyle='-', alpha=1)

savefig(pdb_file[:-4] + ".png", format="png")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"n = {len(valleys)}, Time = {elapsed_time:.2f} seconds")
