import argparse
from matplotlib.pyplot import figure, subplot, imshow, scatter, plot, axis, axvline, savefig
from numpy import array, convolve, ones, arange, abs, linalg
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
    verticies = []
    with open(pdb, "r") as file:
        for line in file:
            if line.startswith('ATOM'):
                temp_atom = line.split()
                atom_name = temp_atom[2]
                if atom_name in valid_atoms[atoms] and temp_atom[4] == chain:
                    atom_vertex = [float(temp_atom[6]), float(temp_atom[7]), float(temp_atom[8])]
                    verticies.append(atom_vertex)
    return array(verticies)


def moving_average(data, window_size):
    window = ones(window_size) / window_size
    smoothed_data = convolve(data, window, mode='valid')
    return smoothed_data


def find_large_valleys(data, prominence):
    inverted_data = -data  # Invert the data to treat valleys as peaks
    valleys, _ = find_peaks(inverted_data, width=verticies.shape[0], prominence=prominence)

    return valleys
    

verticies = extract_vertices(pdb_file, atoms, chain)

vertex_count = verticies.shape[0]

distance_matrix = linalg.norm(verticies[:, None] - verticies, axis=-1)

signal = distance_matrix.ravel()

wavelet = 'morl'
scales = arange(1, resolution)

coefficients, frequencies = pywt.cwt(signal, scales, wavelet)

data = abs(coefficients[-1])
data = moving_average(data, vertex_count * 10)
valleys = find_large_valleys(data, prominence)

cmap = 'turbo'

fig = figure(figsize=(16, 8), constrained_layout=True)
spec = fig.add_gridspec(2, 4)

ax0 = fig.add_subplot(spec[0, 0])
ax0.scatter(verticies[:, 0], verticies[:, 2], alpha=1, s=5)
ax0.plot(verticies[:, 0], verticies[:, 2], alpha=1, linewidth=1)
ax0.axis('equal')
ax0.set_xlabel('X')
ax0.set_ylabel('Y')
ax0.set_title('Spatial Arrangement')

ax1 = fig.add_subplot(spec[0, 1])
ax1.imshow(distance_matrix, cmap=cmap)
ax1.set_xlabel('Index')
ax1.set_ylabel('Index')
ax1.set_title('Distrance Matrix')

ax2 = fig.add_subplot(spec[0, 2])
im1 = ax2.imshow(abs(coefficients), aspect='auto', cmap=cmap, extent=[0, len(signal), frequencies[-1], frequencies[0]])
ax2.set_xlabel('Index')
ax2.set_ylabel('Scale/Frequency (ω)')
ax2.set_title('Wavelet Transform Spectrogram')

ax3 = fig.add_subplot(spec[0, 3])
ax3.plot(data)
ax3.margins(0)
ax3.set_xlabel('Convolved Index')
ax3.set_ylabel('Magnitude')
ax3.set_title(f"ω = {round(frequencies[-1], 4)}, n = {len(valleys)}")

for valley in valleys:
    ax3.axvline(x=valley, color='orange', linestyle='-', alpha=1)

ax4 = fig.add_subplot(spec[1, :])
ax4.plot(signal, linewidth=1)
ax4.margins(0)
ax4.set_xlabel('Index')
ax4.set_ylabel('Distance')
ax4.set_title('Flattened Distance Matrix')

savefig(pdb_file[:-4] + ".png", format="png")

end_time = time.time()
elapsed_time = end_time - start_time

print(f"n = {len(valleys)}, Time = {elapsed_time:.2f} seconds")
