#  rmsd_calculation.py -i input_pdb -c1 A -r reference_pdb -c2 A 
import argparse
from numpy import sum as npsum
from numpy import array, mean, sqrt
from scipy.linalg import svd, orthogonal_procrustes
from scipy.spatial.transform import Rotation
import sys
import time

start_time = time.time()

parser = argparse.ArgumentParser(
                    prog='Count Number of Secondary Structures',
                    description='Takes in a PDB file, calculates the Morlet Wavelet Transform' 
                    'of the flattened distance matrix, and finds local minima to estimate'
                    'the number of secondary structures.')

parser.add_argument('-i0', '--input0', dest='reference_pdb', required=True, help='Input PDB file.')
parser.add_argument('-i1', '--input1', dest='test_pdb', required=True, help='Input PDB file.')
parser.add_argument('-a', '--atoms', dest='atoms', default="alpha-carbons",
                    choices=["all", "backbone", "backbone-no-carbonyl", "alpha-carbons"],
                    help='Atoms to consider.')
parser.add_argument('-c0', '--chain0', dest='reference_chains', required=True, type=list, default="A",
                    help='Which chain to consider for the reference pdb.')
parser.add_argument('-c1', '--chain1', dest='test_chains', required=True, type=list, default="A",
                    help='Which chain to consider for the test pdb.')
parser.add_argument('-m', '--alignment-method', dest='alignment_method', required=True, default="kabsch",
					choices=["kabsch", "lsf", "svd", "quaternion"], help='Which chain to consider for the test pdb.')

args = parser.parse_args()

reference_pdb = args.reference_pdb
test_pdb = args.test_pdb
atoms = args.atoms
reference_chains = args.reference_chains
test_chains = args.test_chains
alignment_method = args.alignment_method

valid_atoms = {
    "all": {"CA", "C", "N", "O"},
    "backbone": {"CA", "C", "N", "O"},
    "backbone-no-carbonyl": {"CA", "C", "N"},
    "alpha-carbons": {"CA"}
}


def extract_vertices(pdb, atoms, chains):
    verticies = []
    with open(pdb, "r") as file:
        for line in file:
            if line.startswith('ATOM'):
                temp_atom = line.split()
                atom_name = temp_atom[2]
                if atom_name in valid_atoms[atoms] and temp_atom[4] in chains:
                    atom_vertex = [float(temp_atom[6]), float(temp_atom[7]), float(temp_atom[8])]
                    verticies.append(atom_vertex)
    return array(verticies)

def calculate_rmsd(protein1, protein2):
	return sqrt(mean(npsum((protein1 - protein2)**2, axis=1)))


def svd_alignment(protein1, protein2):
    protein1 = array(protein1)
    protein2 = array(protein2)
    
    center1 = mean(protein1, axis=0)
    center2 = mean(protein2, axis=0)
    protein1_centered = protein1 - center1
    protein2_centered = protein2 - center2
    
    U, _, Vt = svd(protein2_centered.T @ protein1_centered)
    rotation_matrix = Vt.T @ U.T
    
    protein2_aligned = protein2_centered @ rotation_matrix
    
    return protein1_centered, protein2_aligned


def lsf_alignment(protein1, protein2):
    protein1 = array(protein1)
    protein2 = array(protein2)

    center1 = mean(protein1, axis=0)
    center2 = mean(protein2, axis=0)
    protein1_centered = protein1 - center1
    protein2_centered = protein2 - center2

    rotation_matrix, _ = orthogonal_procrustes(protein2_centered, protein1_centered)
    
    protein2_aligned = protein2_centered @ rotation_matrix
    
    return protein1_centered, protein2_aligned


def quaternion_alignment(protein1, protein2):
    protein1 = array(protein1)
    protein2 = array(protein2)

    center1 = mean(protein1, axis=0)
    center2 = mean(protein2, axis=0)
    protein1_centered = protein1 - center1
    protein2_centered = protein2 - center2

    covariance_matrix = protein2_centered.T @ protein1_centered

    rotation_quaternion = Rotation.from_matrix(covariance_matrix).as_quat()
    rotation_matrix = Rotation.from_quat(rotation_quaternion).as_matrix()

    protein2_aligned = protein2_centered @ rotation_matrix

    return protein1_centered, protein2_aligned


def kabsch_alignment(protein1, protein2):
    protein1 = array(protein1)
    protein2 = array(protein2)

    center1 = mean(protein1, axis=0)
    center2 = mean(protein2, axis=0)
    protein1_centered = protein1 - center1
    protein2_centered = protein2 - center2

    covariance_matrix = protein2_centered.T @ protein1_centered

    U, _, Vt = svd(covariance_matrix)

    optimal_rotation_matrix = U @ Vt

    protein2_aligned = protein2_centered @ optimal_rotation_matrix

    return protein1_centered, protein2_aligned


reference_verticies = extract_vertices(reference_pdb, atoms, reference_chains)
test_verticies = extract_vertices(test_pdb, atoms, test_chains)

reference_vertex_count = reference_verticies.shape[0]
test_vertex_count = test_verticies.shape[0]

if reference_vertex_count != test_vertex_count:
	print(f"Number of verticies in reference ({reference_vertex_count}) is not equal to"
		  f"number of verticies in test ({test_vertex_count}).")

if alignment_method == 'kabsch':
	p1, p2 = kabsch_alignment(reference_verticies, test_verticies)
elif alignment_method == 'lsf':
	p1, p2 = lsf_alignment(reference_verticies, test_verticies)
elif alignment_method == 'svd':
	p1, p2 = svd_alignment(reference_verticies, test_verticies)
elif alignment_method == 'quaternion':
	p1, p2 = quaternion_alignment(reference_verticies, test_verticies)


rmsd = calculate_rmsd(p1, p2)
end_time = time.time()
elapsed_time = end_time - start_time

print(f"RMSD = {rmsd}, Time = {elapsed_time:.2f} seconds")