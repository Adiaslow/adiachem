import argparse
from numpy import sqrt, mean, array, empty
from numpy import sum as npsum
from numpy.linalg import norm
import time

start_time = time.time()

parser = argparse.ArgumentParser(
                    prog='Calculate P-Near Score.',
                    description='Given a score file containing RMSD and Interface Score.')
parser.add_argument('-i', '--input', dest='input_pdb', required=True)
parser.add_argument('-c', '--centroid-atoms', dest='centroid_atoms', type=str)
parser.add_argument('-b', '--binder', dest='binder', type=str)
args = parser.parse_args()

def extract_atoms_coordinates(pdb_file_path, centroid_atoms_string):
    centroid_atoms = []
    centroid_atoms_list = centroid_atoms_string.split(":")
    for atom_info in centroid_atoms_list:
        # print(atom_info.split(","))
        chain_id, residue_serial, atom_name = atom_info.split(",")
        residue_serial = int(residue_serial)
        centroid_atoms.append((chain_id, residue_serial, atom_name))

    vertices = []
    with open(pdb_file_path, "r") as pdb_file:
        for line in pdb_file:
            record_type = line[0:6].strip()

            if record_type == "ATOM":
                atom_serial_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21].strip()
                residue_serial_number = int(line[22:26].strip())
                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())

                for centroid_atom in centroid_atoms:
                    if (chain_id, residue_serial_number, atom_name) == centroid_atom:
                        vertices.append([x_coord, y_coord, z_coord])
                        break

    if len(vertices) == len(centroid_atoms):
        vertices = array(vertices)
        centroid = mean(vertices, axis=0)
        # print("Centroid:", centroid)
        return centroid
    else:
        print("Not all atoms found.")
        return empty()

def extract_chain_vertices(pdb_file_path, target_chain_id):
    chain_atoms = []

    with open(pdb_file_path, "r") as pdb_file:
        for line in pdb_file:
            record_type = line[0:6].strip()

            if record_type == "ATOM":
                atom_serial_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21].strip()
                residue_serial_number = int(line[22:26].strip())
                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())

                if chain_id == target_chain_id and atom_name == "CA":
                    chain_atoms.append([x_coord, y_coord, z_coord])

    if chain_atoms:
        # print(f"Vertices for chain {target_chain_id}:")
        # print(chain_atoms)
        return(chain_atoms)
    else:
        print(f"No atoms found for chain {target_chain_id}.")
        return empty()

def calculate_distances(vertices_array, given_vertex):
    vertices_array = array(vertices_array)
    given_vertex = array(given_vertex)

    # Calculate the distance between the given vertex and all vertices in the array
    distances = norm(vertices_array - given_vertex, axis=1)

    return distances

centroid = extract_atoms_coordinates(args.input_pdb, args.centroid_atoms)
binder_vertices = extract_chain_vertices(args.input_pdb, args.binder)
distances = calculate_distances(centroid, binder_vertices)
rmsd = sqrt(mean(npsum((centroid - binder_vertices)**2, axis=1)))
n_term_score = round(distances[0]/rmsd,3)
c_term_score = round(distances[-1]/rmsd,3)

end_time = time.time()
elapsed_time = end_time - start_time

print(f"N Term Score = {n_term_score}, C Term Score = {c_term_score}")
