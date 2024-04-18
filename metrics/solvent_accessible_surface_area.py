import argparse
import numpy as np
from scipy.spatial import cKDTree
import os
import sys

path = os.path.dirname(os.path.realpath("../atlas.py"))
sys.path.append(path)
import atlas

path = os.path.dirname(os.path.realpath("../utilities/pdb_parser.py"))
sys.path.append(path)
from pdb_parser import parse_pdb_to_dataframe, dataframe_to_pdb

atomic_radii = atlas.atomic_radii

def rolling_ball_sasa(coordinates, radii, probe_radius):
    if len(coordinates) != len(radii):
        raise ValueError("Number of coordinates and radii should be the same.")

    num_atoms = len(coordinates)

    # Convert radii list to numpy array
    radii = np.array(radii)

    atom_list = np.array(coordinates)
    total_sasa = 0.0
    sasas = []

    for i in range(num_atoms):
        atom_i = atom_list[i]
        radius_i = radii[i] + probe_radius  # Effective radius considering probe size

        tree = cKDTree(atom_list)

        # Find atoms within the sum of radii + probe radius distance from atom_i
        neighbors = tree.query_ball_point(atom_i, radius_i)

        sasa_i = 0.0
        for neighbor in neighbors:
            if neighbor != i:
                atom_j = atom_list[neighbor]
                distance = np.linalg.norm(atom_i - atom_j)

                if distance < radius_i:
                    # Calculate exposed area between atom_i and atom_j
                    d_i = radii[i]  # Atom_i's radius
                    d_j = radii[neighbor]  # Atom_j's radius

                    h = (d_i + d_j + probe_radius) / 2 - distance
                    exposed_area = 2 * np.pi * h * (radius_i - h)

                    sasa_i += exposed_area
        if sasa_i < 0:
            sasa_i = 0
        sasas.append(sasa_i)
        total_sasa += sasa_i

    return total_sasa, sasas

def get_sasa(pdb, solvent_probe_radius):
    df = parse_pdb_to_dataframe(pdb)
    coordinates = df[['X', 'Y', 'Z']].values
    atom_types = df['AtomType'].values
    radii = [atlas.atomic_radii[atom_type] for atom_type in atom_types]

    total_accessible_area, solvent_accessible_areas = rolling_ball_sasa(coordinates, radii, solvent_probe_radius)

    df['SASA'] = solvent_accessible_areas
    return total_accessible_area, solvent_accessible_areas, df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Calculates Solvent Accessible Surface Area (SASA).')
    parser.add_argument('-i', '--input', dest='pdb', required=True)
    parser.add_argument('-r', '--radius', dest='solvent_probe_radius', required=True, type=float)
    parser.add_argument('-v', '--verbose', dest='verbose', default=False, type=bool)
    args = parser.parse_args()
    pdb = args.pdb
    solvent_probe_radius = args.solvent_probe_radius
    verbose = args.verbose

    total_accessible_area, solvent_accessible_areas, df = get_sasa(pdb, solvent_probe_radius)
    dataframe_to_pdb(df, "sasa.pdb")
    print(f"Total Accessible Surface Area: {total_accessible_area}")
    if verbose:
        print(solvent_accessible_areas)
