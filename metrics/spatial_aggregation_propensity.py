import argparse
import numpy as np
import os
import sys


path = os.path.dirname(os.path.realpath("../atlas.py"))
sys.path.append(path)
import atlas

from metrics.sasa import get_sasa

class Atom:
	def __init__(self, coordinates, hydrophobicity, sasa):
		self.coordinates = coordinates
		self.hydrophobicity = hydrophobicity
		self.sasa = sasa

def map_range(value, x_min, x_max, y_min, y_max):
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    # Avoid division by zero
    if x_range == 0:
        return y_min if y_range > 0 else y_max
    
    # Map the value to the new range
    mapped_value = ((value - x_min) / x_range) * y_range + y_min
    
    # Ensure the mapped value is within the bounds of the new range
    return min(max(mapped_value, y_min), y_max)

def get_sap(scale, sap_radius, sasa_radius, iterations):
	total_accessible_area, solvent_accessible_areas, df = get_sasa(pdb, sasa_radius)
	coordinates = df[['X', 'Y', 'Z']].values

	hydrophobicities = [atlas.hydrophobicity_scale[scale][h] for h in df['ResidueName'].to_numpy()]
	hydro_min = np.min(hydrophobicities)
	hydro_max = np.max(hydrophobicities)
	normed_hydrophobicities = [map_range(hydrophobicity, hydro_min, hydro_max, -0.5, 0.5) for hydrophobicity in hydrophobicities]

	sasa_of_fully_exposed_res = [atlas.sasa_of_fully_exposed_residues[h] for h in df['ResidueName'].to_numpy()]

	atoms = [Atom(coordinates[i], hydrophobicities[i], solvent_accessible_areas[i]) for i in range(0, len(coordinates))]
	
	return "Hi" # sap, saps


if __name__ == "__main__":
	parser = argparse.ArgumentParser(prog='Calculates Solvent Accessible Surface Area (SASA).')
	parser.add_argument('-i', '--input', dest='pdb', required=True)
	parser.add_argument('-r0', '--sap-radius', dest='sap_radius', default=5.0, type=float)
	parser.add_argument('-r1', '--sasa-radius', dest='sasa_radius', default=1.4, type=float)
	parser.add_argument('-c', '--hydrophobicity-scale', dest='hydrophobicity_scale_choice', default='Black-Mould',
						choices=['Bandyopadhyay-Mehler', 'Black-Mould', 'Eisenberg', 'Kyte-Doolittle', 
			 					 'Meek', 'Rose', 'Wimley-White', 'Jain', 'Miyazawa'], type=str)
	parser.add_argument('-t', '--iterations', dest='iterations', type=int)
	parser.add_argument('-v', '--verbose', dest='verbose', default=False, type=bool)

	args = parser.parse_args()

	pdb = args.pdb
	sap_radius = args.sap_radius
	sasa_radius = args.sasa_radius
	hydrophobicity_scale_choice = args.hydrophobicity_scale_choice
	iterations = args.iterations
	verbose = args.verbose

	sap, saps = get_sap(hydrophobicity_scale_choice, sap_radius, sasa_radius, iterations)

	print(f"Total Spatial Agragation Propensity: {sap}")
	if verbose:
	    print(saps)
