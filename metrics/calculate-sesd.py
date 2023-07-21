import argparse
import numpy as np
import sys
import time

start_time = time.time()

parser = argparse.ArgumentParser(
                    prog='Calculate Smallest Enclosing Sphere Diameter',
                    description='Takes in a PDB file and estimates how compact a protein is by '
                                'approximating the diameter of a sphere that encloses all atoms '
                                'of the protein.')

parser.add_argument('-i', '--input', dest='pdb_file', required=True, help='Input PDB file.')
parser.add_argument('-c', '--chain', dest='chain', required=True, help='Which chain to consider.')
parser.add_argument('-n', '--num-samples', dest='num_samples', type=int, required=True,
                    help='Number of samples to average.')
parser.add_argument('-z', '--z-threshold', dest='z_threshold', required=True, 
                    help='Z score threshold for rejecting outliers.')
parser.add_argument('-a', '--atoms', dest='atoms', required=True,
                    choices=["all", "backbone", "backbone-no-carbonyl", "alpha-carbons"],
                    help='Atoms to consider.')
parser.add_argument('-o', '--outfile', dest='output_file', default='sesd.txt',
                    help='Output file name.')

args = parser.parse_args()

pdb_file = args.pdb_file
num_samples = args.num_samples
atoms = args.atoms
z_threshold = args.z_threshold
chain = args.chain

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
    return np.array(verticies)


def welzl(points):
    def point_inside_sphere(point, center, radius):
        distance = np.linalg.norm(point - center)
        return distance < radius

    def find_sphere_center(p1, p2, p3):
        A = np.vstack((p1, p2, p3))
        b = np.sum(np.square(A), axis=1)
        center = 0.5 * np.linalg.solve(A, b)
        return center

    def point_inside_triangle(point, p1, p2, p3):
        v1 = p3 - p1
        v2 = p2 - p1
        v3 = point - p1

        dot11 = np.dot(v1, v1)
        dot22 = np.dot(v2, v2)
        dot12 = np.dot(v1, v2)
        dot13 = np.dot(v1, v3)
        dot23 = np.dot(v2, v3)

        inv_denom = 1.0 / (dot11 * dot22 - dot12 * dot12)
        u = (dot22 * dot13 - dot12 * dot23) * inv_denom
        v = (dot11 * dot23 - dot12 * dot13) * inv_denom

        return np.logical_and.reduce([0 <= u, u <= 1, 0 <= v, v <= 1, u + v <= 1])


    def sphere_from_points(points):
        center = np.mean(points, axis=0)
        radius = np.max(np.linalg.norm(points - center, axis=1))
        return center, radius

    def sphere_from_4_points(boundary_points):
        center, radius = sphere_from_points(boundary_points[:3])

        for i in range(3, len(boundary_points)):
            point = boundary_points[i]
            if not point_inside_triangle(point, *boundary_points[:3]):
                center = find_sphere_center(*boundary_points[:3], point)
                radius = np.linalg.norm(point - center)

        return center, radius

    def welzl_helper(points, boundary_points):
        if len(boundary_points) == 4:
            return sphere_from_4_points(boundary_points)

        if len(points) == 0 or len(boundary_points) == 3:
            if len(boundary_points) == 0:
                return np.zeros(3), 0.0
            elif len(boundary_points) == 1:
                return boundary_points[0], 0.0
            elif len(boundary_points) == 2:
                center = np.mean(boundary_points, axis=0)
                radius = np.linalg.norm(boundary_points[0] - center)
                return center, radius
            else:
                center = np.mean(boundary_points, axis=0)
                radius = np.linalg.norm(boundary_points[0] - center)
                return center, radius

        p = points[0]
        center, radius = welzl_helper(points[1:], boundary_points)

        if point_inside_sphere(p, center, radius):
            return center, radius

        return welzl_helper(points[1:], boundary_points + [p])

    points = np.array(points)
    np.random.shuffle(points)  # Shuffle points for better performance
    center, radius = welzl_helper(points, [])
    return center, radius

def calculate_sesd(pep_verts):
    center, radius = welzl(pep_verts)
    sesd = radius * 2
    return center, radius, sesd


def reject_outliers(data, threshold=3):
    """
    Reject outliers from a list of floats using the Z-score method with median.
    
    Args:
        data (list): A list of float values.
        threshold (float): The threshold for considering a value as an outlier.
                           Default is 3.
                           
    Returns:
        list: A new list with outliers removed.
    """
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    z_scores = 0.6745 * (data - median) / mad
    filtered_data = [x for x, z in zip(data, z_scores) if np.all(np.abs(z) < threshold)]
    
    return filtered_data


def calculate_average_sesd(pdb, n, a, z_threshold):
    verticies = extract_vertices(pdb, a)
    sys.setrecursionlimit(10**6)

    center_list = []
    radius_list = []
    sesd_list = []

    for _ in range(n):
        try:
            center, radius, sesd = calculate_sesd(verticies)
            center_list.append(center)
            radius_list.append(radius)
            sesd_list.append(sesd)
        except:
            print("Max recursive depth of 10^6 was reached. Try with less atoms. setting value to -1.")
            center_list.append(999)
            radius_list.append(999)
            sesd_list.append(999)

    filtered_sesd = reject_outliers(sesd_list, float(z_threshold))

    average_center = np.mean(center_list, axis=0)
    average_radius = np.mean(radius_list)
    average_sesd = np.mean(filtered_sesd)

    sesd_standard_deviation = np.std(filtered_sesd)

    final_n = len(filtered_sesd)

    return average_center, average_radius, average_sesd, sesd_standard_deviation, final_n


center, radius, sesd, sesd_standard_deviation, final_n = calculate_average_sesd(pdb_file, num_samples, atoms, z_threshold)
end_time = time.time()
elapsed_time = end_time - start_time

print(f"SESD = {sesd}, Standard Deviation = {sesd_standard_deviation} ({final_n}/{num_samples}), Time = {elapsed_time:.2f} seconds")