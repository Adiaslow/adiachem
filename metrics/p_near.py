import argparse
import numpy as np
import pandas as pd
import time

start_time = time.time()

parser = argparse.ArgumentParser(
                    prog='Calculate P-Near Score.',
                    description='Given a score file containing RMSD and Interface Score.')
parser.add_argument('-i', '--input', dest='score_file', required=True,
                    help='Input score file.')
parser.add_argument('-l', '--lambda', dest='lamb', type=float, default=1.0,
                    help='Lambda value.')
parser.add_argument('-b', '--boltzman', dest='k_B', type=float, default=1.0,
                    help='Boltzman k value.')
parser.add_argument('-t', '--tempurature', dest='T', type=float, default=1.0,
                    help='Tempurature value.')

args = parser.parse_args()
score_file = args.score_file
lamb = args.lamb
k_B = args.k_B
T = args.T


def parse_score_file(file):
  f = open(file, "r")
  data = f.read()
  f.close()

  data = [[element for element in row.split("\t")[1:] if element != ""] for row in data.split('\n')]
  df = pd.DataFrame(data[1:-1],columns=data[0])

  rmsd = [float(element) for element in df["RMSD"]]
  interface_score = [float(element) for element in df["Interface score"]]
  return rmsd, interface_score


def calculate_p_near(rmsd, interface_score, lamb, k_B = 1, T = 1) -> float:
  numerator_sum = 0
  denominator_sum = 0
  for i in range(0, len(rmsd)):
    numerator_sum += np.exp(-rmsd[i]**2 / lamb**2) * np.exp(-interface_score[i] / k_B * T)
    denominator_sum += np.exp(-interface_score[i] / k_B * T)
  return numerator_sum / denominator_sum

rmsd, interface_score = parse_score_file(score_file)
p_near = calculate_p_near(rmsd, interface_score, lamb, k_B, T)

end_time = time.time()
elapsed_time = end_time - start_time

print(f"P-Near = {p_near}, Time = {elapsed_time:.2f} seconds")