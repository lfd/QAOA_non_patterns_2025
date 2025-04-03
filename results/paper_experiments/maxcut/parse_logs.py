import os, sys

sys.path.append('../../../src')
from plot_scan import parse_optimised_param_results

path = "./optimised"

if __name__ == "__main__":

     for dir_name in [path + "/" + dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
      if dir_name.endswith("symmetry"):
         qubits = dir_name[22:24]
         logpath = f"{path}/logopt{qubits}qubits.txt"
         parse_optimised_param_results(dir_name, logpath, 7)