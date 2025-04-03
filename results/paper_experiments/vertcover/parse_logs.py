import os, sys

sys.path.append('../../../src')
from plot_scan import parse_optimised_param_results

path = "./optimised"

if __name__ == "__main__":

     for dir_name in [path + "/" + dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
      if dir_name.endswith("symmetry"):
         graphseed = dir_name[17:23]
         qubits = dir_name[34:36]
         logpath = f"{path}/log{qubits}qubits{graphseed}.txt"
         parse_optimised_param_results(dir_name, logpath, 7)