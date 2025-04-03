import os, sys
import numpy as np
import json

sys.path.append('../../src')

from plot_scan import read_existing_result_data

problems = ["maxcut","max3sat/hard","max3sat/easy","vertcover","higher_depth_maxcut"]
methods = ["optimised","linear_ramp","linear_ramp_neg_mixer","sequential","fixed"]

def save_only_if_higher(save_val,problem,method,qubits,seed):
   save_path = f"{problem}/max_energies/n{qubits}_seed{seed}.json"
   if os.path.isfile(save_path)==True:#
      print("exists")
      with open(save_path) as f:
          data = json.load(f)
          prev_max = data["max_energy"]
      if save_val>prev_max:
         save(save_path,save_val,method)
   else:
      save(save_path,save_val,method)
   
def save(path,val,method):
   save_data = {"max_energy":val,"found_by":method}
   with open(path, 'w') as f:
      json.dump(save_data, f)

def read_and_update_max_vals(dir_name,out_name,meth,qubits,seed):
   mat, scan, p, resolution, best_params = read_existing_result_data(dir_name,add_offset=False)
   max = np.amax(mat)

   save_only_if_higher(max,out_name,meth,qubits,seed)
   
if __name__ == "__main__":
     
   for prob in problems:
      for meth in methods:
         print("----")
         if meth=="fixed" and prob!="maxcut" and prob!="higher_depth_maxcut":
            continue
         if prob=="higher_depth_maxcut" and meth!="linear_ramp" and meth!="optimised":
            continue
         path = f"./{prob}/{meth}"

         for dir_name in [path + "/"+ dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
            if dir_name.endswith("symmetry"):
               if prob=="maxcut" or prob=="vertcover" or prob=="higher_depth_maxcut":
                  qubits = dir_name[len(path)+23:len(path)+25]
                  seed = dir_name[len(path)+6:len(path)+12]
               else:
                  qubits = dir_name[len(path)+11:len(path)+13]
                  seed = 123456

               mat, scan, p, resolution, best_params = read_existing_result_data(dir_name,add_offset=False)
               max = np.amax(mat)

               read_and_update_max_vals(max,prob,meth,qubits,seed)
            
            

            