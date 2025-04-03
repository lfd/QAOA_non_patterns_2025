from problem_util import load_problem

import json,os

from qiskit.primitives import Sampler
from qiskit_algorithms import QAOA
from qiskit_algorithms.utils import algorithm_globals

from functools import partial
from scipy.optimize import minimize

def read_best_existing_params(path_to_json, single_file=False):
    if not single_file:
        json_res = []
        for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
            with open(path_to_json + "/" + file_name) as json_file:
                if file_name!="offset.json" and file_name!="eigenvalues.json" and file_name!="opt_data.json":
                    data = json.load(json_file)
                    json_res.append(data["best"]["params"])
        init_params = json_res[-1][0]
    else:
        with open(path_to_json) as json_file:
                data = json.load(json_file)
                init_params = data['opt_params']
            
    return init_params

def optimize_from_inits(problemtype, n, opt_method, problem_path, param_path, seed=None):
    if seed:
        print(f"setting algorithm globals seed to {seed}...")
        algorithm_globals.random_seed = seed
    
    print("loading problem...")
    pauli_op, offset = load_problem(problemtype,n,problem_path)
    print(f"Problem: {pauli_op}, offset: {offset}")

    print("instantiating sampler...")
    sampler = Sampler()

    print(f"instantiating {opt_method} optimizer...")
    optimizer = partial(minimize,method=opt_method,options={'disp': True})

    print(f"loading best parameters for highest p found in {problem_path}...")
    if param_path.endswith('.json'):
        initial_params = read_best_existing_params(param_path,True)
    else:
        initial_params = read_best_existing_params(param_path)
    p = len(initial_params)//2

    
    print("-----------------------------------------------------------------------------------------------------------------------------------------------")
    print(f"Optimizing QAOA parameters for {p} layers using {opt_method} starting from initial point {initial_params}...")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------")

    print(f"instantiating QAOA wrapper...")
    qaoa = QAOA(sampler, optimizer, reps=p, initial_point = initial_params)

    print("running optimization...")
    result = qaoa.compute_minimum_eigenvalue(pauli_op)
    result_dict = {"n": n,
                    "p": p,
                    "init_params": initial_params,
                    "optimizer"
                    "optimizer_evals": result.cost_function_evals,
                    "opt_params": result.optimal_point,
                    "opt_value": result.optimal_value,
                    "optimizer_time":result.optimizer_time}
    
    print("Optimization finished")
    for key in result_dict:
       print(f"{key}: {result_dict[key]}")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------")

    return result_dict