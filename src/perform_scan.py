
from problem_util import load_problem
from optimize_from_scan import read_best_existing_params

import os
import json
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from qiskit.primitives import Estimator
from qiskit.circuit.library import QAOAAnsatz

from qiskit_algorithms.utils import algorithm_globals

def make_scan_params(beta_bounds,gamma_bounds,resolution,fixed_params=None):
    print(f"creating cartesian grid over {beta_bounds}x{gamma_bounds} with resolution {resolution}x{resolution} and fixed parameters {fixed_params}...")
    params = []
    if fixed_params!=None:
        fixed_beta = fixed_params[:len(fixed_params)//2]
        fixed_gamma = fixed_params[len(fixed_params)//2:]
    beta_interval = np.flip(np.linspace(beta_bounds[0],beta_bounds[1],resolution))
    gamma_interval = np.flip(np.linspace(gamma_bounds[0],gamma_bounds[1],resolution))
    for i in beta_interval:
        for j in gamma_interval:
                if fixed_params!=None:
                    params.append(fixed_beta + [i] + fixed_gamma + [j])
                else:
                    params.append([i,j])
    return params

def get_best_params(res, scan, p):
    print("evaluating best parameters...")
    #evaluate best params and save in dict
    min_results = res.values.argsort()[:10]
    best_params = {"params":[],"expectation":[]}
    print(f"Best params for depth {p}:")
    for i in min_results:
        best_params["params"].append(scan[i])
        best_params["expectation"].append(res.values[i])
        print(f'params {scan[i]} with expectation {res.values[i]}')
    return best_params

def scan_qaoa_landscape_by_depth(estimator, pauli_op, depth, n, beta_bounds, gamma_bounds, resolution, outpath, shots=None, fixed_parameters=None):
    #create QAOA ansatz for current depth
    print(f"creating QAOAAnsatz of depth {depth}...")
    ansatz = QAOAAnsatz(cost_operator=pauli_op, reps=depth)

    # create scan parameters for current depth using fixed parameters for previous depths
    print("creating scan params...")
    scan_params = make_scan_params(beta_bounds=beta_bounds,gamma_bounds=gamma_bounds,resolution=resolution,fixed_params=fixed_parameters)
    #print(f"scan parameters: {scan_params}")
    # create list of ansätze and observables of same length as scan parameters
    circuits = [ansatz]*len(scan_params)
    observables = [pauli_op]*len(scan_params)

    print("running estimator with scan parameters...")
    # run estimator for ansatz with the different parameters in scan params
    result = estimator.run(circuits,observables,scan_params,shots=shots).result()

    #evaluate best params and save in dict
    best_params = get_best_params(result,scan_params,depth)
    fixed_params = best_params["params"][0]
    
    #prepare data for plotting and json dump
    scan_grid = np.array(scan_params).reshape(-1,resolution,len(fixed_params))
    result_grid = result.values.reshape(-1,resolution)

    #save data as json
    print("saving results as json...")
    data = {"scan_grid":scan_grid.tolist(),"expecation_values":result_grid.tolist(),"best":best_params,"metadata":result.metadata}
    with open(f"{outpath}/estimator_data_{resolution}x{resolution}scan_p{depth}_{n}qubits.json", 'w') as f:
       json.dump(data, f)

    #return fixed parameters
    return fixed_params

def run_scan_with_iterative_parameter_fixing(problemtype,problem_kwargs,methodname,problem_path,n,target_p,resolution=3,beta_bounds =[-np.pi,np.pi],gamma_bounds =[-np.pi,np.pi],fixed_params_list=None,shots=None,seed=None,file_affix="",pstep=1,pstart=1,resume=None,resume_startparams=False):
    '''
    problem_path: path to the problem instance
    n: problem size
    target_p: qaoa depth up to which a layer-by-layer based scan is performed
    resolution: grid resolution (results in grid with resolution * resolution points and step size len(bounds)/resolution)
    beta_bounds: beta grid bounds
    gamma_bounds: gamma grid bounds
    shots: number of shots for the estimator (none = calculate exact expectation value)
    seed: seed for statevector simulator
    fixed_params_list: optional list of parameters to which lower p params are fixed to when searching for params of QAOA with depth>1
    file_affix: string to add to end of filename
    '''
    p = pstart
    if resume!=None:
        if resume_startparams==True:
            fixed_beta = fixed_params_list[:p-1]
            fixed_gamma = fixed_params_list[target_p:target_p+p-1]
            fixed_params = np.concatenate([fixed_beta,fixed_gamma]).tolist()
        else:
            fixed_params = read_best_existing_params(resume)
    else:
        fixed_params = None
    
    algorithm_globals.random_seed = seed

    print("-----------------------------------------------------------------------------------------------------------------------------------------------")
    print(f"Running qaoa with {n} qubits from p={p} to p={target_p} within bounds {beta_bounds}x{gamma_bounds} with resolution {resolution}x{resolution}")
    print("-----------------------------------------------------------------------------------------------------------------------------------------------")

    print("creating result directory...")
    index = 1
    if not seed:
        seeddata = ""
    else:
        seeddata = f"seed{seed}_"
    if problemtype=="maxcut" or problemtype=="vertcover" or problemtype=="MIS":
        result_path = f"./results/reproduction/{problemtype}/{methodname}/graph{problem_kwargs['seed']}_{resolution}x{resolution}scan_{n}qubits_p{p}to{target_p}_{seeddata}{file_affix}"
    elif problemtype=="max3sat" and problem_kwargs is not None:
        if problem_kwargs["alpha"]>4.9 or problem_kwargs["alpha"]<=3.5:
            result_path = f"./results/reproduction/{problemtype}/easy/{methodname}/{resolution}x{resolution}scan_{n}qubits_p{p}to{target_p}_{seeddata}{file_affix}"
        else:
            result_path = f"./results/reproduction/{problemtype}/hard/{methodname}/{resolution}x{resolution}scan_{n}qubits_p{p}to{target_p}_{seeddata}{file_affix}"
    else:
        result_path = f"./results/reproduction/{problemtype}/{methodname}/{resolution}x{resolution}scan_{n}qubits_p{p}to{target_p}_{seeddata}{file_affix}"
    os.makedirs(result_path, exist_ok=True)

    print("loading problem...")
    pauli_op, offset = load_problem(problemtype,n,problem_path)

    print(f"Problem: {pauli_op}, offset:{offset}")
    offset_data={"offset": offset}
    with open(f"{result_path}/offset.json", 'w') as f:
       json.dump(offset_data, f)

    print("instantiating estimator...")
    estimator = Estimator()

    print("beginning scan...")
    while p <= target_p:
        print("  ")
        print(f"-------Depth {p}-------")
        if fixed_params_list is None:
            fixed_params = scan_qaoa_landscape_by_depth(estimator=estimator,pauli_op=pauli_op,depth=p,n=n,beta_bounds=beta_bounds,gamma_bounds=gamma_bounds,resolution=resolution,outpath=result_path,shots=shots,fixed_parameters=fixed_params)
        else:
            fixed_beta = fixed_params_list[:p-1]
            fixed_gamma = fixed_params_list[target_p:target_p+p-1]
            fixed_params = np.concatenate([fixed_beta,fixed_gamma]).tolist()
            if p > 1:
                print(f"fixing parameters to imported parameters {fixed_params} for beta{p},gamma{p} and increasing depth...")
            search_result_params = scan_qaoa_landscape_by_depth(estimator=estimator,pauli_op=pauli_op,depth=p,n=n,beta_bounds=beta_bounds,gamma_bounds=gamma_bounds,resolution=resolution,outpath=result_path,shots=shots,fixed_parameters=fixed_params)
        p += pstep
        

    print("-----------------------------------------------------------------------------------------------------------------------------------------------")
    print("Scan finished")

    return result_path
    
    
    

