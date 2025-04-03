import os, sys

sys.path.append('../../src')
from problem_util import create_problem, load_problem, get_min_energy_eigenvalue, load_offset

graphsizes = [10, 12, 14, 16]
graphseeds = {10: [123456, 123457, 123459, 123461, 123464, 123465, 123466, 123468, 123474, 123475], 
              12: [123456, 123457, 123458, 123459, 123460, 123462, 123463, 123464, 123465, 123466], 
              14: [123456, 123457, 123458, 123459, 123460, 123461, 123462, 123463, 123464, 123466], 
              16: [123456, 123457, 123458, 123459, 123460, 123461, 123462, 123463, 123464, 123465]}
max3satsizes = [14, 15, 16, 17, 19, 20, 21, 22, 23, 24]
alphas = {"easy":[2.5, 2.75, 3, 3.25, 5.3333, 5.6666, 6, 6.3333, 6.6666, 5], "hard":[3.6666, 4, 4.3333, 4.6666, 3.75, 4, 4.25, 4.5, 4.75, 3.8]}

def execute(outpath,problem,n,degree,alpha,seed):
    problem_path = create_problem(problem,n,kwargs={"degree":degree,"alpha":alpha,"seed":seed})
    pauli_op, offset = load_problem(problem,n,problem_path)
    get_min_energy_eigenvalue(pauli_op,outpath,fname=f"min_eigenvalue_{seed}_{n}qubits")

def get_offset(outpath,problem,n,degree,alpha,seed):
    problem_path = create_problem(problem,n,kwargs={"degree":degree,"alpha":alpha,"seed":seed})
    load_offset(problem,problem_path,n,outpath)


def solve(problem="maxcut",add_info=None):
    if problem!="max3sat":
        outpath = f'./{problem}/eigensolver'
        for n in graphsizes:
            seeds = graphseeds[n]
            for seed in seeds:
                execute(outpath,problem,n,3,None,seed)
    else:
        outpath = f'./{problem}/{add_info}/eigensolver'
        for n, alpha in zip(max3satsizes, alphas[add_info]):
                execute(outpath,problem,n,None,alpha,123456)

def offsets(problem="maxcut",add_info=None):
    if problem!="max3sat":
        for n in graphsizes:
            seeds = graphseeds[n]
            for seed in seeds:
                if problem=="maxcut":
                    outpaths = [f"./maxcut/linear_ramp/graph{seed}_32x32scan_{n}qubits_p1to7_None_optimized_maxcut_symmetry",
                                f"./maxcut/linear_ramp_neg_mixer/graph{seed}_32x32scan_{n}qubits_p1to7_None_optimized_maxcut_symmetry",
                                f"./maxcut/sequential/graph{seed}_32x32scan_{n}qubits_p1to7_maxcut_symmetry",
                                f"./maxcut/fixed/graph{seed}_32x32scan_{n}qubits_p1to7_None_optimized_maxcut_symmetry",
                                f"./maxcut/optimised/graph{seed}_32x32scan_{n}qubits_p1to7_COBYLA_optimized_maxcut_symmetry",
                                f"./higher_depth_maxcut/linear_ramp/graph{seed}_32x32scan_{n}qubits_p1to21_None_optimized_maxcut_symmetry",
                                f"./higher_depth_maxcut/linear_ramp_neg_mixer/graph{seed}_32x32scan_{n}qubits_p1to21_None_optimized_maxcut_symmetry",
                                f"./higher_depth_maxcut/optimised/graph{seed}_32x32scan_{n}qubits_p1to21_COBYLA_optimized_maxcut_symmetry",
                                f"./higher_depth_maxcut/sequential/graph{seed}_32x32scan_{n}qubits_p1to21_maxcut_symmetry"]
                elif problem=="vertcover":
                     outpaths = [f"./vertcover/linear_ramp/graph{seed}_32x32scan_{n}qubits_p1to7_None_optimized_general_symmetry",
                                 f"./vertcover/linear_ramp_neg_mixer/graph{seed}_32x32scan_{n}qubits_p1to7_None_optimized_general_symmetry",
                                f"./vertcover/sequential/graph{seed}_32x32scan_{n}qubits_p1to7_general_symmetry",
                                f"./vertcover/optimised/graph{seed}_32x32scan_{n}qubits_p1to7_COBYLA_optimized_general_symmetry"]
                     
                get_offset(outpaths,problem,n,3,None,seed)
    else:
        for n, alpha in zip(max3satsizes, alphas[add_info]):
                outpaths = [f"./{problem}/{add_info}/sequential/32x32scan_{n}qubits_p1to7_general_symmetry",
                            f"./{problem}/{add_info}/linear_ramp/32x32scan_{n}qubits_p1to7_None_optimized_general_symmetry",
                            f"./{problem}/{add_info}/linear_ramp_neg_mixer/32x32scan_{n}qubits_p1to7_None_optimized_general_symmetry",
                            f"./{problem}/{add_info}/optimised/32x32scan_{n}qubits_p1to7_COBYLA_optimized_general_symmetry"]
                get_offset(outpaths,problem,n,None,alpha,123456)

def make_config_combinations(problem="maxcut"):
     counter=0
     if problem!="max3sat":
        for n in graphsizes:
            seeds = graphseeds[n]
            for seed in seeds:
                print(f"[instance{counter}]")
                print(f"qubits={n}")
                print(f"seed={seed}")
                print("degree=3")
                print("alpha=false\n")
                counter+=1
     else:
        for add_info in ["easy","hard"]:
            for n, alpha in zip(max3satsizes, alphas[add_info]):
                print(f"[instance{counter}]")
                print(f"qubits={n}")
                print("seed=123456")
                print("degree=false")
                print(f"alpha={alpha}\n")
                counter+=1


if __name__ == "__main__":
    #runs in ~20mins on my device
    #solve("maxcut")
    #solve("vertcover")
    #solve("max3sat","easy")
    #solve("max3sat","hard")


    #offsets("maxcut")
    #offsets("vertcover")
    #offsets("max3sat","hard")
    #offsets("max3sat","easy")
    #print("edit comments in main to run")

    make_config_combinations("maxcut")
    make_config_combinations("max3sat")