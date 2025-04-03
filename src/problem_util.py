import numpy as np
import networkx as nx
import json
import itertools
import os
from qiskit_optimization.applications import Maxcut, VertexCover, StableSet
from qiskit_optimization.problems import QuadraticProgram
from qiskit_optimization.converters import QuadraticProgramToQubo
from qiskit_algorithms.eigensolvers import NumPyEigensolver

def create_problem(problemtype, problemsize, kwargs):
    if problemtype=="maxcut" or problemtype=="vertcover" or problemtype=="MIS":
        graphpath=create_nx_Rgraph(problemsize, kwargs["degree"], kwargs["seed"])
        return graphpath
    elif problemtype=="max3sat":
        rng = np.random.default_rng(seed=kwargs["seed"]) #if None, then fresh, unpredictable entropy will be pulled from the OS, else gets set to specified seed
        qubopath=create_3sat_formula(rng,problemsize,kwargs["alpha"],3,None)
        return qubopath
    else:
        print(f"sorry, the problem type {problemtype} has not been implemented yet")
        return None

def load_problem(problemtype, problemsize, problempath):
    print(problempath)
    if problemtype=="maxcut" or problemtype=="vertcover" or problemtype=="MIS":
        G, op, off = load_graph_as_ising(problempath,problemsize,problemtype)
        return op, off
    elif problemtype=="max3sat":
        op, off = load_qubo_as_ising(problempath)
        return op, off
    else:
        print(f"sorry, the problem type {problemtype} has not been implemented yet")
        return None

def create_nx_Rgraph(order, degree, seed):
    '''
    seed: numpy rng seed for creating nx graph
    degree: degree of graph
    order: order of graph/number of qubits
    '''

    G = nx.random_regular_graph(degree, order, seed)

    w = np.zeros([order, order])
    for i in range(order):
        for j in range(order):
            temp = G.get_edge_data(i, j, default=0)
            if temp != 0:
                w[i, j] = 1

    graph_outpath = f"./graphs/u{degree}R_{order}n_{seed}.gpickle"
    nx.write_gpickle(G,graph_outpath)
    return graph_outpath

def create_3sat_formula(rng,problemsize, alpha, vars_per_clause=3, custom_sat_formula=None):
    #problemsize==number of qubits(m+n) m=number of clauses, n=number of valiables
    if custom_sat_formula!=None:
        sat_formula =custom_sat_formula
        num_clauses = len(sat_formula)
        num_vars = problemsize-num_clauses
        possible_alpha = num_clauses/num_vars
    else:
        num_vars = max(3,int(problemsize/(alpha+1)))
        num_clauses = problemsize - num_vars
        possible_alpha = num_clauses/num_vars # alpha may be rounded to the closest possible alpha if assignment is impossible, resulting in new alpha
        if possible_alpha!=alpha:
            print(f"alpha {alpha} not possible, using alpha={possible_alpha} instead  with {num_vars} variables and {num_clauses} clauses")
        else:
            print(f"using alpha={alpha} with {num_vars} variables and {num_clauses} clauses")
        #create random sat formula
        sat_formula=[]
        for _ in range(num_clauses):
            literals = rng.choice(range(0,num_vars), size=vars_per_clause, replace=False)
            truth_val = rng.choice([True,False], size=vars_per_clause, replace=True)
            sat_formula.append(tuple(zip(literals,truth_val)))
    print(f"SAT Formula:\n{sat_formula}")

    # the following has been adapted from https://github.com/lfd/qsw_23_codesign_scalability/
    # several changes were made however to create a standard mapping
    QUBO=np.zeros((problemsize,problemsize))
    # goal: sum_i=1^m​((1+wi​)(yi1​+yi2​+yi3​)−yi1​yi2​−yi1​yi3​−yi2​yi3​−2wi)"
    # sum over m (that is, the number of clauses)
    for i, clause in enumerate(sat_formula):
            c_i = num_vars + i
            QUBO[c_i][c_i] += 2

            #part one: (1+wi​)(yi1​+yi2​+yi3​)
            for j in clause:
                l_ij =j[0]
                val_l_ij = j[1]
                if val_l_ij: # if var j in clause i is True
                    QUBO[l_ij][l_ij] -= 1.
                    QUBO[l_ij][c_i] -= .5
                    QUBO[c_i][l_ij] -= .5
                if not val_l_ij: # if var j in clause i is False
                    QUBO[l_ij][l_ij] += 1.
                    QUBO[c_i][c_i] -= 1.
                    QUBO[l_ij][c_i] += .5
                    QUBO[c_i][l_ij] += .5
            
            #part two: −yi1​yi2​−yi1​yi3​−yi2​yi3​−2wi
            for (item1, item2) in itertools.combinations(clause, 2): #iterate though variables pairwise (i,j)
                idx1 = item1[0]
                idx2 = item2[0]
                val1 = item1[1]
                val2 = item2[1]

                if val1 and val2: # if i=j, i,j>0
                    QUBO[idx1][idx2] += .5
                    QUBO[idx2][idx1] += .5
                if not val1 and val2: # if i!=j, i<0<j
                    QUBO[idx2][idx2] += 1.
                    QUBO[idx1][idx2] -= .5
                    QUBO[idx2][idx1] -= .5
                if val1 and not val2: # if i!=j, i>0>j
                    QUBO[idx1][idx1] += 1.
                    QUBO[idx1][idx2] -= .5
                    QUBO[idx2][idx1] -= .5
                if not val1 and not val2: # if i=j, i,j<0
                    QUBO[idx1][idx2] += 1.
                    QUBO[idx2][idx1] += 1.
                    QUBO[idx1][idx1] -= 1.
                    QUBO[idx2][idx2] -= 1.

    print(f"SAT formula converted to QUBO: \n{QUBO}")
    alphastr = str(possible_alpha).replace(".","")
    outpath = f"./qubos/{vars_per_clause}sat_{problemsize}qubits_n{num_vars}_m{num_clauses}_alpha{alphastr}"
    np.save(outpath,QUBO)
    return outpath

def load_qubo_as_ising(path):
    qubo = np.load(f"{path}.npy")
    Q = qubo / np.max(np.abs(qubo)) #normalize qubo
    n = Q.shape[0]
    qp = QuadraticProgram()
    linear_terms=[]
    for i in range(n):
        qp.binary_var(name=f"x_{i}")
        linear_terms.append(0.5 * Q[i, i] + 0.25 * np.sum([(Q[i, j] + Q[j, i]) for j in range(n) if j != i]))
    qp.minimize(linear=linear_terms,quadratic=0.25*Q)
    ising = qp.to_ising()
    pauli_op, offset = ising
    return pauli_op, offset

def load_graph_as_ising(graph_path,n,problem_type="maxcut"):
    # load graph
    G = nx.read_gpickle(graph_path)

    # quadratic problem from weight matrix
    if problem_type == "maxcut":
        #get weight matrix
        w = np.zeros([n, n])
        for i in range(n):
            for j in range(n):
                temp = G.get_edge_data(i, j, default=0)
                if temp != 0:
                    w[i, j] = 1

        max_cut = Maxcut(w)
        qp = max_cut.to_quadratic_program()
        ising = qp.to_ising()
    elif problem_type == "vertcover":
        vert_cover = VertexCover(G)
        qp = vert_cover.to_quadratic_program()
        qp2qubo = QuadraticProgramToQubo()
        qubo = qp2qubo.convert(qp)
        ising = qubo.to_ising()
    elif problem_type == "MIS":
        mis = StableSet(G)
        qp = mis.to_quadratic_program()
        qp2qubo = QuadraticProgramToQubo()
        qubo = qp2qubo.convert(qp)
        ising = qubo.to_ising()
    pauli_op, offset = ising
    return G, pauli_op, offset

def get_energy_eigenvalues(n,op,outpath):
    eig_solver = NumPyEigensolver(k=2**n)
    print(f"computing eigenvalues for operator {op}")
    eig_values = eig_solver.compute_eigenvalues(op).eigenvalues
    print(f"Eigenvalues:{eig_values}; Min:{float(eig_values.min())}; Max:{float(eig_values.max())}")
    print(f"saving results as json in {outpath}...")
    data = {"energy_eigenvalues":[float(val) for val in eig_values], "min_energy":float(eig_values.min()), "max_energy":float(eig_values.max())}
    with open(f"{outpath}/eigenvalues.json", 'w') as f:
       json.dump(data, f)

def get_min_energy_eigenvalue(op,outpath,fname="min_eigenvalue"):
    eig_solver = NumPyEigensolver(k=1)
    print(f"computing minimum eigenvalue for {op}")
    res = eig_solver.compute_eigenvalues(op)
    eig_val = res.eigenvalues
    eig_state = res.eigenstates
    print(f"Eigenvalue:{eig_val}; Eigenstate: {eig_state}")
    print(f"saving results as json in {outpath}...")
    data = {"min_energy":float(eig_val.real),"eigenstate":str(eig_state)}
    os.makedirs(outpath, exist_ok=True)
    with open(f"{outpath}/{fname}.json", 'w') as f:
       json.dump(data, f)

def save_only_if_higher(save_val,problem,method,qubits,seed):
   save_path = f"{problem}/max_energies"
   if os.path.isdir(save_path)==False:
       os.makedirs(save_path, exist_ok=True)
   fname=f"n{qubits}_seed{seed}.json"
   path=f"{save_path}/{fname}"
   if os.path.isfile(path)==True:
      print("exists")
      with open(path) as f:
          data = json.load(f)
          prev_max = data["max_energy"]
      if save_val>prev_max:
         save(path,save_val,method)
   else:
      save(path,save_val,method)
   
def save(path,val,method):
   save_data = {"max_energy":val,"found_by":method}
   with open(path, 'w') as f:
      json.dump(save_data, f)

def read_and_update_max_vals(dir_name,out_name,meth,qubits,seed,resolution):
    json_res = []

    for file_name in [file for file in os.listdir(dir_name) if file.endswith('.json')]:
        with open(dir_name + "/" + file_name) as json_file:
            if file_name!="offset.json" and file_name!="eigenvalues.json" and file_name!="opt_data.json":
                data = json.load(json_file)
                json_res.append(data["expecation_values"])

    mat = np.ndarray((len(json_res),resolution,resolution))

    for i in range(len(json_res)):
        expectation_mat = np.asmatrix(json_res[i])
        mat[i] = expectation_mat
    max = np.amax(mat)

    save_only_if_higher(max,out_name,meth,qubits,seed)

def check_isomorphism_of_seeds(degree,order,seeds):
    for i, seed in enumerate(seeds):
        for j, other_seed in enumerate(seeds):
            G1 = nx.random_regular_graph(degree, order, seed)
            G2 = nx.random_regular_graph(degree, order, other_seed)
            if i!=j and nx.is_isomorphic(G1,G2):
                print(f"graphs with seeds {seeds[i]} and {seeds[j]} are isomorphic")
                return False
    print(f"graphs with seeds {seeds} are not isomorphic to oneanother")
    return True

def makeLRschedule(p,delta_beta,delta_gamma,beta_neg=False,gamma_neg=False):
    gamma = []
    beta = []
    for i in range(p):
        beta_i = (1-i/p)*delta_beta
        gamma_i = ((i+1)/p)*delta_gamma
        beta.append(beta_i)
        gamma.append(gamma_i)
    neg_ramp = ""
    if beta_neg:
        beta = [-1*b for b in beta]
        neg_ramp = neg_ramp + "_neg_beta"
    if gamma_neg:
        gamma = [-1*g for g in gamma]
        neg_ramp = neg_ramp + "_neg_gamma"
    params = beta+gamma
    data = {"opt_params":params}

    rampname = f"linear_ramp_p{p}{neg_ramp}.json"
    with open(f"./fixed_params/{rampname}", 'w') as f:
       json.dump(data, f)

def load_offset(problemtype,problempath,problemsize,outpaths):
    pauli_op, offset = load_problem(problemtype,problemsize,problempath)
    offsetdata={"offset":offset}
    pauli_list = pauli_op.to_list(array=True)
    print(f"saving offsets and operator for problem corresponding to {outpaths}...")
    for path in outpaths:
        with open(f"{path}/offset.json", 'w') as f:
            json.dump(offsetdata, f)
        np.save(f"{path}/pauli_op", pauli_list)

    








