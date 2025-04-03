from problem_util import create_problem, load_problem
from perform_scan import run_scan_with_iterative_parameter_fixing

from plot_scan import meshplot_series, voxelplot, fourier_landscape_series, lineplot, parse_opt_expectation_values
from optimize_from_scan import optimize_from_inits, read_best_existing_params
from problem_util import get_energy_eigenvalues, get_min_energy_eigenvalue, read_and_update_max_vals

import os
import argparse
import numpy as np

parser = argparse.ArgumentParser()
# general parameters
parser.add_argument("-pr","--problemtype",type=str,default="maxcut",help="problem type; currently supported: maxcut, max3sat, vertcover, MIS")
parser.add_argument("-n","--problemsize",type=int,help="[required] graph order / problem size")
parser.add_argument("-meth","--methodname",type=str,default="any",help="method name; results conducted with the same method are stored in a subdirectory with this name. Aside from this, this argument has no effect")
# scan parameters
parser.add_argument("-d","--degree",type=int,default=3,help="problem graph degree. not necessary when using -plt / -pltavg")
parser.add_argument("-a","--alpha",type=float,default=None,help="sat alpha (clauses/valiables). not necessary when using -plt / -pltavg")
parser.add_argument("-p","--target_p",type=int, default=1,help="qaoa depth up to which a layer-by-layer based scan is performed. not necessary when using -plt / -pltavg")
parser.add_argument("-pstep","--stepsize_p",type=int, default=1,help="amount by which to increase the circuit depth for each scan iteration. defaults to 1. not necessary when using -plt / -pltavg, but may lead to wrong plotlabvelling when ommited")
parser.add_argument("-pstart","--starting_p",type=int, default=1,help="qaoa depth from which to start. defaults to 1. not necessary when using -plt / -pltavg, but may lead to wrong plotlabvelling when ommited")
parser.add_argument("-resume","--resume_from",type=str, default=None,help="path to results from which to continue sequential scan, must be used with --starting_p")
parser.add_argument("-resume_stp","--resume_startparams",action=argparse.BooleanOptionalAction,help="to resume a sequential or optimised scan, add this arg and pass the respective params to the -stp arg, in conjunction with the -resume arg and the -startp arg")
parser.add_argument("-r","--resolution",type=int, default=32,help="grid resolution (results in grid with resolution * resolution points and step size len(bounds)/resolution). not necessary when using -plt / -pltavg")
parser.add_argument("-sy","--symmetry",type=str,default="no_symmetry",help="string, either 'no_symmetry' (beta range [-pi,pi], gamma range [-pi,pi]), 'maxcut_symmetry' (beta range [-pi/4,pi/4], gamma range [-pi/2,pi/2]), 'general_symmetry' (beta range [-pi/2,pi/2], gamma range [-pi,pi]) 'both_maxcut' (no_symmetry and maxcut_symmetry), 'both_general' (no_symmetry and maxcut_symmetry). not necessary when using -plt / -pltavg")
parser.add_argument("-sh","--shots",type=int,default=None,help="[optional] number of shots for the estimator (none = calculate exact expectation value)")
parser.add_argument("-prs","--problemseed",type=int,default=None,help="[optional] seed for generating the graph/problem")
parser.add_argument("-prp","--problempath",default=None,help="[optional] path to graph/problem for which to scan/optimize qaoa landscape; alternative to creating a graph with -n -d and -gs flags")
parser.add_argument("-sis","--simseed",type=int,default=None,help="[optional] seed for statevector simulator")
# optimize from scan parameters
parser.add_argument("-opt","--optimizer",default=None, help="scypi optimizer name as string (e.g 'COBYLA', 'BFGS', etc; see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html). If not None optimizes from the best parameters found in scan (or existing scan data if -stp is not None) and then scans and plots the landscape around the optimized params")
parser.add_argument("-stp","--startparams",type=str,default=None,help="[optional] path to json results containing fixed parameters of the form {'opt_params':[params]} for which to scan the landscape, e.g. params found by a optimizer or best params found by scan. If you choose this option instead of performing a scan, make sure to specify the graph using -grp as well.")
# eigensolver
parser.add_argument("-eig","--eigensolver",action=argparse.BooleanOptionalAction,help="exactly solve the problem hamiltonian and save in a json file")
parser.add_argument("-upmax","--update_max",action=argparse.BooleanOptionalAction,help="uüpdate max found energies and save in a json file")
parser.add_argument("-eigall","--eigensolverall",type=str,help="Warning: computationally expensive, dont execute unless you are prepared for that. compute every eigenstate and eigenvalue of the problem hamiltonian and save in a json file; string to folder where results should be saved")
# plot parameters
parser.add_argument("-mskg","--maskgeq",type=float,default=None,help="[optional] mask expectation values greater or equal to this in second voxelplot")
parser.add_argument("-mskl","--maskleq",type=float,default=None,help="[optional] mask expectation values less or equal to this in second voxelplot")
parser.add_argument("-azh","--azimuth",type=float,default=-120,help="[optional] azimuthal angle of voxelplot")
parser.add_argument("-elv","--elevation",type=float,default=20,help="[optional] elevation / polar angle of voxelplot")
parser.add_argument("-plt","--plotexistingresults",type=str,nargs="*",default=None,help="[optional] multiple strings from which to plot results. If None: performs scans; if not None: will skip the scans")
parser.add_argument("-pltdir","--plotexistingdirectory",type=str,default=None,help="string which designates a directory from which to plot all results contained in this dir. When using this argument, -n has no effect, however -sy must be specified and the same for all results in the directory. If None: performs scans; if not None: will skip the scans")
parser.add_argument("-pltst","--plotmultipleresultsstats",type=str,default=None,help="[optional] path to dir containg multiple results. will calculate average of expecation/max_expecation of all result directories contained in this dir. can only be used with already existing results")
parser.add_argument("-pltsti","--plotmultipleresultsstatsinfo",type=str,default="unspecified",help="[optional] information about problem sizes to display in average plot title and filename when plotting multiple results")
parser.add_argument("-fft","--fourierplot", action=argparse.BooleanOptionalAction, help="[optional] creates fourier plot from scanned or loaded params")
parser.add_argument("-pltpd","--plotparamdata", action=argparse.BooleanOptionalAction, help="[optional] creates lineplots with arrangement and quality of the parameters given by the method")
parser.add_argument("-fftu","--fouriermaxfrequ",type=int, default=6,help="[optional] maximum u (transformed betas) frequency components to show in plot")
parser.add_argument("-fftv","--fouriermaxfreqv",type=int, default=12,help="[optional] maximum v (transformed gammas) frequency components to show in plot")
parser.add_argument("-skp","--skipplot", action=argparse.BooleanOptionalAction, help="[optional] skips plotting, unless -plt / -pltavg is not None")
parser.add_argument("-skpv","--skipvoxel", action=argparse.BooleanOptionalAction, help="[optional] skips voxelplot")
parser.add_argument("-png","--savepng", action=argparse.BooleanOptionalAction, help="[optional] saves plot in a png file in addition to the normal pdf")
parser.add_argument("-svg","--savesvg", action=argparse.BooleanOptionalAction, help="[optional] saves plot in a svg file in addition to the normal pdf")
parser.add_argument("-pgf","--savepgf", action=argparse.BooleanOptionalAction, help="[optional] saves plot in a pgf file in addition to the normal pdf")
parser.add_argument("-tikz","--savetikz", action=argparse.BooleanOptionalAction, help="[optional] saves plot in a tikz tex file in addition to the normal pdf")
parser.add_argument("-ps","--saveps", action=argparse.BooleanOptionalAction, help="[optional] saves plot in a postscript file in addition to the normal pdf")

if __name__ == '__main__':
    args = parser.parse_args()
    # set up bounds
    if args.symmetry=="no_symmetry":
        bounds_dict={"no_symmetry":[[-np.pi,np.pi],[-np.pi,np.pi]]}
    elif args.symmetry=="maxcut_symmetry":
        bounds_dict={"maxcut_symmetry":[[-np.pi/4,np.pi/4],[-np.pi/2,np.pi/2]]}
    elif args.symmetry=="general_symmetry":
        bounds_dict={"general_symmetry":[[-np.pi/2,np.pi/2],[-np.pi,np.pi]]}
    elif args.symmetry=="both_maxcut":
        bounds_dict={"no_symmetry":[[-np.pi,np.pi],[-np.pi,np.pi]],"maxcut_symmetry":[[-np.pi/4,np.pi/4],[-np.pi/2,np.pi/2]]}
    elif args.symmetry=="both_general":
        bounds_dict={"no_symmetry":[[-np.pi,np.pi],[-np.pi,np.pi]],"general_symmetry":[[-np.pi/2,np.pi/2],[-np.pi,np.pi]]}
    elif args.plotexistingresults or args.plotmultipleresultsstats:
        bounds_dict={"no_symmetry":[[-np.pi,np.pi],[-np.pi,np.pi]],"maxcut_symmetry":[[-np.pi/4,np.pi/4],[-np.pi/2,np.pi/2]],"general_symmetry":[[-np.pi/2,np.pi/2],[-np.pi,np.pi]]}


    scan_paths = None
    if args.problempath:
        problem_path = args.problempath
    elif (not (args.plotmultipleresultsstats or args.plotexistingresults)) or (args.plotexistingresults and args.optimizer):
        if not args.problemsize:
            print(f"please provide the path (-prp) to the problem problem arguments (-n or other problem specific parameters) for the scan in {args.plotexistingresults}")
        prob_kwargs = {"degree":args.degree,"alpha":args.alpha,"seed":args.problemseed}
        problem_path = create_problem(args.problemtype,args.problemsize,kwargs=prob_kwargs)


    if args.plotexistingresults:
        print("Skipping scan and loading existing results...")
        scan_paths = args.plotexistingresults
    elif args.plotexistingdirectory or args.startparams:
            scan_paths = []
            if args.plotexistingdirectory:
                paths = os.listdir(args.plotexistingdirectory)
                for path in paths:
                    scan_paths.append(args.plotexistingdirectory + "/" + path)
    elif args.startparams:
        scan_paths = []
    elif not (args.plotmultipleresultsstats or args.eigensolver or args.eigensolverall):
        print("scanning")
        scan_paths = []
        for key, values in bounds_dict.items():
            beta_bounds = values[0]
            gamma_bounds= values[1]
            scan_result_path = run_scan_with_iterative_parameter_fixing(args.problemtype,prob_kwargs,args.methodname,problem_path,args.problemsize,args.target_p,args.resolution,beta_bounds,gamma_bounds,None,args.shots,args.simseed,key,args.stepsize_p,args.starting_p,args.resume_from,args.resume_startparams)
            scan_paths.append(scan_result_path) 

    if (args.optimizer or args.startparams) and not args.plotmultipleresultsstats:
        initial_params = []
        if args.startparams:
            if args.optimizer:
                results = optimize_from_inits(args.problemtype, args.problemsize,args.optimizer,problem_path,args.startparams,args.simseed)
                initial_params.append(results["opt_params"])
                
            else:
                params = read_best_existing_params(args.startparams, single_file=True)
                initial_params.append(params)
        else:
            for path in scan_paths:
                results = optimize_from_inits(args.problemsize,args.optimizer,problem_path,path,args.simseed)
                initial_params.append(results["opt_params"])
        for key, values in bounds_dict.items():
                beta_bounds = values[0]
                gamma_bounds= values[1]
                for point in initial_params:
                    scan_path = run_scan_with_iterative_parameter_fixing(args.problemtype,prob_kwargs,args.methodname,problem_path,args.problemsize,args.target_p,args.resolution,beta_bounds,gamma_bounds,point,args.shots,args.simseed,f"{args.optimizer}_optimized_{key}",args.stepsize_p,args.starting_p,args.resume_from,args.resume_startparams)
                    scan_paths.append(scan_path)
                    if args.optimizer:
                        parse_opt_expectation_values(scan_path,results["opt_params"].tolist(),results["opt_value"].tolist(),args.target_p,args.stepsize_p)
    
    if args.eigensolver or args.update_max:
        if args.problemtype=="max3sat":
            if (args.alpha>4.9 or args.alpha<=3.5):
                subdir="max3sat/easy"
            else:
                subdir="max3sat/hard"
        else:
            subdir=args.problemtype
        if args.eigensolver:
            pauli_op, offset = load_problem(args.problemtype,args.problemsize,problem_path)
            get_min_energy_eigenvalue(pauli_op,f"./results/reproduction/{subdir}/eigensolver",fname=f"min_eigenvalue_{args.problemseed}_{args.problemsize}qubits")
        if args.update_max:
            for path in scan_paths:
                read_and_update_max_vals(path,f"./results/reproduction/{subdir}",args.methodname,args.problemsize,args.problemseed,args.resolution)

        
    elif args.eigensolverall:
        pauli_op, offset = load_problem(args.problemtype,args.problemsize,problem_path)
        get_energy_eigenvalues(args.problemsize,pauli_op,args.eigensolver)

    if not args.skipplot:
        print("-----------------------------------------------------------------------------------------------------------------------------------------------")
        print(f"Plotting results for scans contained in {scan_paths if scan_paths else args.plotmultipleresultsstats}...")
        print("-----------------------------------------------------------------------------------------------------------------------------------------------")
        
        saveformats = ["pdf"]
        otherformats = ["png","svg","pgf","tex","ps"]
        formatargs = [args.savepng,args.savesvg,args.savepgf,args.savetikz,args.saveps]
        for i in range(5):
            if formatargs[i]:
                saveformats.append(otherformats[i])

        if scan_paths:
            print(f"creating single instance plots...")
            plot_best = False
            for i, path in enumerate(scan_paths):
                if path.endswith("no_symmetry"):
                    symmetry = "no_symmetry"
                elif path.endswith("maxcut_symmetry"):
                    symmetry = "maxcut_symmetry"
                elif  path.endswith("general_symmetry"):
                    symmetry = "general_symmetry"
                if path.endswith(f"optimized_{symmetry}"):
                    title_base = f" for QAOA with increasing p for a {str(args.problemsize)} qubit problem instance, parameters optimized using the {args.optimizer} method"
                    title = "Expectation Values" + title_base
                    title_fourier = "Expectation Values" + title_base
                    plot_best = True
                else:
                    title = None
                    title_fourier = None
                beta_bounds = bounds_dict[symmetry][0]
                gamma_bounds= bounds_dict[symmetry][1]
                print(f"---Meshplot ({i+1}/{len(scan_paths)})---")
                meshplot_series(path,beta_bounds,gamma_bounds,args.problemsize,"viridis",title,format=saveformats,pstep=args.stepsize_p,pstart=args.starting_p)
                if not args.maskgeq:
                    args.maskgeq = -((args.problemsize/4)**2)*0.35
                if not args.skipvoxel:
                    print(f"---Voxelplot ({i+1}/{len(scan_paths)})---")
                    voxelplot(path,beta_bounds,gamma_bounds,args.problemsize,args.maskleq,args.maskgeq,args.elevation,args.azimuth,"viridis",title,format=saveformats,plot_best=plot_best)
                if args.fourierplot:
                    print(f"---Fourier Landscape plot ({i+1}/{len(scan_paths)})---")
                    fourier_landscape_series(path,args.problemsize,args.fouriermaxfrequ,args.fouriermaxfreqv,"viridis",title_fourier,format=saveformats)
        
        if args.plotmultipleresultsstats:
            paths = os.listdir(args.plotmultipleresultsstats)
            if paths[0].endswith("no_symmetry"):
                sym = "no_symmetry"
            elif paths[0].endswith("maxcut_symmetry"):
                sym = "maxcut_symmetry"
            elif paths[0].endswith("general_symmetry"):
                sym = "general_symmetry"
            else:
                print(f"Error: cannot average the results due to different symmetries; please make sure all results contained in {args.plotmultipleresultsstats} were scanned using the same symmetry range")
            beta_bounds = bounds_dict[sym][0]
            gamma_bounds= bounds_dict[sym][1]

            print("aggregating results and creating average/stdev plots from resulting data")
            title = "Average of normalized Expectation values for QAOA with increasing p over " + str(len(paths)) + " problem instances of size(s) " + args.plotmultipleresultsstatsinfo
            print("---Meshplot (1/2)---")
            meshplot_series(args.plotmultipleresultsstats,beta_bounds,gamma_bounds,args.plotmultipleresultsstatsinfo,"plasma",title,True,format=saveformats,pstep=args.stepsize_p,pstart=args.starting_p)
            if not args.maskleq:
                    args.maskleq = 0.8
            if not args.skipvoxel:
                print("---Voxelplot (1/2)---")
                voxelplot(args.plotmultipleresultsstats,beta_bounds,gamma_bounds,args.plotmultipleresultsstatsinfo,args.maskleq,None,args.elevation,args.azimuth,"plasma",title,True,format=saveformats)
            print("creating 2 plots from standard-deviation data...")
            title = "Standard deviation of normalized Expectation values for QAOA with increasing p over " + str(len(paths)) + " problem instances of size(s) " + args.plotmultipleresultsstatsinfo
            print("---Meshplot (2/2)---")
            meshplot_series(args.plotmultipleresultsstats,beta_bounds,gamma_bounds,args.plotmultipleresultsstatsinfo,"cividis",title,False,True,format=saveformats,pstep=args.stepsize_p,pstart=args.starting_p)
            if not args.maskgeq:
                    args.maskgeq = 0.05
            if not args.skipvoxel:
                print("---Voxelplot (2/2)---")
                voxelplot(args.plotmultipleresultsstats,beta_bounds,gamma_bounds,args.plotmultipleresultsstatsinfo,None,args.maskgeq,args.elevation,args.azimuth,"cividis",title,False,True,format=saveformats)

            if args.plotparamdata:
                print("aggregating results and calculating average/stdev of only the parameters given by the method")
                lineplot(args.plotmultipleresultsstats,args.problemtype,args.methodname,args.target_p,args.starting_p,args.stepsize_p,add_offset=True,format=saveformats)




