# QAOAscan Framework: Advanced usage - passing arguments to main.py directly

## Run a simple scan

To run a scan, run main.py with the desired arguments; the only strictly required argument is -n (the number of nodes/qubits), but you will most likely want to specify some further parameters, as otherwise defaults, which may not align with your requirements will be used.
every scan requires the arguments:
 - -n, --order (number of nodes/qubits)
 - -d, --degree (degree of the problem graph, default 3)
 - -p. --target_p (target QAOA depth to perform scan up to, default 1)
 - -r, --resolution (resolution of the scan; creates a grid of size resolution x resolution, default 32)
 - -sy, --symmetry (type of symmetry, this determines the parameter bounds. must be either "maxcut_symmetry,  "no_symmetry" or "both,  default "no_symmetry")

For example, to run a 32x32 pixel scan of QAOA of depth 5 on a 4-node 3-regular graph without accounting for symmetries (that is, gamma and beta are between -pi and pi), run the following command:
```
python src/main.py -n 4 -d 3 -p 5 -r 32 -sy "no_symmetry"
```
This produces scan data which results in the following plots:
![example_simple_scan_no_symmetry_meshplot](results/examples/example1.png)
![example_simple_scan_no_symmetry_voxelplot](results/examples/example2.png)

running the same scan with "maxcut_symmetry" as -sy argument, accounts for symmetries for maxcut on regular graphs (that is, gamma is between -pi/2 and pi/2 and beta is between -pi/4 and pi/4)
```
python src/main.py -n 4 -d 3 -p 5 -r 32 -sy "maxcut_symmetry"
```
yields the following results:
![example_simple_scan_maxcut_symmetry_meshplot](results/examples/example3.png)
![example_simple_scan_maxcut_symmetry_voxelplot](results/examples/example4.png)

### Run a scan from an exising .gpickle nx graph
If you want to specify a certain graph instead of creating a new one, you can do this using the -prp (--problempath) argument
```
python src/main.py -n 4 -prp "graphs/my_graph.gpickle" -p 5 -r 32 -sy "maxcut_symmetry"
```

### Optional scan arguments 
You can specify some further parameters for the scan, though generally this is not needed:
 - -sh, --shots (number of shots for the estimator; if None calculate exact expectation value. default None)
 - -ps, --problemseed (rng seed for generating the graph. default 123456)
 - -sim, --simseed (seed for qiskit.algorithm_globals. default None)

## Plot existing scans
If you already have results from a previous scan, you can (re)plot these results by passing the name of the directory, which containes these results, to the -plt (--plotexistingresults) arguments.
if this argument is not None, any scans will be skipped completely
Example:
```
python src/main.py -n 4 -plt "results/my_results"
```
you can also (re)plot all results in a directory with the -pltdir (--plotdirectory) argument. When doing this, make sure that all scans in the directory have the same symmetry and specify this using the -sy argument.
Example:
```
python src/main.py -n 4 -sy "maxcut_symmetry" -pltdir "results"
```

### Optional plot arguments
You can specify some further parameters for the plots, though generally this is not needed:
 - -skp, --skipplot (boolean argument: if True, skips plotting, unless -plt / -pltavg is not None. default False)
 - -mskg, --maskgeq (mask expectation values greater or equal to this in second voxelplot. default None)
 - -mskl, --maskleq (mask expectation values less or equal to this in second voxelplot. default None)
 - -azh, --azimuth (azimuthal angle of voxelplot. default -120)
 - -elv, --elevation (elevation / polar angle of voxelplot. default 20)
 - -png, --savepng (boolean argument: if True, saves plot in a png file in addition to the normal pdf. default False)
 - -svg, --savesvg (boolean argument: if True, saves plot in a svg file in addition to the normal pdf. default False)
 - -ps, --saveps (boolean argument: if True, saves plot in a postscript file in addition to the normal pdf. default False)

## Optimize Scan parameters and scan QAOA for the resulting optimized parameters
In order to optimize starting from the best parameters found in a scan and then scan QAOA again for the resulting optimized parameters, specify an optimizer to use with the -opt (--optimizer) argument.
this can be any string accepted by the "method" parameter of the [scipy minimize()](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html) function (e.g "COBYLA", "BFGS")
for example:
```
python src/main.py -n 4 -d 3 -p 5 -r 32 -sy "maxcut_symmetry" -opt "COBYLA"
```
produces the results
![example_optimizer_scan_meshplot](results/examples/example5.png)
![example_optimizer_scan_voxelplot](results/examples/example6.png)

### Optimize from existing best params
Instead of running a scan from which to optimize, you can also specify the starting parameters directly using the -stp (--startparams) arguments:
this argument msut be a path to json results containing fixed parameters of the form {'opt_params':[params]} for which to scan the landscape, e.g. params found by a optimizer or best params found by scan. If you choose this option instead of performing a scan, make sure to specify the graph using -prp as well.
```
python src/main.py -n 4 -p 5 -prp "graphs/u3R_4n_123456.gpickle" -stp "fixed_params/fixed_parameters_p5.json" -opt "COBYLA"
```

## Plot Averages and Standard deviations of multiple existing scans
In order to plot the results of mutiple scans, you can pass the directory name (where these scans are contained) to the -pltst (--plotmultipleresultsstats) argument, which will will calculate average and standard deviation of expectation/max_expectation for the scans in this directory.
Important: Make sure the scans have the same symmetries and resolutions, otherwise the results may be plotted erroneously.
Using the -pltsti (--plotmultipleresultsstatsinfo) argument, you can manually specify the range of problem sizes, which will be included in the plot title.
```
python src/main.py -n -pltst "results" -pltsti "4to16"
```
produces plots for the average
![example_avg_scan_meshplot](results/examples/example7.png)
![example_avg_scan_voxelplot](results/examples/example9.png)
and the standard deviation
![example_std_scan_meshplot](results/examples/example8.png)
![example_std__scan_voxelplot](results/examples/example10.png)

## Plot Fourier-Landcape scan of scan results
In order to plot a Fourier-Landscape Plot (see [Stęchły et al.](https://arxiv.org/abs/2305.13594)), you can add the -fft (--fourierplot) argument, which creates fourier plot from scanned or loaded params, by applying a fourier transformation to the parameter grid and plotting the expectation values of the reulting frequency components; beta parameters get transoformed to u and gamma parameters get tranformed to v in the frequency domain representation.

```
python src/main.py -n 4 -d 3 -p 5 -r 32 -sy "maxcut_symmetry" -fft
```
or (when plotting from existing scan results)
```
python src/main.py -n 4 -plt "results/my_results" -fft
```
results in
![example_fourier__scan_meshplot](results/examples/example11.png)

Other arguments (use as needed):
 - -fftu, --fouriermaxfrequ (maximum u (transformed betas) frequency components to show in plot)
 - -fftv, --fouriermaxfreqv (maximum v (transformed gammas) frequency components to show in plot)

 ## Check isomorphism of a set of graphs

 In order to ensure that the graphs you are running scans for are not isomorphic to oneanother, you can check for isomorphisms by passing the graph degree, order and a list graph-seeds you want to check to the check_isomorphism_of_seeds function of the problem_util module.

 For example, to check if the 10-vertex 3-regular graph generated by the seeds 123456 and 123457 are isomorphic to oneanother, run the command
 ```
 python -c "from src.problem_util import check_isomorphism_of_seeds; check_isomorphism_of_seeds(3,10,[123456,123457])"
 ```
 this results in the output "graphs with seeds [123456, 123457] are not isomorphic to oneanother"