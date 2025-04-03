# QAOAscan Framework

A Framework for visualizing the QAOA cost landscape, which performs scans by creating a grid over values of beta and gamma, and, for p>1, fixes the parameters of previous layers (to either the best scanned parameters, the best parameters found by an optimizer, or explicitly specified parameters.) In this way we can get a glimpse/cross-section of the landscape near these parameters, in order to draw conclusions about the characteristics of QAOA parameters, such as parameter concentrations, patterns and symmetries.

## Installation

##### 1. Clone repo

via SSH
```
git clone ssh://git@gitlab.oth-regensburg.de:IM/lfd/dissemination/ma/2024_eichenseher_vincent/qaoascan_framework.git
```
via HTTPS
```
git https://gitlab.oth-regensburg.de:IM/lfd/dissemination/ma/2024_eichenseher_vincent/qaoascan_framework.git
```

##### 2. (optional) Create virtual environment
```
python -m venv .venv
source .venv/bin/activate
```
In our experiments we used python 3.9.5

##### 3. Install requirements from requirements.txt
```
pip install -r requirements.txt
```

The general usage of this reproduction package is as follows:

## Run Experiments
To run experiments, call the bash script `run_from_config.sh` and pass the path(s) of one or more config file(s) as arguments
```
bash run_from_config.sh conf.config
```
```
bash run_from_config.sh conf1.config conf2.config conf3.config    
```

## Plot results
The results for each individual instance are plotted automatically after each experiment is run; To plot the average results over multiple instances call the bash script `plot_from_config.sh` and pass the path(s) of the same config files that were previously given as arguments to run_from_config.sh
```
bash plot_from_config.sh conf.config
```
```
bash plot_from_config.sh conf1.config conf2.config conf3.config    
```

## Executing both plot and run
when executing both the `plot_from_config.sh` plot `run_from_config.sh` scripts on one line, the operators ‘;’, ‘&&’, or ‘||’ may be used, but executing the commands asynchronously with '&' will cause the plot script to run before the results are present.
One possible way of running both commands is:
```
bash run_from_config.sh conf.config ; bash plot_from_config.sh "$_"
```
or
```
bash run_from_config.sh conf1.config conf2.config conf3.config ; bash plot_from_config.sh "$_"
```

### Simple Example
To familiarise yourself with the framework and parameter selection methods, you can run simple example experiment sets with few small problem instances, which do not take long to execute
To run one such example for a specific problem, use one of these commands:
```
bash run_from_config.sh config/test_maxcut.config
bash run_from_config.sh config/test_vertcover.config
bash run_from_config.sh config/test_max3sat.config
bash run_from_config.sh config/test_maxcut_higher_p.config 
```
or to execute the examples for all problems in one command:
```
bash run_from_config.sh config/test_maxcut.config config/test_vertcover.config config/test_maxcut_higher_p.config config/test_max3sat.config
```

To plot the average results, run the same command with `plot_from_config.sh` instead of `run_from_config.sh`as explained [here](#plot-results):
```
bash plot_from_config.sh config/test_maxcut.config
bash plot_from_config.sh config/test_vertcover.config
bash plot_from_config.sh config/test_max3sat.config
bash plot_from_config.sh config/test_maxcut_higher_p.config 
```
or
```
bash plot_from_config.sh config/test_maxcut.config config/test_vertcover.config config/test_maxcut_higher_p.config config/test_max3sat.config
```

To both run and plot, as explained [here](#executing-both-plot-and-run):
```
bash run_from_config.sh config/test_maxcut.config ; bash plot_from_config.sh "$_" 
bash run_from_config.sh config/test_vertcover.config ; bash plot_from_config.sh "$_"
bash run_from_config.sh config/test_max3sat.config ; bash plot_from_config.sh "$_"
bash run_from_config.sh config/test_maxcut_higher_p.config ; bash plot_from_config.sh "$_"
```
or:
```
bash run_from_config.sh config/test_maxcut.config config/test_vertcover.config config/test_maxcut_higher_p.config config/test_max3sat.config ; bash plot_from_config.sh "$_"
```

to check if you were able to reproduce the examples, compare results/reproduction with results/simple_examples, which contains the results of these simple example experiments run on our device(s)

### Reproduce Paper Experiments
To reproduce all the results obtained in the paper experiments, run the following command:
```
bash run_from_config.sh config/paper_maxcut.config config/paper_vertcover.config config/paper_maxcut_higher_p.config config/paper_max3sat.config
```
and subsequently plot the average results with
```
bash plot_from_config.sh config/paper_maxcut.config config/paper_vertcover.config config/paper_maxcut_higher_p.config config/paper_max3sat.config 
```
or execute both these commands simultaneously:
```
bash run_from_config.sh config/paper_maxcut.config config/paper_vertcover.config config/paper_maxcut_higher_p.config config/paper_max3sat.config ; bash plot_from_config.sh "$_"
```
**Be aware that running all experiments takes quite some time (especially the max3sat instances)**, so if are only interested in some results it is recommended to create a modified version of the config data, which contains only the methods or instances of interest, and running that instead. Depending on the device you intend on running the experiments on, it may also be advisable to divide the configs into smaller chunks to run sequentially instead of running all experiments at once (the following section goes into how to create custom configuration files)

### creating custom Config-files to run your own Experiments
A config file is generally structured in the following way:
 - parameters pertaining to the combinatorial problem and the landscape scan are given in the form `{param:value}`
 - the methods to be used are specified with `(method#)` (#=0,1,2,... increasing with each additional method) and the method parameters are given in the form `param=value` on following lines and a newline following the parameter specification
 - the methods to be used are specified with `[instance#]` (#=0,1,2,... increasing with each additional method) and the instance parameters are given in the form `param=value` on following lines and a newline following the parameter specification 

below is an example using generalised values (for terms of the form value1|value2|value3 choose one of these options):
```
{problem:maxcut|vertexcover|max3sat|MIS}
{symmetry:no_symmetry|general_symmetry|maxcut_symmetry}
{depth:int}
{pstart:int}
{stepsize:int}
{resolution:int}

(method0)
dir_name=string
type=string|lr|opt

[instance0]
qubits=int
seed=int
degree=int
alpha=int

```
note that 
- when the parameter `type` of a method is set to any string that is not `lr` or `opt`,the sequential method will be run as default.
- when `{problem:maxcut}`,`{problem:vertcover}` or `{problem:MIS}` is specified, it will not matter what you specify for `alpha=int` in the individual instances.
- when `{problem:max3sat}` is specified, it will not matter what you specify for `degree=int` in the individual instances.

additially, when using the options `lr` or `opt` for the method parameter `type`, specify the following additional method parameters:
```
(method1)
dir_name=string
type=opt
optimiser=COBYLA|string
init=string
init_type=lr|string

(method2)
dir_name=string
type=lr
delta_beta=float
delta_gamma=float
negative_mixer=true|string

```
**note that for plotting to work all result data in a given directory has to have the same resolution, depth and symmetry, otherwise calculating averages and measures of spread will fail** for this reason, it is inadvisable to run custom configs which save results of different resolution, depth or symmetry to the same directories (`dir_name`). When running multiple configs, make sure this does not happen by specifying a different `dir_name` for experiments which differ in these parameters.
Furthermore, make sure that the parameters are set to sensible values, seeing as contradictory parameter specification, e.g when `pstart` is greater than `depth` may cause the script fail to skip those methods/instances. the parameters are explained in the section [Explanation of parameters](#explanation-of-parameters)

##### creating custom configs for plotting

Generally, you can just reuse the custom config file used to run the experiments as a config for the plotting script, however if for some reason you need a custom config (for example if you canged the basic structure or directory names of the results), you can specify a custom config used only for plotting. when doing this, you only need to specify `problem`,`symmetry`,`depth`,`pstart`,`stepsize` and the different methods to plot, each of which has to have at least the `dir_name` parameter specified. The general structure of such a config is:
```
{problem:maxcut|vertexcover|max3sat|MIS}
{symmetry:no_symmetry|general_symmetry|maxcut_symmetry}
{depth:int}
{pstart:int}
{stepsize:int}

(method0)
dir_name=string

(method1)
dir_name=string

```
the plot script will then plot the average and standard deviation plots for all instances contained in the directory "results/reproduction/`problem`/`dir_name`/" (in the case of `problem:max3sat` both "max3sat/easy/`dir_name`/" and "max3sat/hard/`dir_name`/" will be plotted separately)

##### Explanation of parameters

Problem and Scan parameters:
 - `problem`: the combinatorial optimisation problem that he instances in this config file belong to
 - `symmetry`: the type of symmetry we to restrict the parameter bounds to
 - `depth`: the target depth, up to which we scan the landscapes given by individual parameter-components 
 - `pstart`: the initial depth at which to start the scan and iteratively increase up to `depth`
 - `stepsize`: the amoiút by which we increase the QAOA circuit depth in each scan iteration
 - `resolution`: resolution of the grid over the parameters bounds, whcih we use to scan the landscapes given by individual parameter-components 

General Method Parameters
 - `dir_name`: the name of the subdirectory in results/reproduction where the scan results for his method witll be saved
 - `type`: method to use for fixing parameetrs of lower depth QAOA layers, see paper

 `type=opt` Method Parameters:
  - `optimiser`: optimiser to use; specified as scypi optimizer name (e.g 'COBYLA', 'BFGS', etc; see https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.minimize.html)
  - `init`: path to the directory from which parameters will be loaded and used as an initialisation point for the optimiser
  - `init_type`: type of method used to determine initialisation parameters; necessary for correct loading of lr parameters, can be any string otherwise

 `type=lr` Method Parameters:
  - `delta_beta`: slope of the linear ramp used for the mixer
  - `delta_gamma`: slope of the linear ramp used for the phase operator
  - `negative_mixer` whether or not to use a negative ramp for the mixer; evaluates to false when set to anything but `true`

Problem Instance parameters: 
 - `qubits`: problemsize and qubits used for the qaoa circuit. for maxcut,vertexcover and MIS, the problemsize corresponds to the number of nodes in the graph of the problem instance, while for max3sat, the problemsize corresponds to m+n, where m is the number of clauses and n is the number of instances
 - `seed`: seed to use when generating random numbers during problem initialisation, set zo ensure reproducibility
 - `degree`: for maxcut,vertexcover and MIS, the maximum degree/valency of the graph, i.e the maximum number of edges incident to a given vertex
 - `alpha`: for max3sat, the ratio of clauses to variables

##### Troubleshooting
- if the script fails due to certain variables being empty even though they are correctly specified in the config, make sure the config file ends in a newline 
- if `plot_from_config.sh` fails to plot the results in a directory, this may be due to the directory containing data of differing dimensions. make sure the directory of a given method contains results with the same depth, resolution and symmetry, otherwise computing  averages and measures of spread will fail
- if `plot_from_config.sh` fails due to not finding a directory, make sure the this directory is present in results/reproduction before trying to plot the average data. if necessary, call `plot_from_config.sh` after the results have finished.
