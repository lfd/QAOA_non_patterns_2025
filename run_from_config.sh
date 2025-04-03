#!/bin/bash

function read_config_and_run_experiments () {
    #parse config and create arrays with experiment parameters for each instance
    instances=()
    methods=()
    while read line; do 
        if [[ $line =~ ^"{"(.+)":"(.*)"}" ]]; then
            local ${BASH_REMATCH[1]}="${BASH_REMATCH[2]}"
        elif [[ $line =~ ^"["(.+)"]" ]]; then 
            arrname="${BASH_REMATCH[1]}"
            instances+=($arrname)
            declare -A $arrname
        elif [[ $line =~ ^"("(.+)")" ]]; then 
            arrname="${BASH_REMATCH[1]}"
            methods+=($arrname)
            declare -A $arrname
        elif [[ $line =~ ^(.+)"="(.*) ]]; then 
           declare ${arrname}[${BASH_REMATCH[1]}]="${BASH_REMATCH[2]}"
        fi
    done < $1
    echo "=====Loading Config======"
    echo "problem: $problem"
    echo "symmetry $symmetry"
    echo "starting depth: $pstart"
    echo "target depth: $depth"
    echo "stepsize: $stepsize"
    echo "resolution: $resolution"
    echo "-------------------------"
    #for each problem instance parsed from the config file, run the following
    for ((i=0; i<${#instances[@]}; i++)); do
        declare -n inst="${instances[i]}" #assign nameref so we can call each instance with "inst" instead of by its assigned arrayname
        count=$((i+1))
        echo "=====loading instance $count/${#instances[@]}====="
        echo "qubits: ${inst[qubits]}"
        echo "seed: ${inst[seed]}"
        echo "degree: ${inst[degree]}"
        echo "clause-to-var ratio (alpha): ${inst[alpha]}"

        #determine values depending only on the problemtype
        if [[ "$problem" = "max3sat" ]]; then
            if [[ (( "${inst[alpha]} > 4.9")) || (("${inst[alpha]} =< 3.5")) ]]; then #easy range 
                echo "alpha satisfiability: easy"
            else #hard range
                echo "alpha satisfiability: hard"
            fi
            eigsolvercmd="python ./src/main.py -pr $problem -prs ${inst[seed]} -n ${inst[qubits]} -a ${inst[alpha]} -eig -skp" #set command for eigensolver, we will use this later
        else
            eigsolvercmd="python ./src/main.py -pr $problem -prs ${inst[seed]} -n ${inst[qubits]} -d ${inst[degree]} -eig -skp"
        fi

        #for each instance in this problem instance array, run every method once
        for ((j=0; j<${#methods[@]}; j++)); do
            count=$((j+1))
            declare -n met="${methods[j]}" #assign nameref so we can call each instance with "met" instead of by its assigned arrayname
            echo "Param selection method $count/${#methods[@]}: ${met[dir_name]} ..."
            echo "method type: ${met[type]}"

            # base command, with args that are used by every method
            cmd="python ./src/main.py -pr $problem -p $depth -pstart $pstart -r $resolution -sy $symmetry -meth ${met[dir_name]} -prs ${inst[seed]} -n ${inst[qubits]} -skpv -upmax"
            # determine values depending on both the problemtype and the method
            if [[ "$problem" = "max3sat" ]]; then
                cmd="$cmd -a ${inst[alpha]}" #append problem specific args to base command
                inst_outpath="${resolution}x${resolution}scan_${inst[qubits]}qubits_p${pstart}to${depth}" #name of instance dir that will be generated as a subdir of $outpath and contains results for this instance
                log_outpath="results/reproduction/logs/log_${problem}_n${inst[qubits]}_alpha${inst[alpha]/.}"
            else
                cmd="$cmd -d ${inst[degree]}"
                graph="graph${inst[seed]}"
                inst_outpath="${graph}_${resolution}x${resolution}scan_${inst[qubits]}qubits_p${pstart}to${depth}"
                log_outpath="results/reproduction/logs/log_${problem}_n${inst[qubits]}_graph${inst[seed]}"
            fi

            # if condition for different methods; constructs command accordingly by adding method specific arguments to base command
            if [[ "${met[type]}" =~ "opt" ]]; then #optimised parameters
                if [ ! -z "${met[init]}" ]; then #if startparameters are given
                    echo "startparams: ${met[init]}"
                    echo "startparams type: ${met[init_type]}"
                    if  [[ "${met[init_type]}" =~ "lr" ]]; then
                        cmd="$cmd -pstep $stepsize -opt ${met[optimiser]} -stp ./${met[init]}/${inst_outpath}_None_optimized_${symmetry}"
                    else
                        cmd="$cmd -pstep $stepsize -opt ${met[optimiser]} -stp ./${met[init]}/${inst_outpath}_${symmetry}" #startparameters are associated by instance dir name, which is generated identically for identical instances
                    fi
                else #if startparameters are not given
                    cmd="$cmd -pstep $stepsize -opt ${met[optimiser]}"
                fi
            elif [[ "${met[type]}" =~ "lr" ]]; then #linear ramp parameters
                echo "Delta_beta: ${met[delta_beta]}"
                echo "Delta_beta: ${met[delta_gamma]}"
                echo "mixer_negative: ${met[negative_mixer]}"
                if [[ "${met[negative_mixer]}" =~ "true" ]]; then
                    python -c"import src.problem_util as putil; putil.makeLRschedule($depth, ${met[delta_beta]}, ${met[delta_gamma]}, beta_neg=True)"
                    cmd="$cmd -pstep $stepsize -stp ./fixed_params/linear_ramp_p$depth\_neg_beta.json"
                else
                    python -c"import src.problem_util as putil; putil.makeLRschedule($depth, ${met[delta_beta]}, ${met[delta_gamma]})"
                    cmd="$cmd -pstep $stepsize -stp ./fixed_params/linear_ramp_p$depth.json"
                fi
            fi #sequential method does not need any further arguments

            #execute constructed command
            eval "$cmd >> ${log_outpath}_p${depth}_${met[type]}.out &"

        done

        #run eigensolver once
        echo "Running eigensolver..."
        eval "$eigsolvercmd >> ${log_outpath}_eigensolver.out &"

    done
}

mkdir -p ./results/reproduction/logs

#read specified configs
for arg # for loops through positional input parameters by default, no need to specify in "$@"
do read_config_and_run_experiments $arg
done

if [[ $# -eq 0 ]] ; then #if no input parameters were given
    echo "no configuration file was found, please pass one or more config files to this script as commandline arguments; see README for more details"
fi