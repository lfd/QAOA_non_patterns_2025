#!/bin/bash

function read_config_and_plot_experiments () {
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

    for ((i=0; i<${#methods[@]}; i++)); do
        declare -n met="${methods[i]}"

        cmd="python ./src/main.py -p $depth -pstart $pstart -pstep $stepsize -sy $symmetry -meth ${met[dir_name]} -pltpd -skpv"

        if [[ "$problem" = "max3sat" ]]; then
            res_dir1="./results/reproduction/${problem}/easy/${met[dir_name]}"
            res_dir2="./results/reproduction/${problem}/hard/${met[dir_name]}"
            cmd="$cmd -pltst $res_dir1 2>/dev/null; $cmd -pltst $res_dir2 2>/dev/null"
        else
            res_dir="./results/reproduction/${problem}/${met[dir_name]}"
            cmd="$cmd -pltst $res_dir"
        fi
        eval "$cmd"

    done

}

for arg # for loops through positional input parameters by default, no need to specify in "$@"
do read_config_and_plot_experiments $arg
done

if [[ $# -eq 0 ]] ; then #if no input parameters were given
    echo "no configuration file was found, please pass one or more config files to this script as commandline arguments; see README for more details"
fi