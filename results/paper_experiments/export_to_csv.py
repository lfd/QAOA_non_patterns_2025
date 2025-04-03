import json
import csv
from statistics import mean, stdev
import math

import sys, os

sys.path.append('../../src')
from plot_scan import read_multiple_results_and_calc_stats, read_best_param_results, read_fixed_param_results, read_optimised_param_results, read_ground_truth

def export_single_instances(writer, param_writer, files, max_p=7, p_step=1, nesting_level=1, add_offset=False):
        for res in files:
            offset = 0
            if add_offset==True:
                offset_file = res[:-42] + "offset.json"
                with open(offset_file, "r") as f:
                        data = json.load(f)
                        offset = data["offset"]
            
            for i in range(0,max_p,p_step):

                file_name = res.format(i+1)
                problem, method, additional_info = parse_name(file_name,nesting_level=nesting_level)
                qubits = file_name[-13:-11]

                with open(file_name, "r") as f:
                    data = json.load(f)

                for x, row in enumerate(data["scan_grid"]):
                    for y, col in enumerate(row):
                        #if(i == 1):
                            #print(i, x, y, col, data["expecation_values"][x][y])
                        j = i/p_step
                        writer.writerow([problem, method, qubits, i+1, x, y, col[i], col[i*2 + 1], data["expecation_values"][x][y] + offset, additional_info])
            export_param_quality_data(param_writer,res,problem,method,additional_info,max_p,p_step,add_offset=add_offset)

def export_multiple_instances_stats(writer, param_writer, dirs,max_p=7,p_step=1, nesting_level=1,add_offset=False):
    for dir in dirs:
        
        problem, method, additional_info = parse_name(dir,nesting_level=nesting_level)

        avg_mat, std_mat, scan, p, resolution = read_multiple_results_and_calc_stats(dir,"mesh",add_offset,load_energies_from_file=True,pstep=p_step)
        for i in range(0,max_p,p_step):
            j = int(i/p_step)

            for x in range(resolution):
                for y in range(resolution):
                    writer.writerow([problem,f"{method}_average", "multiple", i+1, x, y, scan[x][y][1], scan[x][y][0], avg_mat[j][x][y], additional_info])

            for x in range(resolution):
                for y in range(resolution):
                    writer.writerow([problem, f"{method}_standard_deviation","multiple", i+1, x, y, scan[x][y][1], scan[x][y][0], std_mat[j][x][y], additional_info])
        export_multiple_param_quality_data(param_writer,dir,problem,method,additional_info,max_p,p_step,add_offset=add_offset)

def export_multiple_instances_all_params(param_writer, dirs,max_p=7,p_step=1, nesting_level=1,add_offset=False):
    for path in dirs:
        
        problem, method, additional_info = parse_name(path,nesting_level=nesting_level)
        for res in [path + dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
            print(res)

            export_param_quality_data(param_writer,res,problem,method,additional_info,max_p,p_step,add_offset=add_offset,id_instance=True)

def export_param_quality_data(writer, path, problem, method="sequential", additional_info=None, maxp=7,p_step=1,add_offset=False,id_instance=False):
    #only solution quality of the method used and best parameters that were found, not complete landscape
    if id_instance==False:
        path = path[:-43]
    offset = 0
    if add_offset==True:
        offset_file = f"{path}/offset.json"
        with open(offset_file, "r") as f:
                data = json.load(f)
                offset = data["offset"]
    if additional_info is None:
        if problem=="max3sat":
            return
        else:
            parentdir = f"{problem}/"   
    else:
        parentdir = f"{problem}/{additional_info}/"
    short_path = path[(len(parentdir)+len(method)+1):]
    print(additional_info)
    print(path)

    min_energy, max_energy, qubits = read_ground_truth(short_path,parentdir,offset)
    if method=="sequential" or method=="no_symmetry" or method=="symmetry":
        params, energy_vals, best_found_params, best_found_energy_vals = read_best_param_results(path)
    elif method=="fixed":
        params, energy_vals, best_found_params, best_found_energy_vals = read_fixed_param_results(path,f"fixed_params/fixed_parameters_p{maxp}.json",maxp,p_step)
    elif method=="linear_ramp":
        params, energy_vals, best_found_params, best_found_energy_vals = read_fixed_param_results(path,f"fixed_params/linear_ramp_p{maxp}.json",maxp,p_step)
    elif method=="linear_ramp_neg_mixer":
        params, energy_vals, best_found_params, best_found_energy_vals = read_fixed_param_results(path,f"fixed_params/linear_ramp_p{maxp}_neg_beta.json",maxp,p_step)
    elif method=="optimised":
        params, energy_vals, best_found_params, best_found_energy_vals = read_optimised_param_results(path)
    else:
        print("this method has not been implemented yet")
    
    approx = [((e+offset)-min_energy)/(max_energy-min_energy) for e in energy_vals]
    best_found_approx = [((e+offset)-min_energy)/(max_energy-min_energy) for e in best_found_energy_vals]
    if method=="sequential":
        p_step=1 # override pstep for sequential method, which always requires pstep of 1
    for i in range(0,maxp,p_step):
        j = int(i/p_step)
        if id_instance==False:
            writer.writerow([problem, method, qubits, i+1, params[i], params[maxp+i], approx[j], energy_vals[j]+offset, best_found_params[j][:i+1], best_found_params[j][i+1:], best_found_approx[j], best_found_energy_vals[j]+offset, additional_info])
        else:
            if additional_info==None:
                string=short_path
            else:
                string = short_path + additional_info
            writer.writerow([problem, method, qubits, i+1, params[i], params[maxp+i], approx[j], energy_vals[j]+offset, best_found_params[j][:i+1], best_found_params[j][i+1:], best_found_approx[j], best_found_energy_vals[j]+offset, string])

def export_multiple_param_quality_data(writer, path, problem, method="sequential", additional_info=None, maxp=7,p_step=1,add_offset=False):
    energies = []
    best_found_energies = []
    approx = []
    best_found_approx = []
    for dir_name in [path + dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
        offset = 0
        if add_offset==True:
            offset_file = f"{dir_name}/offset.json"
            with open(offset_file, "r") as f:
                    data = json.load(f)
                    offset = data["offset"]

        if additional_info is None:
            parentdir = f"{problem}/"
        else:
            parentdir = f"{problem}/{additional_info}/"
        short_path = dir_name[(len(parentdir)+len(method)+1):]
        min_energy, max_energy, qubits = read_ground_truth(short_path,parentdir,offset)
        
        if method=="sequential":
            params, energy, best, best_energy = read_best_param_results(dir_name)
        elif method=="fixed":
            params, energy, best, best_energy = read_fixed_param_results(dir_name,f"fixed_params/fixed_parameters_p{maxp}.json",maxp,p_step)
        elif method=="linear_ramp":
            params, energy, best, best_energy = read_fixed_param_results(dir_name,f"fixed_params/linear_ramp_p{maxp}.json",maxp,p_step)
        elif method=="linear_ramp_neg_mixer":
            params, energy, best, best_energy = read_fixed_param_results(dir_name,f"fixed_params/linear_ramp_p{maxp}_neg_beta.json",maxp,p_step)
        elif method=="optimised":
            params, energy, best, best_energy = read_optimised_param_results(dir_name)
        
        energy_approx = [((e+offset)-min_energy)/(max_energy-min_energy) for e in energy]
        best_energy_approx = [((e+offset)-min_energy)/(max_energy-min_energy) for e in best_energy]

        energy_offset_list = [offset]*len(energy)
        energies.append(energy+energy_offset_list)
        approx.append(energy_approx)
        best_energy_offset_list = [offset]*len(best_energy)
        best_found_energies.append(best_energy + best_energy_offset_list)
        best_found_approx.append(best_energy_approx)

    for i in range(0,maxp,p_step):
        j = int(i/p_step)
        approx_current_p = [item[j] for item in approx]
        energy_current_p = [item[j] for item in energies]
        best_approx_current_p = [item[j] for item in best_found_approx]
        best_energy_current_p = [item[j] for item in best_found_energies]
        #calc mean
        mean_approx = mean(approx_current_p)
        mean_energy = mean(energy_current_p)
        mean_best_approx = mean(best_approx_current_p)
        mean_best_energy = mean(best_energy_current_p)
        writer.writerow([problem, f"{method}_average", "multiple", i+1, None, None, mean_approx, mean_energy, None, None, mean_best_approx, mean_best_energy, additional_info])
        #calc stdev
        std_approx = stdev(approx_current_p)
        std_energy = stdev(energy_current_p)
        std_best_approx = stdev(best_approx_current_p)
        std_best_energy = stdev(best_energy_current_p)
        writer.writerow([problem, f"{method}_standard_deviation", "multiple", i+1, None, None, std_approx, std_energy, None, None, std_best_approx, std_best_energy, additional_info])

def parse_name(dirname,nesting_level=1):

    problem = dirname[:dirname.find("/")] # first directory in path
    # the following 9 line of code are not pretty, but they get the job done
    if nesting_level==1:
        method = dirname[dirname.find("/")+1:dirname.find("/",dirname.find("/")+1)] # second directory in path
        additional_distinction = None
    if nesting_level==2:
        additional_distinction = dirname[dirname.find("/")+1:dirname.find("/",dirname.find("/")+1)] # second directory in path
        temp = dirname[(len(problem)+len(additional_distinction)+2):] #slice first two directories from path so code does not get too verbose
        method = temp[:temp.find("/")] # third directory in path
    elif nesting_level>2:
        print("WARNING: currently only nesting levels 1 and 2 are supported. Try restructuring your result directories.")

    return problem, method, additional_distinction

if __name__ == "__main__":

    sys.path.append('../../results')

    csv_f = open("qaoa_data.csv", "w", newline='')
    writer = csv.writer(csv_f)
    writer.writerow(["problem", "method", "num_qubits", "p", "x", "y", "gamma", "beta", "expectation", "additional_info"])

    csv_f2 = open("param_data.csv", "w", newline='')
    writer2 = csv.writer(csv_f2)
    writer2.writerow(["problem", "method", "num_qubits", "p", "gamma", "beta", "approx", "expectation", "best_found_gamma", "best_found_beta", "best_found_approx", "best_found_expectation", "additional_info"])

    csv_f3 = open("all_param_data.csv", "w", newline='')
    writer3 = csv.writer(csv_f3)
    writer3.writerow(["problem", "method", "num_qubits", "p", "gamma", "beta", "approx", "expectation", "best_found_gamma", "best_found_beta", "best_found_approx", "best_found_expectation", "additional_info"])


    # single instance plots
    result_files = {
        "max3sat/no_symmetry/32x32scan_16qubits_p1to7_no_symmetry_hard/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/no_symmetry/graph123456_32x32scan_16qubits_p1to7_no_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "max3sat/symmetry/32x32scan_16qubits_p1to4_general_symmetry_hard/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/symmetry/graph123456_32x32scan_16qubits_p1to4_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/sequential/graph123456_32x32scan_16qubits_p1to7_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/optimised/graph123456_32x32scan_16qubits_p1to7_COBYLA_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/fixed/graph123456_32x32scan_16qubits_p1to7_None_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/linear_ramp/graph123456_32x32scan_16qubits_p1to7_None_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "maxcut/linear_ramp_neg_mixer/graph123456_32x32scan_16qubits_p1to7_None_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "vertcover/sequential/graph123456_32x32scan_16qubits_p1to7_general_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "vertcover/linear_ramp/graph123456_32x32scan_16qubits_p1to7_None_optimized_general_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "vertcover/linear_ramp_neg_mixer/graph123456_32x32scan_16qubits_p1to7_None_optimized_general_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "vertcover/optimised/graph123456_32x32scan_16qubits_p1to7_COBYLA_optimized_general_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
    }
    # single instance plots in nested directorys with more than one level of subdirectories
    nested_result_files = {
        "max3sat/easy/sequential/32x32scan_21qubits_p1to7_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/easy/optimised/32x32scan_21qubits_p1to7_COBYLA_optimized_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/easy/linear_ramp/32x32scan_21qubits_p1to7_None_optimized_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/easy/linear_ramp_neg_mixer/32x32scan_21qubits_p1to7_None_optimized_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/hard/sequential/32x32scan_21qubits_p1to7_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/hard/optimised/32x32scan_21qubits_p1to7_COBYLA_optimized_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/hard/linear_ramp/32x32scan_21qubits_p1to7_None_optimized_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
        "max3sat/hard/linear_ramp_neg_mixer/32x32scan_21qubits_p1to7_None_optimized_general_symmetry/estimator_data_32x32scan_p{}_21qubits.json",
    }
    #plots which are average/standard deviation on multiple instances:
    results_dirs = {
        "maxcut/sequential/",
        "maxcut/optimised/",
        "maxcut/fixed/",
        "maxcut/linear_ramp/",
        "maxcut/linear_ramp_neg_mixer/",
        "vertcover/sequential/",
        "vertcover/optimised/",
        "vertcover/linear_ramp/",
        "vertcover/linear_ramp_neg_mixer/",
    }
    #plots which are average/standard deviation on multiple instances in nested directorys with more than one level of subdirectories
    nested_results_dirs = {
        "max3sat/easy/sequential/",
        "max3sat/easy/optimised/",
        "max3sat/easy/linear_ramp/",
        "max3sat/easy/linear_ramp_neg_mixer/",
        "max3sat/hard/sequential/",
        "max3sat/hard/optimised/",
        "max3sat/hard/linear_ramp/",
        "max3sat/hard/linear_ramp_neg_mixer/",
    }

    writer3_dirs = {
        "maxcut/sequential/",
        "maxcut/optimised/",
        "vertcover/sequential/",
        "vertcover/optimised/",
    }

    writer3_higher_depth = {
        "higher_depth_maxcut/optimised/",
        "higher_depth_maxcut/sequential/",
    }

    nested_writer3_dirs = {
        "max3sat/easy/sequential/",
        "max3sat/easy/optimised/",
        "max3sat/hard/sequential/",
        "max3sat/hard/optimised/",
    }

    #plots with p_step=2
    twostep_result_files = {
        "higher_depth_maxcut/linear_ramp/graph123456_32x32scan_16qubits_p1to21_None_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "higher_depth_maxcut/linear_ramp_neg_mixer/graph123456_32x32scan_16qubits_p1to21_None_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "higher_depth_maxcut/optimised/graph123456_32x32scan_16qubits_p1to21_COBYLA_optimized_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
        "higher_depth_maxcut/sequential/graph123456_32x32scan_16qubits_p1to21_maxcut_symmetry/estimator_data_32x32scan_p{}_16qubits.json",
    }
    twostep_results_dirs = {
        "higher_depth_maxcut/linear_ramp/",
        "higher_depth_maxcut/linear_ramp_neg_mixer/",
        "higher_depth_maxcut/optimised/",
        "higher_depth_maxcut/sequential/",
    }

    print("---exporting single instance plots-----")
    export_single_instances(writer,writer2,result_files,add_offset=True)
    export_single_instances(writer,writer2,nested_result_files,nesting_level=2,add_offset=True)
    export_single_instances(writer,writer2,twostep_result_files,max_p=21,p_step=2,add_offset=True)

    print("---exporting multiple instace plots-----")
    export_multiple_instances_stats(writer,writer2,results_dirs,add_offset=True)
    export_multiple_instances_stats(writer,writer2,nested_results_dirs,nesting_level=2,add_offset=True)
    export_multiple_instances_stats(writer,writer2,twostep_results_dirs,max_p=21,p_step=2,add_offset=True)

    print("---exporting params-----")
    export_multiple_instances_all_params(writer3,writer3_dirs,add_offset=True)
    export_multiple_instances_all_params(writer3,writer3_higher_depth,max_p=21,p_step=2,add_offset=True)
    export_multiple_instances_all_params(writer3,nested_writer3_dirs,nesting_level=2,add_offset=True)
    
   # csv_f.close()
    #csv_f2.close()
    csv_f3.close()

