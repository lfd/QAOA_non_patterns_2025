from mpl_toolkits.mplot3d import Axes3D, art3d
from matplotlib.ticker import FuncFormatter, MultipleLocator, MaxNLocator
import matplotlib.pyplot as plt

import os, json, re
import numpy as np
from statistics import mean, stdev
from textwrap import wrap
from operator import sub

def read_existing_result_data(path_to_json,plot_type="voxel",add_offset=False):
    json_res = []
    scan_grid = None
    offset = 0

    if add_offset==True:
      with open(path_to_json + "/offset.json") as offset_file:
          data = json.load(offset_file)
          offset = data["offset"]

    for file_name in [file for file in os.listdir(path_to_json) if file.endswith('.json')]:
      with open(path_to_json + "/" + file_name) as json_file:
        if file_name!="offset.json" and file_name!="eigenvalues.json" and file_name!="opt_data.json":
          data = json.load(json_file)
          if scan_grid is None:
            scan_grid = data["scan_grid"]
          json_res.append(data["expecation_values"])
          best_params = data["best"]["params"]

    p = len(json_res)
    resolution = int(len(json_res[0]))
    mat = np.ndarray((p,resolution,resolution))
    scan = np.ndarray((p,resolution,resolution))
    scan_grid = np.array([p[0] for p in scan_grid])

    for i in range(p):
        expectation_mat = np.asmatrix(json_res[i])
        if plot_type=="voxel":
          mat[i] = np.rot90(expectation_mat + offset)
        else:
          mat[i] = expectation_mat + offset
        for m, row in enumerate(scan_grid):
          for n, scalar in enumerate(row):
            scan[i][m][n] = scalar

    return mat, scan, p, resolution, best_params

def read_multiple_results_and_calc_stats(path_to_dirs,plot_type="voxel",add_offset=False,load_energies_from_file=False,pstep=1):
    print(f"reading results from {path_to_dirs} ...")
    json_dict = {}
    json_res = []
    json_res_list = []
    scan_grid = None
    scan = None

    for dir_name in [path_to_dirs + "/" + dir for dir in os.listdir(path_to_dirs) if os.path.isdir(path_to_dirs + "/" + dir)]:
      if dir_name.endswith("symmetry"):

        offset = 0
        if add_offset==True:
          with open(dir_name+ "/offset.json") as offset_file:
              data = json.load(offset_file)
              offset = data["offset"]

        for file_name in [file for file in os.listdir(dir_name) if file.endswith('.json')]:
          with open(dir_name + "/" + file_name) as json_file:
            data = json.load(json_file)
            if file_name=="eigenvalues.json":
              continue
              #eigenvals_found = True
              #min_energy = data["min_energy"]
              #max_energy =data["max_energy"]
            elif file_name=="opt_data.json" or file_name=="offset.json":
              continue
            else:
              file_grid = data["scan_grid"]
              file_p = len(file_grid[0][0])/2
              if file_p==1: # basically if p=1
                scan_grid = file_grid
              if pstep==1 or file_p%pstep == 1: #assumes starting from p=1, we want file_p%pstep to be equal to p_start%pstep (pstart is not passed to the function currently)
                json_dict[file_p] = data["expecation_values"]
        
        for i in sorted(json_dict.keys()):
          json_res.append(json_dict[i])

        p = len(json_res)
        resolution = int(len(json_res[0]))
        mat = np.ndarray((p,resolution,resolution))
        for i in range(p):
          expectation_mat = np.asmatrix(json_res[i])
          if plot_type=="voxel":
            mat[i] = np.rot90(expectation_mat + offset)
          else:
            mat[i] = expectation_mat + offset

        if load_energies_from_file:
          split_dirs = dir_name.split("/")
          problem = split_dirs[0]

          if problem=="max3sat":
            problem = f"{split_dirs[0]}/{split_dirs[1]}"
          subdir = split_dirs[-1]

          if subdir[0:5]=="graph":
            qubits = subdir[22:24]
            seed = subdir[5:11]
          else:
            qubits = subdir[10:12]
            seed = 123456
          #load min and max val
          min_path = f"{problem}/eigensolver/min_eigenvalue_{seed}_{qubits}qubits.json"
          max_path = f"{problem}/max_energies/n{qubits}_seed{seed}.json"
          with open(min_path) as f:
            data = json.load(f)
            x_min = data["min_energy"]+offset
          with open(max_path) as f:
            data = json.load(f)
            x_max = data["max_energy"]+offset

        else:
          # normalize data with best found values
          x_min = mat.min(axis=(1, 2), keepdims=True)
          x_max = mat.max(axis=(1, 2), keepdims=True) #offset already added in this case
        mat_norm = (mat-x_min)/(x_max-x_min)

        json_res_list.append(mat_norm)
        json_res.clear()

        if scan is None and plot_type=="voxel":
          scan = np.ndarray((p,resolution,resolution))
          scan_grid = np.array([p[0] for p in scan_grid])
          for i in range(p):
            for m, row in enumerate(scan_grid):
              for n, scalar in enumerate(row):
                scan[i][m][n] = scalar
        elif scan is None:
          scan = scan_grid
    
    if not all(x.shape[0]==json_res_list[0].shape[0] for x in json_res_list):
       raise Exception(f"Error: cannot average scans due to differing dimensions: ensure that all results in {path_to_dirs} have the same resolution, depth and symmetry")

    avg_mat = sum(json_res_list)/len(json_res_list)

    for i, json_res in enumerate(json_res_list):
      json_res_list[i] = (json_res - avg_mat)**2
    std_mat = np.sqrt(sum(json_res_list)/len(json_res_list))

    return avg_mat, std_mat, scan, p, resolution

def save_and_close(path_to_data, plot_type_str, resolution, p, n, affix="", format=["png"],pstep=1,pstart=1):
    print(f"saving plot: creating file for each format in {format}...")
    for ext in format:
      if resolution is not None:
        fname = f"{path_to_data}/{plot_type_str}plot_{resolution}x{resolution}scan_p{p*pstep}_{n}qubits{affix}.{ext}"
      else:
        fname = f"{path_to_data}/{plot_type_str}plot_p{p}_{n}qubits{affix}.{ext}"
      plt.savefig(fname)
      print(f"saved as {fname}")
    plt.close()

def meshplot_series(path_to_json,beta_bounds,gamma_bounds,n,cmap_theme="viridis",title=None,plot_avg=False,plot_std=False,format=["png"],pstep=1,pstart=1):  
    '''
    path_to_json: path to experiment data
    beta_bounds: beta grid bounds
    gamma_bounds: gamma grid bounds
    n: problem oder/number of qubits
    cmap_theme: matplotlib colormap theme for the plot
    title: plot suptitle
    '''
    if plot_avg:
      savefile_affix="_norm_avg"
      cbar_label=r"Average residual Energy $\frac{F(|\beta,\gamma\rangle) - E_0}{E_{max} - E_0}$"
      mat, std_mat, scan, p, resolution = read_multiple_results_and_calc_stats(path_to_json,"mesh")
    elif plot_std:
      savefile_affix="_norm_std"
      cbar_label=r"Standard deviation of residual Energy $\frac{F(|\beta,\gamma\rangle) - E_0}{E_{max} - E_0}$"
      avg_mat, mat, scan, p, resolution = read_multiple_results_and_calc_stats(path_to_json,"mesh")
    else:
      savefile_affix=""
      cbar_label=r"Expectation value F($|\beta$,$\gamma\rangle$)"
      mat, scan, p, resolution, best_params = read_existing_result_data(path_to_json,"mesh")
    print("Plotting results as colormesh...")
    
    beta_indices = np.flip(np.linspace(beta_bounds[0],beta_bounds[1],resolution))
    gamma_indices = np.flip(np.linspace(gamma_bounds[0],gamma_bounds[1],resolution))#
    num_rows = int(np.ceil(p/4))
    num_cols = min(p, 4)
    fig, axs = plt.subplots(num_rows, num_cols, sharey=True, sharex=True)
    if num_rows > 1:
      ax = axs[0][0]
    elif p==1:
      ax = axs
    else:
      ax = axs[0]
    ax.tick_params(axis='y', which='major', labelsize=7)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda val,pos: '{:.2g}$\pi$'.format(val/np.pi) if val !=0 else '0'))
    ax.yaxis.set_major_locator(MultipleLocator(base=np.pi/4))

    for i, item in enumerate(mat):
      if num_rows > 1:
        ax = axs[int(np.floor(i/4))][i%4]
      elif p==1:
        ax = axs
      else:
        ax = axs[i] 
      ax.set_box_aspect(1)
      ax.pcolormesh(gamma_indices, beta_indices, item, vmin=mat.min(), vmax=mat.max(), cmap=cmap_theme, shading="auto")
      if i==0:
        ax.set_title(f"p={pstart}",fontsize=7)
      ax.set_title(f"p={pstep*i+pstart}",fontsize=7)
      if int(np.ceil((i+1)/4))==num_rows:
        ax.set_xlabel(r'$\gamma$')
      if i%4 == 0:
        ax.set_ylabel(r'$\beta$')
      ax.xaxis.set_major_formatter(FuncFormatter(lambda val,pos: '{:.2g}$\pi$'.format(val/np.pi) if val !=0 else '0'))
      ax.xaxis.set_major_locator(MultipleLocator(base=np.pi/2))
      ax.tick_params(axis='x', which='major', labelsize=7, rotation=30)
      ax.tick_params(axis='y', which='major', labelsize=7)

    for i in range(0,num_rows*num_cols-p):
      ax = axs[num_rows-1][num_cols-(i+1)]
      ax.axis("off")
      ax2 = axs[num_rows-2][num_cols-(i+1)]
      ax2.set_xlabel(r'$\gamma$')
      ax2.xaxis.set_tick_params(labelbottom=True)

    fig.subplots_adjust(wspace=0.3,top=0.85,bottom=0.35)

    cmap = plt.get_cmap(cmap_theme)
    norm = plt.Normalize(mat.min(), mat.max())

    if "tex" in format:
      print(f"colorbar - min:{mat.min()},max:{mat.max()}")
    if num_rows > 1:
      p0 = axs[num_rows-1][0].get_position().get_points().flatten()
      p1 = axs[num_rows-1][num_cols-1].get_position().get_points().flatten()
      ax_cbar = fig.add_axes([p0[0], 0.55-num_rows*0.2, p1[2]-p0[0], 0.04])
    elif p==1:
      p0 = axs.get_position().get_points().flatten()
      p1 = axs.get_position().get_points().flatten()
      ax_cbar = fig.add_axes([p0[0], 0.2, p1[2]-p0[0], 0.04])
    else:
      p0 = axs[num_rows-1].get_position().get_points().flatten()
      p1 = axs[num_cols-1].get_position().get_points().flatten()
      ax_cbar = fig.add_axes([p0[0], 0.2, p1[2]-p0[0], 0.04])
  
    m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    m.set_array([])

    plt.colorbar(m, cax=ax_cbar, orientation="horizontal", label=cbar_label)

    if not (plot_avg or plot_std or num_rows>1):
      subtitle = "\n".join(wrap(r"${\bf{Best \: \beta}}$: "+str(best_params[0][:p])+r" ${\bf{Best \: \gamma}}$: "+str(best_params[0][p:]),70))
      plt.figtext(0.1, 0, subtitle, fontsize=8)

    if title is None:
      fig.suptitle("\n".join(wrap("Expectation values for QAOA with increasing p for a "+str(n)+" qubit problem instance",70)))
    else:
      fig.suptitle("\n".join(wrap(title,70)))

    save_and_close(path_to_json,"mesh",resolution,p,n,savefile_affix,format,pstep,pstart)

def voxelplot(path_to_json,beta_bounds,gamma_bounds,n,mask_vals_leq=None,mask_vals_geq=0,elevation=25,azimuth=-120,cmap_theme="viridis",title=None,plot_avg=False,plot_std=False,format=["png"],plot_best=False):
  '''
  path_to_json: path to experiment data
  beta_bounds: beta grid bounds
  gamma_bounds: gamma grid bounds
  n: problem oder/number of qubits
  mask_vals_geq: mask values greater or equal to this in second plot
  mask_vals_leq: mask values less or equal to this in second plot
  elevation: voxelplot elevation angle 
  azimuth: voxelplot azimuthal angle
  cmap_theme: matplotlib colormap theme for the plot
  title: plot suptitle
  '''
  if plot_avg:
    savefile_affix="_norm_avg"
    cbar_label=r"Average residual Energy $\frac{F(|\beta,\gamma\rangle) - E_0}{E_{max} - E_0}$"
    mat, std_mat, scan, p, resolution = read_multiple_results_and_calc_stats(path_to_json,"voxel")
  elif plot_std:
    savefile_affix="_norm_std"
    cbar_label=r"Standard deviation of residual Energy $\frac{F(|\beta,\gamma\rangle) - E_0}{E_{max} - E_0}$"
    avg_mat, mat, scan, p, resolution = read_multiple_results_and_calc_stats(path_to_json,"voxel")
  else:
    savefile_affix=""
    cbar_label=r"Expectation value F($|\beta$,$\gamma\rangle$)"
    mat, scan, p, resolution, best_params = read_existing_result_data(path_to_json,"voxel")
  print("Plotting results as voxelplot...")
  scan2 = scan
  if mask_vals_geq:
    scan2 = np.ma.masked_where(mat >= mask_vals_geq, scan2)
  if mask_vals_leq:
    scan2 = np.ma.masked_where(mat <= mask_vals_leq, scan2)

  cmap = plt.get_cmap(cmap_theme)
  norm = plt.Normalize(mat.min(), mat.max())

  fig, (ax0, ax1) = plt.subplots(1, 2, subplot_kw=dict(projection='3d'))
  ax0.set_aspect('auto')
  ax0.voxels(scan, facecolors=cmap(norm(mat)))
  ax1.set_aspect('auto')
  if plot_best:
    ax1.voxels(scan2, facecolors=cmap(norm(mat)),alpha=0.8)
  else:
    ax1.voxels(scan2, facecolors=cmap(norm(mat)))

  beta_ticks = ['{:.2f}$\pi$'.format(x) for x in np.linspace(beta_bounds[0]/np.pi,beta_bounds[1]/np.pi,5)]
  gamma_ticks = ['{:.2f}$\pi$'.format(x) for x in np.flip(np.linspace(gamma_bounds[0]/np.pi,gamma_bounds[1]/np.pi,5))]

  for ax in (ax0, ax1):
    ax.view_init(elev=elevation, azim=azimuth)
    ax.set_xlabel("p")
    ax.set_ylabel(r'$\gamma$')
    ax.set_zlabel(r'$\beta$')
    ax.xaxis.set_major_locator(MultipleLocator(base=1))
    ax.yaxis.set_major_locator(MultipleLocator(base=resolution/4))
    ax.zaxis.set_major_locator(MultipleLocator(base=resolution/4))
    ax.yaxis.set_ticks(ticks=ax.get_yticks()[1:-1],labels=gamma_ticks,va="center",fontsize=8)
    ax.zaxis.set_ticks(ticks=ax.get_zticks()[1:-1],labels=beta_ticks, va="center",fontsize=8)

  if plot_best==True:
    beta = [resolution-(x/(beta_bounds[1]-beta_bounds[0]))*resolution for x in list(map(sub,best_params[0][:p],[beta_bounds[0]]*p))]
    gamma = [resolution-(x/(gamma_bounds[1]-gamma_bounds[0]))*resolution for x in list(map(sub,best_params[0][p:],[gamma_bounds[0]]*p))]
    filled=np.zeros((p,resolution,resolution))
    for i in range(p):
      filled[i][int(gamma[i])][int(beta[i])] = 1
      #the following 4 lines look wierd, but are here to avoid wrong plotting incase optimizer finds values out of bounds
      if gamma[i] < 0:
        gamma[i] = resolution+gamma[i]%resolution
      if beta[i] < 0:
        beta[i] = resolution+beta[i]%resolution
    ax1.voxels(filled=filled,facecolors="magenta")
    #ax1.plot(np.arange(0.5,p+0.5,1),gamma,beta,linestyle="--",color="magenta")

  fig.subplots_adjust(wspace=0.25,top=0.95,bottom=0.35)

  p0 = ax0.get_position().get_points().flatten()
  p1 = ax1.get_position().get_points().flatten()
  ax_cbar = fig.add_axes([p0[0], 0.3, p1[2]-p0[0], 0.04])
  m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  m.set_array([])
  plt.colorbar(m, cax=ax_cbar, orientation="horizontal", label=cbar_label)

  if not (plot_avg or plot_std):
    subtitle = "\n".join(wrap(r"${\bf{Best \: \beta}}$: "+str(best_params[0][:p])+r" ${\bf{Best \: \gamma}}$: "+str(best_params[0][p:]),95))
    plt.figtext(0.1, 0.05, subtitle, fontsize=8)

  if title is None:
      fig.suptitle("\n".join(wrap("Expectation values for QAOA with increasing p for a "+str(n)+" qubit problem instance",70)))
  else:
      fig.suptitle("\n".join(wrap(title,70)))

  save_and_close(path_to_json,"voxel",resolution,p,n,savefile_affix,format)

def fourier_landscape_series(path_to_json, n, max_freq_u=6, max_freq_v=12, cmap_theme="viridis",title=None,format=["png"]):
  mat, scan, p, resolution, best_params = read_existing_result_data(path_to_json,"fourier")
  halfway = int(resolution/2)
  max_x = int(halfway + max_freq_v)
  min_y = int(halfway - max_freq_u)
  max_y = int(halfway + max_freq_u)

  print("transforming parameters to frequency domain...")
  fmat = np.ndarray((p,max_freq_u*2,max_freq_v))
  for i in range(p):
    transformed = np.fft.fftshift(np.fft.fftn(mat[i], norm="forward"))
    fmat[i] = np.abs(transformed)[min_y:max_y, halfway:max_x]

  x_axis = np.arange(0, max_freq_v)
  y_axis = np.arange(-max_freq_u, max_freq_u)
  X, Y = np.meshgrid(x_axis, y_axis)

  fig, axs = plt.subplots(1, p, sharey=True)
  ax = axs[0]
  ax.set_ylabel(r'$\vec{u}$')
  ax.tick_params(axis='y', which='major', labelsize=7)

  for i, item in enumerate(fmat):
    ax = axs[i]
    ax.set_box_aspect(1)
    ax.pcolormesh(X,Y,item, cmap=cmap_theme, shading="nearest")
    ax.set_title(f"p={i+1}",fontsize=7)
    ax.set_xlabel(r'$\vec{v}$')
    ax.tick_params(axis='x', which='major', labelsize=7)

  fig.subplots_adjust(wspace=0.3,top=1.0,bottom=0.35)

  cmap = plt.get_cmap(cmap_theme)
  norm = plt.Normalize(fmat.min(), fmat.max())

  p0 = axs[0].get_position().get_points().flatten()
  p1 = axs[p-1].get_position().get_points().flatten()
  ax_cbar = fig.add_axes([p0[0], 0.4, p1[2]-p0[0], 0.04])
  m = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
  m.set_array([])

  plt.colorbar(m, cax=ax_cbar, orientation="horizontal", label=r"Expectation value F($\vec{u}$,$\vec{v}$)")
  
  if title is None:
      fig.suptitle("\n".join(wrap(r"Fourier Landscapes of QAOA with increasing p for a "+str(n)+" qubit problem instance",70)))
  else:
      fig.suptitle("\n".join(wrap(title,70)))

  save_and_close(path_to_json,"fourier_landscape_",resolution,p,n,"",format)

def lineplot(path,problem,method,maxp,p_start,p_step,add_offset=True,format=["png"]):
  energies = []
  approx = []
  params_arr = []
  p = []
  avg = []
  std = []
  for dir_name in [path + "/" + dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
      offset = 0
      if add_offset==True:
          offset_file = f"{dir_name}/offset.json"
          with open(offset_file, "r") as f:
                  data = json.load(f)
                  offset = data["offset"]

      parentdir = path[:-(len(method))]
      short_path = dir_name[(len(parentdir)+len(method)+1):]
      min_energy, max_energy, qubits = read_ground_truth(short_path,parentdir,offset)
      params, energy, best, best_energy = read_params(method,dir_name,maxp,p_step)

      energy_approx = [((e+offset)-min_energy)/(max_energy-min_energy) for e in energy]

      energy_offset_list = [offset]*len(energy)
      energies.append(energy+energy_offset_list)
      approx.append(energy_approx)
      params_arr.append(params)

  for i in range(0,maxp,p_step):
    j = int(i/p_step)
    approx_current_p = [item[j] for item in approx]
    #calc mean
    mean_approx = mean(approx_current_p)
    #calc stdev
    std_approx = stdev(approx_current_p)
    p.append(p_start+i)
    avg.append(mean_approx)
    std.append(std_approx)

  print("---Lineplot (1/2)---")

  fig, ax = plt.subplots()
  ax.errorbar(p,avg,yerr=std,capsize=5,marker='o')
  ax.set_ylabel(r"average residual Energy $\frac{F(|\beta,\gamma\rangle) - E_0}{E_{max} - E_0}$")
  ax.set_xlabel(r"depth $p$")
  ax.xaxis.set_major_locator(MaxNLocator(integer=True))
  fig.suptitle(f"Parameter quality for {method} applied to {problem}")

  save_and_close(path,"param_avg_quality_",None,maxp,"unspecified","",format)

  beta = params[:maxp]
  gamma = params[maxp:]

  print("---Lineplot (2/2)---")

  fig, (ax0, ax1) = plt.subplots(1, 2)

  for params in params_arr:
    beta = params[:maxp]
    gamma = params[maxp:]
    ax0.plot(range(maxp),beta,marker='o')
    ax1.plot(range(maxp),gamma,marker='o')
  ax0.set_ylabel(r'$\beta$')
  ax1.set_ylabel(r'$\gamma$')
  ax0.xaxis.set_major_locator(MaxNLocator(integer=True))
  ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
  fig.supxlabel(r"depth $p$")
  fig.suptitle(f"Parameter arrangement for {method} applied to {problem}")

  save_and_close(path,"param_arrangement_",None,maxp,"unspecified","",format)

def read_params(method,dir_name,maxp,p_step):
  if "/" in method:
    method = method.split("/")[1]
  if method=="fixed":
      params, energy, best, best_energy = read_fixed_param_results(dir_name,f"fixed_params/fixed_parameters_p{maxp}.json",maxp,p_step)
  elif method=="linear_ramp" or method=="linear_ramp_pos_mixer":
      params, energy, best, best_energy = read_fixed_param_results(dir_name,f"fixed_params/linear_ramp_p{maxp}.json",maxp,p_step)
  elif method=="linear_ramp_neg_mixer":
      params, energy, best, best_energy = read_fixed_param_results(dir_name,f"fixed_params/linear_ramp_p{maxp}_neg_beta.json",maxp,p_step)
  elif method=="optimised":
      params, energy, best, best_energy = read_optimised_param_results(dir_name)
  else:
    params, energy, best, best_energy = read_best_param_results(dir_name)

  return params, energy, best, best_energy

def read_ground_truth(respath,parentdir,offset=0):
  if respath[:5]=="graph":
    seed =respath[5:11]
    qubits = respath[22:24]
    if qubits.isdigit() is False:
      qubits = qubits[:-1]
  else:
    seed = 123456
    qubits = respath[10:12]
    if qubits.isdigit() is False:
      qubits = qubits[:-1]
  min_path = f"{parentdir}eigensolver/min_eigenvalue_{seed}_{qubits}qubits.json"
  max_path = f"{parentdir}max_energies/n{qubits}_seed{seed}.json"
  with open(min_path) as f:
    data = json.load(f)
    min = data["min_energy"]+offset
  with open(max_path) as f:
    data = json.load(f)
    max = data["max_energy"]+offset
  return min, max, qubits

def read_best_param_results(path):
  energy_vals = []
  best_params = []
  for file_name in sorted([file for file in os.listdir(path) if (file.endswith('.json') and file!="eigenvalues.json" and file!="offset.json")], key=lambda x: int(x.split("scan_p")[1].split("_")[0])):
    with open(path + "/" + file_name) as json_file:
      data = json.load(json_file)
      best = data["best"]["params"][0]
      expectation = data["best"]["expectation"][0]
      energy_vals.append(expectation)
      best_params.append(best)
  return best, energy_vals, best_params, energy_vals #best_params, energy_vals returned twice since the sequential parameters and the best found parameters are identical

def read_fixed_param_results(path, path_to_fixed, maxp = 7, pstep=1):
  with open(path_to_fixed) as json_file:
    fixed_data = json.load(json_file)
  init_params = fixed_data['opt_params']
  energy_vals = []
  best_energy_vals = []
  best_params = []
  p = 0
  for file_name in sorted([file for file in os.listdir(path) if (file.endswith('.json') and file!="eigenvalues.json" and file!="offset.json")], key=lambda x: int(x.split("scan_p")[1].split("_")[0])):
    p += pstep
    with open(path + "/" + file_name) as json_file:
      data = json.load(json_file)
      #best found energies
      best = data["best"]["params"][0]
      expectation = data["best"]["expectation"][0]
      best_params.append(best)
      best_energy_vals.append(expectation)

      #fixed param energies
      beta = init_params[:p]
      gamma = init_params[maxp:(maxp+p)]
      current_p_params = beta + gamma
      i,j = find_nearest_in_scan_grid(data["scan_grid"],current_p_params)
      fixed_energy = data["expecation_values"][i][j]
      energy_vals.append(fixed_energy)

  return init_params, energy_vals, best_params, best_energy_vals

def read_optimised_param_results(path):
  opt_datapath = path + "/opt_data.json"
  if os.path.exists(opt_datapath):
    with open(opt_datapath) as json_file:
      data = json.load(json_file)
      return data["params"], data["expectation_values"], data["best"], data["best_expectation_values"]
  else:
    raise Exception("Sorry, opt_data.json must be present in the directory. use the method parse_optimised_param_results(path, path_to_log, maxp = 7) to create such a file ")

def parse_optimised_param_results(path, path_to_log, maxp = 7, pstep=1):
  opt_params_pattern = re.compile(r'(?<=opt_params:\s\[)(.*)(?=\]\nopt_value)',re.DOTALL)
  opt_val_pattern = re.compile(r'(?<=opt_value:\s)(.*)(?=optimizer_time)',re.DOTALL)
  with open(path_to_log) as log_file:
    print(f"loading {path_to_log}...")
    log = log_file.read()
    log_file.close
  for match in opt_params_pattern.finditer(log):
      opt_params = list(map(float, match.group().split()))
  for match_val in opt_val_pattern.finditer(log):
      opt_val = float(match_val.group())
  parse_opt_expectation_values(path,opt_params,opt_val,maxp,pstep)
  
def parse_opt_expectation_values(path, opt_params, opt_val, maxp = 7, pstep=1):
  energy_vals = []
  best_energy_vals = []
  best_params = []
  p = 0
  for file_name in sorted([file for file in os.listdir(path) if (file.endswith('.json') and file!="eigenvalues.json" and file!="opt_data.json"  and file!="offset.json")], key=lambda x: int(x.split("scan_p")[1].split("_")[0])):
    if pstep>1:
      print(file_name)
      print(p)
    p += pstep
    print(f"loading {file_name} and extracting data....")
    with open(path + "/" + file_name) as json_file:
      data = json.load(json_file)
      #best found energies
      best = data["best"]["params"][0]
      expectation = data["best"]["expectation"][0]
      best_params.append(best)
      best_energy_vals.append(expectation)

      #optimised param energies
      if p ==maxp:
        energy_vals.append(opt_val)
      else:
        beta = opt_params[:p]
        gamma = opt_params[maxp:(maxp+p)]
        current_p_params = beta + gamma
        i,j = find_nearest_in_scan_grid(data["scan_grid"],current_p_params)
        fixed_energy = data["expecation_values"][i][j]
        energy_vals.append(fixed_energy)

  #save data as json
  print("saving results as json...")
  opt_data = {"params": opt_params,"expectation_values":energy_vals,"best":best_params,"best_expectation_values":best_energy_vals}
  with open(f"{path}/opt_data.json", 'w') as f:
      json.dump(opt_data, f)

def find_nearest_in_scan_grid(scan_grid,params):
    best = 1000
    idx = [0,0]
    for i in range(len(scan_grid)):
      for j in range(len(scan_grid[0])):
        tuples = zip(params,scan_grid[i][j])
        sum = 0
        for tuple in tuples:
          sum += abs(tuple[0]-tuple[1])
        if sum<best:
          best = sum
          idx = [i,j]
    return idx[0], idx[1]


