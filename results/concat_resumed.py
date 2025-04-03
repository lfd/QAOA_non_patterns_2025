import sys, os
import shutil

def parse_first_half_and_make_outdirs(path,outpath):
    for dir_name in [dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
        print(dir_name)
        outdir = f"{outpath}/{dir_name}"
        os.mkdir(outdir)
        print(f"making dir at {outdir}")
        for file in os.listdir(path + "/" + dir_name):
            shutil.copy(path + "/" + dir_name + "/" + file,outdir)
            print(f"copying {file} to {outdir}")

def parse_second_half_and_add_to_outdirs(path,outpath):
    for dir_name in [dir for dir in os.listdir(path) if os.path.isdir(path + "/" + dir)]:
        for outdir_name in [outdir for outdir in os.listdir(outpath) if os.path.isdir(outpath + "/" + outdir)]:
            if outdir_name[0:30]==dir_name[0:30]:
                matchdir = outdir_name
        for file in os.listdir(path + "/" + dir_name):
            shutil.copy(path + "/" + dir_name + "/" + file,outpath + "/" + matchdir)
            print(f"copying {file} to {outpath + '/' + matchdir}")

def concat_dirs(dir1,dir2,out):
    parse_first_half_and_make_outdirs(dir1,out)
    parse_second_half_and_add_to_outdirs(dir2,out)

if __name__=="__main__":

    print("script for concating resumed landscape scans") 
    #concat_dirs("./paper_experiments/maxcut/sequential","./maxcut_sequential_p8_to_21","./paper_experiments/higher_depth_maxcut/sequential")