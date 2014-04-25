#!/usr/bin/python

import os
import subprocess
from __builtin__ import range
import threading
from math import pi


k_b = 1.3806488e-23
diameter = 5e-9
T = 300.0
viscosity = 8.9e-4
D = k_b*T/(3.0*pi*viscosity*diameter);

time = 1e-6
nout = 500
#run_dir = "/scratch/robinsonm/git/paper_crowding"
run_dir = "/home/mrobins/git/paper_crowding"
#data_dir = "/mi/share/scratch/robinsonm/data/crowding/self_crowding"
data_dir = "/home/mrobins/tmp"

def run(k_s, sl_div_diam, vol_ratio):
    new_dir = str(k_s)+"_"+str(sl_div_diam)+"_"+str(vol_ratio)
    print "creating subdir ",new_dir
    if not os.path.exists(data_dir + "/" + new_dir):
        os.makedirs(data_dir + "/" + new_dir)
        subprocess.call([run_dir+"/self_crowding",data_dir + "/" + new_dir,
                     str(time),str(nout),str(k_s),str(sl_div_diam),str(vol_ratio)])

def run_sweep():
    k_s_sweep = [10**(x-5) for x in range(10)]
    print "k_s = ",k_s_sweep
    sl_div_diam_sweep = [5.0*(x+1)/1000 for x in range(20)]
    print "sl_div_diam = ",sl_div_diam_sweep
    vol_ratio_sweep = [5.0*(x+1)/100 for x in range(10)]
    print "vol_ratio = ",vol_ratio_sweep

    params = []
    for k_s in k_s_sweep:
        for sl_div_diam in sl_div_diam_sweep:
            for vol_ratio in vol_ratio_sweep:
                params.append((k_s,sl_div_diam,vol_ratio))
                
    n_cpu = 4
    num_its = int(len(params)/n_cpu)
    for i in range(num_its):
        threads = [threading.Thread(target = run, args = params[i*n_cpu+j]) for j in range(n_cpu)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

    


if __name__ == '__main__':
    run_sweep()
