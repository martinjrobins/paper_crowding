#!/usr/bin/python

import os
import subprocess
from __builtin__ import range
import threading
from math import pi
import numpy
import csv
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mayavi.mlab as mlab

k_b = 1.3806488e-23
diameter = 5e-9
T = 300.0
viscosity = 8.9e-4
D = k_b*T/(3.0*pi*viscosity*diameter);

time = 1e-6
nout = 500
run_dir = "/scratch/robinsonm/git/paper_crowding"
#run_dir = "/home/mrobins/git/paper_crowding"
data_dir = "/mi/share/scratch/robinsonm/data/crowding/self_crowding"
#data_dir = "/home/mrobins/tmp"

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
            
def calc_average_rdf(k_s, sl_div_diam, vol_ratio):
    subdir = str(k_s)+"_"+str(sl_div_diam)+"_"+str(vol_ratio)
    dir = data_dir + "/" + subdir
    print "reducing dir ",dir
    data = numpy.zeros((100,2,400))
    for i in range(100,500):
        data[:,:,i-100] = numpy.loadtxt(dir+ "/rdf%05d.csv"%i, delimiter=',', usecols = (0,1))        
    return numpy.average(data, 2)

def reduce(k_s, sl_div_diam, vol_ratio):
    subdir = str(k_s)+"_"+str(sl_div_diam)+"_"+str(vol_ratio)
    dir = data_dir + "/" + subdir
    print "reducing dir ",dir
    data = numpy.loadtxt(dir+ "/" + 'msd.csv', delimiter=',')
    valid_range = data[20:,:]
    t = valid_range[:,1]
    msd = valid_range[:,2]
    #flux = valid_range[:,3]
    A = numpy.array([t,numpy.ones(len(t))])
    msd_fit = numpy.linalg.lstsq(A.T, msd)[0]
    #flux_fit = numpy.linalg.lstsq(A.T, flux)[0]
    rdf = numpy.zeros((100,2,400))
    for i in range(100,500):
        rdf[:,:,i-100] = numpy.loadtxt(dir+ "/rdf%05d.csv"%i, delimiter=',', usecols = (0,1))        
    numpy.savetxt(dir+ "/rdf_average.csv", numpy.average(rdf,2), delimiter=',')
    return msd_fit[0]/6.0

    
    
def reduce_sweep():
    k_s_sweep = [10**(x-5) for x in range(10)]
    print "k_s = ",k_s_sweep
    sl_div_diam_sweep = [5.0*(x+1)/1000 for x in range(20)]
    print "sl_div_diam = ",sl_div_diam_sweep
    vol_ratio_sweep = [5.0*(x+1)/100 for x in range(10)]
    print "vol_ratio = ",vol_ratio_sweep

    params = []
    Dmsd = numpy.zeros((10,20,10))
    with open(data_dir + "/D_row.csv", 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='#')
        for k_s,i in zip(k_s_sweep,range(10)):
            for sl_div_diam,j in zip(sl_div_diam_sweep,range(20)):
                for vol_ratio,k in zip(vol_ratio_sweep,range(10)):
                    Dmsd[i,j,k] = reduce(k_s,sl_div_diam,vol_ratio)
                    writer.writerow((k_s,sl_div_diam,vol_ratio,Dmsd[i,j,k]))
    
    numpy.save(data_dir + "/D_numpy.npy", Dmsd)

def showDsweep():
    k_s_sweep = [10**(x-5) for x in range(10)]
    sl_div_diam_sweep = [5.0*(x+1)/1000 for x in range(20)]
    vol_ratio_sweep = [5.0*(x+1)/100 for x in range(10)]    
    
    Dmsd = numpy.load(data_dir + "/D_numpy.npy")
    kDa = 1.660538921e-30;
    mass = 40.0*kDa; 
    viscosity = 8.9e-4; 
    diameter = 5e-9; 
    T = 300.0; 
    Dbase = k_b*T/(3.0*numpy.pi*viscosity*diameter);
    Dmsd = Dmsd/Dbase
    
    mlab.figure(1, size=(800, 800), fgcolor=(1, 1, 1),
                                    bgcolor=(0.5, 0.5, 0.5))
    mlab.clf()
    contours = numpy.arange(0.01,2,0.2).tolist()
    obj = mlab.contour3d(Dmsd,contours=contours,transparent=True,vmin=contours[0],vmax=contours[-1])
    outline = mlab.outline(color=(.7, .7, .7),extent=(0,10,0,20,0,10))
    axes = mlab.axes(outline, color=(.7, .7, .7),
            nb_labels = 5,
            ranges=(k_s_sweep[0], k_s_sweep[-1], sl_div_diam_sweep[0], sl_div_diam_sweep[-1], vol_ratio_sweep[0], vol_ratio_sweep[-1]), 
            xlabel='spring stiffness', 
            ylabel='step length',
            zlabel='volume ratio')
    mlab.colorbar(obj,title='D',nb_labels=5)

    mlab.show()
    

        
def plots():
    k_s_sweep = [10**(x-5) for x in range(10)]
    sl_div_diam_sweep = [5.0*(x+1)/1000 for x in range(20)]
    sl_div_diam_sweep.reverse()
    vol_ratio = 0.3
    plt.figure(figsize=(6,4.5))

    for k_s,i in zip(k_s_sweep,range(10)):
        plt.clf()
        for sl_div_diam,j in zip(sl_div_diam_sweep,range(20)):
            subdir = str(k_s)+"_"+str(sl_div_diam)+"_"+str(vol_ratio)
            dir = data_dir + "/" + subdir
            rdf = numpy.loadtxt(dir+ "/rdf_average.csv", delimiter=',', usecols = (0,1))
            plt.plot(rdf[:,0],rdf[:,1],label="step = %f"%(sl_div_diam),color=(1-j/20.0, 1-j/20.0, 1.0))
        plt.xlabel('$r$')
        plt.ylabel('RDF')
        print "saving file /rdf_%f.pdf"%(k_s)
        plt.savefig(data_dir + "/rdf_%f.pdf"%(k_s))
    

if __name__ == '__main__':
    #reduce_sweep()
    plots()
