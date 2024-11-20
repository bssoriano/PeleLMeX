import numpy as np
import fnmatch
import os
import matplotlib
import matplotlib.pyplot as plt
import re
import h5py
import pandas as pd
import cantera as ct
import yt
from matplotlib.pyplot import cm
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.size'] = 16

def Y2X(comp,gas_filename):
    gas = ct.Solution(gas_filename)
    gas.Y = comp
    return gas.X[gas.species_index('H2O')]

def compute_RH(X_h2o,T,P):
    #Inputs: T [K]; P [kPa]

    #Tabulated H2O saturation pressure (kPa) as a function of temperature (Celsius)
    #Extracted from: https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
    h2o_sat_p = np.array([0.607,0.61165,0.70599,0.81355,1.2282,1.5990,2.0647,2.3393,3.1699,4.2470,5.3251,7.3849,9.1124,12.352,15.022,19.946,31.201,47.414,70.182,87.771,101.42])
    temp_c = np.array([0.0,0.01,2,4,10,14,18,20,25,30,34,40,44,50,54,60,70,80,90,96,100])
    temp_K = temp_c+273.15

    sat_p     = np.interp(T,temp_K,h2o_sat_p)
    partial_p = X_h2o*P

    return partial_p/sat_p

gas_filename = '/home/bsouzas/my_LMeX/PelePhysics/Mechanisms/humid_air/mechanism.yaml'
gas = ct.Solution(gas_filename)
species_names = gas.species_names

istart = 0
iskip  = 4

cwd = os.getcwd()


istart = 0
iskip  = 1
folders=[
        'humidification/',
        ]

fig, ax = plt.subplots()
# ax2 = ax.twinx(nrows=2)

for path in folders:


    #digits=re.findall(r'\d+', path)
    #label = "{}.{}".format(digits[0],digits[1])
    label = ''

    list_output=fnmatch.filter(os.listdir(path), '*.h5')
    list_output.sort()    
    print("Found output   : {}".format(list_output))

    #list_process=list_output[istart:len(list_output):iskip]
    # list_process=[list_output[-14],list_output[-13]]
    # list_process=[list_output[-1]]
    #list_process=list_output
    list_process=['plt234678']
    print("To be processed: {} \n".format(list_process))
    
    n_files = len(list_process)

    for ifile,file in enumerate(list_process):
        print(file)

        digits=re.findall(r'\d+', file)
        ds = yt.load("{}/plt{}".format(path,digits[0]))
        time = float(ds.current_time)

        data = ds.all_data()
        ncells = len(data['z'])

        RH = np.zeros(ncells)
        count = 0
        #ncells = 100
        for i in range(ncells):
            composition = ",".join(['{}:{}'.format(name,float(data["Y({})".format(name)][i])) for name in species_names])
            X_h2o = Y2X(composition,gas_filename)
            RH[i] = compute_RH(X_h2o,data['temp'][i],data['RhoRT'][i]*1.e-3)


        idx = data['volFrac'] > 0.9
        cmap = 'RdBu_r'
        cm = plt.cm.get_cmap(cmap)
        print(RH[idx])

        cb = ax.scatter(data['z'][idx]*1.e+3,RH[idx],c=data['temp'][idx],vmin=273, vmax=300,cmap=cm)
        cbar = fig.colorbar(cb)
        cbar.ax.set_ylabel('T [K]')
        
#ax[1].ticklabel_format(style='sci', axis='y', scilimits=(0,0))

ax.set_xlabel('z [mm]')
ax.set_ylabel('RH [%]')
ax.set_ylim(0.4,1.0)
ax.set_xlim(10.0,745.0)
ax.grid(axis='y')
fig.tight_layout()
# ax.legend()
# plt.subplots_adjust(left=0.18, bottom=0.16, right=0.85, top=0.92)
plt.savefig('humidification/plots/RH_scatter.png',dpi=300)
# plt.show()
    
