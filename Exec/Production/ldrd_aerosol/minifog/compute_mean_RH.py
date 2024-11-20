import numpy as np
from matplotlib.colors import LogNorm
import os 
import fnmatch
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker,cm
import cantera as ct
import yt
from yt.visualization.base_plot_types import get_multi_plot

def Y2X(comp):
    gas = ct.Solution('mechanism.yaml')
    gas.Y = comp
    return gas.X


gas = ct.Solution('mechanism.yaml')
species_names = gas.species_names


orient = "vertical"

plot_size = 1200.e-3
resolution = 1000.e-6 # [m]

nx = int(plot_size/resolution)
ny = nx

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()
path = 'humidification/'

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
list_process=[list_process[-1]]
# list_process=list_process[-4:]
# list_process=['plt01105']
list_process=['plt161197']
print("To be processed: {} \n".format(list_process))

scalars = ['Y(H2O)','mag_vel','temp','density','RhoRT']

Pvs = 3.169 # vapor saturation pressure [kPa]
for i,name in enumerate(list_process):
        # ds = yt.load("{}/{}".format(path,name))


    ds1 = yt.load("{}/{}".format(path,name))  # load data
    time = float(ds1.current_time)
    ds1.print_stats()
    # scalars = [name[1] for name in ds1.field_list]

    ad = ds1.all_data()
    for scalar in scalars:
        if('temp' in scalar or 'mag_vel' in scalar):
            print(scalar)
            print("Min {}: {:1.2e}".format(scalar,np.min(ad[scalar])))
            print("Max {}: {:1.2e}".format(scalar,np.max(ad[scalar])))

    Xs = {}
    for name in species_names:
        Xs[name] = []

    idx_is_fluid = ad['temp'][ad['volFrac'] > 0.1]

    for i in range(len(ad['temp'])):
        if(ad['volFrac'][i] > 0.1):
            composition = ",".join(['{}:{}'.format(name,float(ad["Y({})".format(name)][i])) for name in species_names])
            local_comp_X = Y2X(composition)
            for j,name in enumerate(species_names):
                Xs[name].append(local_comp_X[j])
    
    Pv = Xs['H2O']*ad['RhoRT']
    RH = Pv/Pvs

    print("Mean density: ",np.mean(ad['density']))
    print("Mean H2O massfraction: ",np.mean(ad['Y(H2O)']))
    print("Target H2O massfraction: ",np.max(ad['Y(H2O)']))
    print("Mean relative humidity: ", np.mean(RH))



