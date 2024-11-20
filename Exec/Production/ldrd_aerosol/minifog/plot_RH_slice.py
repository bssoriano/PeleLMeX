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
import glob
import re
from joblib import Parallel, delayed
from yt.visualization.base_plot_types import get_multi_plot
matplotlib.rcParams['font.family'] = 'serif'

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

def compute_Xsat_from_T(T,P,rh):
    #Tabulated H2O saturation pressure (kPa) as a function of temperature (Celsius)
    #Extracted from: https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
    h2o_sat_p = np.array([0.61165,0.70599,0.81355,1.2282,1.5990,2.0647,2.3393,3.1699,4.2470,5.3251,7.3849,9.1124,12.352,15.022,19.946,31.201,47.414,70.182,87.771,101.42])
    temp_c = np.array([0.01,2,4,10,14,18,20,25,30,34,40,44,50,54,60,70,80,90,96,100])
    temp_K = temp_c+273.15

    sat_p     = np.interp(T,temp_K,h2o_sat_p)
    print('saturation pressure:',sat_p)

    return rh*sat_p/P

def process_RH(i,j,Temp,Pressure,Ys,gas_filename):
    gas = ct.Solution(gas_filename)
    Y = np.zeros(gas.n_species)
    for ispec,species in enumerate(gas.species_names):
        Y[ispec] = Ys[j,i,ispec]
    T = Temp[j,i]
    P = Pressure[j,i]
    
    if(np.isnan(T) or np.isnan(P) or T < 100.0):
        X_H2O = 0.0
    else:
        gas.TPY = T, P, Y
        X_H2O = gas.X[gas.species_index('H2O')]

    return compute_RH(X_H2O,T,P)


P = 101.325
T = 35.+273.15
rh = 1.
print('H2O molefraction for inlets:', compute_Xsat_from_T(T,P,rh))

T = 301.5
rh = 0.89
X_h2o = compute_Xsat_from_T(T,P,rh)
print('H2O molefraction for initial condition:', X_h2o)


T = 299.
# X_h2o = 0.04
print('Relative humidity for T={} and X-H2O={}: '.format(T,X_h2o),compute_RH(X_h2o,T,P))

gas_filename = '/home/bsouzas/my_LMeX/PelePhysics/Mechanisms/humid_air/mechanism.yaml'
gas = ct.Solution(gas_filename)

plot_size = 1200.e-3
resolution = 10000.e-6 # [m]
resolution = 5000.e-6 # [m]

nx = int(plot_size/resolution)
ny = nx
print(nx)

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()
path = 'humidification/'
plt_file = 'plt'

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()    
# print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
list_process=[list_process[-1]]
# list_process=list_process[-4:]
# list_process=['']
# list_process=['plt234678']



# find all plt files to be processed
time_save=[]
time_save_path=glob.glob(path + '/' + plt_file + '*') # path+filename
for ii in time_save_path:
  ii=os.path.basename(ii).split(plt_file)[-1]               # extract numbers from filename
  time_save.append(os.path.basename(ii))                  # unsorted
time_save.sort(key=lambda f: int(re.sub('\D', '', f)))    # sorted
filelen_time = len(time_save)
filelen_time = 0
# print("To be processed: {} \n".format(time_save))

scalars = ['Y(H2O)','mag_vel','temp']
# scalars = ['mag_vel']

os.makedirs(path+'/plots/',exist_ok=True)
os.makedirs(path+'/plots/animation/',exist_ok=True)
minvals = []; maxvals = []
for n in range(filelen_time):
        # ds = yt.load("{}/{}".format(path,name))
    

    ds1 = yt.load("{}/{}{}".format(path,plt_file,time_save[n]))  # load data
    time = float(ds1.current_time)
    ds1.print_stats()
    scalars = [name[1] for name in ds1.field_list if('Y(' in name[1] or 'temp' in name[1] or 'RhoRT' in name[1])]
    scalars.append('volFrac')
    print(scalars)

    slc1 = yt.SlicePlot(
            ds1,
            #"z",
            #center=[1.0e-4,0.0,5.e-3],
            "y",
            center=[1.0e-4,1.0e-4,750.e-3/2.],
            #center=[1.0e-4,-174.2e-3,750.e-3/2.],
            #fields=[scalar,'volFrac','density'],
            fields=scalars
        )

    slc_frb1 = slc1.data_source.to_frb((plot_size), nx)

    Temp     = np.array(slc_frb1['temp'])
    Pressure = np.array(slc_frb1['RhoRT'])*1.e-3 # Pa to kPa
    Ys = np.zeros(shape=(nx,nx,3))
    for ispec,species in enumerate(gas.species_names):
        Ys[:,:,ispec] = np.array(slc_frb1["Y({})".format(species)])

    # Compute relative humidity (RH)
    RH = np.zeros(shape=(nx,nx))
    NCPUs = 32
    tmp = Parallel(n_jobs=NCPUs,verbose=100)(delayed(process_RH)(i,j,Temp,Pressure,Ys,gas_filename) for i in range(nx) for j in range(nx))
    print(np.min(tmp),np.max(tmp))

    RH = np.reshape(tmp,(nx,nx))

#    for i in range(nx):
#        for j in range(nx):
#            Y = np.zeros(gas.n_species)
#            for ispec,species in enumerate(gas.species_names):
#                Y[ispec] = np.array(slc_frb1["Y({})".format(species)])[j,i]
#            T = Temp[j,i]
#            P = Pressure[j,i] 
#            if(np.isnan(T)==False and np.isnan(P)==False and np.array(slc_frb1["volFrac"])[j,i] > 0.9):
#                gas.TPY = T, P, Y
#                X_H2O = gas.X[gas.species_index('H2O')]
#        
#                RH[j,i] = compute_RH(X_H2O,T,P)

    #I = np.where(slc_frb1["volFrac"] < 0.9)
    I = np.where(RH < 0.1)
    RH[I] = np.nan 
     
    extent = -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, 0.0, (plot_size)*1.e+3
    xgrid = np.linspace( -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, nx)
    ygrid = np.linspace( 0.0                 , plot_size*1.e+3    , ny)

    fig, ax1 = plt.subplots()
    
    minval = np.nanmin(RH)
    maxval = np.nanmax(RH)
    minval = 0.4
    maxval = 1.0
    cmap = 'turbo'
    cmap = 'RdBu_r'
    CS = ax1.imshow(RH[:,5:-5],vmin=minval,vmax=maxval,cmap=cmap,extent=extent,origin='lower')
    cbar = fig.colorbar(CS, ax=ax1, format='%03.1f', ticks=np.linspace(minval,maxval,5),shrink=0.6,extend='both')

    cbar.ax.set_ylabel('RH [%]')

    ax1.set_xlabel('x [mm]')
    ax1.set_ylabel('z [mm]')
    ax1.set_xlim(-plot_size*1.e+3/2.0,plot_size*1.e+3/2.0)
    ax1.set_ylim(0.,plot_size*1.e+3)
    ax1.axis("off")
    ax1.annotate('t = {:1.1f} [s]'.format(time),xy=(-150.0,1000.),fontsize=14)
    plt.tight_layout()

    fig.savefig(path+'/plots/animation/RH_{:05d}.png'.format(n),dpi=300)
    #fig.savefig('test_RH.png')
# plt.show()
