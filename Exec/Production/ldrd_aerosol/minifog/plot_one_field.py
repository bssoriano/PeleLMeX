import numpy as np
from matplotlib.colors import LogNorm
import os 
import fnmatch
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker,cm

import yt
from yt.visualization.base_plot_types import get_multi_plot


def compute_RH(X_h2o,T,P):
    #Tabulated H2O saturation pressure (kPa) as a function of temperature (Celsius)
    #Extracted from: https://www.engineeringtoolbox.com/water-vapor-saturation-pressure-d_599.html
    h2o_sat_p = np.array([0.61165,0.70599,0.81355,1.2282,1.5990,2.0647,2.3393,3.1699,4.2470,5.3251,7.3849,9.1124,12.352,15.022,19.946,31.201,47.414,70.182,87.771,101.42])
    temp_c = np.array([0.01,2,4,10,14,18,20,25,30,34,40,44,50,54,60,70,80,90,96,100])
    temp_K = temp_c+273.15

    sat_p     = np.interp(T,temp_c,h2o_sat_p)
    partial_p = X_h2o*P

    return partial_p/sat_p

orient = "vertical"

plot_size = 1200.e-3
resolution = 1000.e-6 # [m]

nx = int(plot_size/resolution)
ny = nx

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()
path = '/archive/bsouzas/ldrd_aerosol/minifog/humidification_after_inj/'
# path = '/home/bsouzas/ldrd_aerosol/humidification/'

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
list_process=[list_process[-1]]
# list_process=list_process[-10:]
# list_process=['plt01105']
# list_process=['plt234678']
print("To be processed: {} \n".format(list_process))

scalars = ['Y(H2O)','mag_vel','temp']
scalars = ['Y(H2O)']
scalars = ['temp']

os.makedirs(path+'/plots/',exist_ok=True)
minvals = []; maxvals = []
for i,name in enumerate(list_process):
        # ds = yt.load("{}/{}".format(path,name))


    ds1 = yt.load("{}/{}".format(path,name))  # load data
    time = float(ds1.current_time)
    if(time > 599. and time < 601.):
        print(name,time)
    ds1.print_stats()
    # scalars = [name[1] for name in ds1.field_list]
    
    #ad = ds1.all_data()
    #scalars = []
    #for scalar in scalars:
    #    if('temp' in scalar or 'mag_vel' in scalar):
    #        print(scalar)
    #        print("Min {}: {:1.2e}".format(scalar,np.min(ad[scalar])))
    #        print("Max {}: {:1.2e}".format(scalar,np.max(ad[scalar])))

    for iscal,scalar in enumerate(scalars):

        os.makedirs(path+'/plots/'+scalar,exist_ok=True)

        slc1 = yt.SlicePlot(
            ds1,
            #"z",
            #center=[1.0e-4,0.0,5.e-3],
            "y",
            center=[1.0e-4,1.0e-4,750.e-3/2.],
            # center=[1.0e-4,-174.2e-3,750.e-3/2.],
            fields=[scalar,'volFrac','density'],
        )


        slc_frb1 = slc1.data_source.to_frb((plot_size), nx)


        slc_temp1   = np.array(slc_frb1[scalar])
        slc_vfrac   = np.array(slc_frb1["volFrac"])
        slc_rho     = np.array(slc_frb1["density"])
        

        if('rhoMixFrac' in scalar):
            I = np.where(slc_rho < 1.e-6)
            slc_rho[I] = 1.
            slc_temp1 = slc_temp1/slc_rho
            print(scalar)
            print(slc_temp1)
        
        I = np.where(slc_frb1["volFrac"] < 0.9)
        slc_temp1[I] = np.nan 
     
        extent = -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, 0.0, (plot_size)*1.e+3
        xgrid = np.linspace( -plot_size*1.e+3/2.0, plot_size*1.e+3/2.0, nx)
        ygrid = np.linspace( 0.0                 , plot_size*1.e+3    , ny)

        fig, ax1 = plt.subplots()
        if(scalar == 'temp'):
            minval = 300.

        if(i == 0 and scalar != 'HeatRelease'):
            cmap = 'turbo'
            cmap = 'RdBu_r'
            if(scalar == 'temp'):
                minvals.append(273.)
            else:
                minvals.append(np.nanmin(slc_temp1))
            maxvals.append(np.nanmax(slc_temp1))
            #maxvals.append(470.)
        elif(i == 0 and scalar == 'HeatRelease'):
            minvals.append(1.e+4)
            maxvals.append(1.e+10)
            idx = slc_temp1 < 1.e+4
            slc_temp1[idx] = 1.e+4
            # minvals.append(np.nanmin(slc_temp1))
            # maxvals.append(np.nanmax(slc_temp1))
            cmap = 'dusk'

        print(np.shape(slc_temp1))
        # cmap='terrain'
        # CS = ax1.contourf(xgrid,ygrid,slc_temp1,levels=np.linspace(0.,8.e-3,20),cmap=cm.turbo)
        if(scalar == 'HeatRelease'):
            CS = ax1.imshow(slc_temp1,norm=LogNorm(vmin=minvals[iscal], vmax=maxvals[iscal]),cmap=cmap,extent=extent,origin='lower')
            cbar = fig.colorbar(CS, ax=ax1, format='%03.1e', ticks=np.geomspace(minvals[iscal],maxvals[iscal],5),shrink=0.8)
        else:
            CS = ax1.imshow(slc_temp1.T,vmin=minvals[iscal],vmax=maxvals[iscal],cmap=cmap,extent=extent,origin='lower')
            cbar = fig.colorbar(CS, ax=ax1, format='%03.1e', ticks=np.linspace(minvals[iscal],maxvals[iscal],5),shrink=0.5)

        cbar.ax.set_ylabel(scalar)

        ax1.set_xlabel('x [mm]')
        ax1.set_ylabel('z [mm]')
        ax1.set_xlim(-plot_size*1.e+3/2.0,plot_size*1.e+3/2.0)
        ax1.set_ylim(0.,plot_size*1.e+3)
        ax1.axis("off")
        ax1.annotate('t = {:1.1f} [s]'.format(time),xy=(-150.0,1000.),fontsize=14)
        plt.tight_layout()

        #fig.savefig(path+'/plots/'+scalar+'/Inst_{}_t_{:1.1f}ms.png'.format(scalar,time*1.0e+3),dpi=300)
plt.show()
