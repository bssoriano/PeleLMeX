import sys
sys.path.append('/home/bsouzas/Python_Scripts/Python_modules/')
sys.path.append('/home/bsouzas/Python_Scripts/mod_otros/pyevtk/')
from pyevtk.hl import gridToVTK

import yt
import fnmatch
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import h5py
import h5py as h5
import hdf5plugin
import time
from mpi4py import MPI

def load_data(cube,scalar_name,idxl,idxr,myid):

    chunk_size = int((idxr[2]-idxl[2])/NCPUs)
    idxZlo = chunk_size*myid+idxl[2]
    idxZhi = idxZlo+chunk_size

    idxZlo, idxZhi = get_mpi_load(idxr[2]-idxl[2],myid)

    # print(cube[scalar_name][idxl[0]:idxr[0],idxl[1]:idxr[1],0])

    fld = np.array(cube[scalar_name][idxl[0]:idxr[0],idxl[1]:idxr[1],idxZlo:idxZhi])
    return fld

def get_mpi_load(data_size,myid):
    chunk_size = int(data_size/NCPUs) 
    idxZlo = chunk_size*myid
    idxZhi = idxZlo + (chunk_size)
    return idxZlo, idxZhi


comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()
NCPUs = comm.Get_size()

if(rank == 0):print("Running on {} processors".format(NCPUs))

precision = 'f4' # Precision of h5 arrays

field_list = ['volFrac','temp','Y(H2O)','mag_vel']

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()
path = 'humidification/'

list_output=fnmatch.filter(os.listdir(path), 'plt*')
list_output.sort()
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
list_process=[list_process[-2],list_process[-1]]
list_process=['plt184011']

write_vtk = True
write_h5  = False

t_start = time.time()

# %% ############################################################################
# Find times in files
for filename in list_process:
    # filename = os.getcwd()+"/plt"+("%05d"%n)+"_scaldis"

    ds = yt.load(path+filename)
    if(rank == 0):
        print(filename)
        ds.print_stats()
        print()

    # field_list = [v[1] for v in ds.field_list]

    level = 0 #ds.max_level
    dims = ds.domain_dimensions * ds.refine_by ** level

    if (ds.domain_dimensions[2] == 1):
        dims[2] = 1

    # t1 = time.time()
    # cube = ds.covering_grid(
    #     level,
    #     left_edge=ds.domain_left_edge,
    #     dims=dims,
    #     )

    # if(rank == 0):print('Time in covering_grid: {:1.3f} [s]'.format(time.time() - t1))
    
    t1 = time.time()
    dx = (ds.domain_right_edge[0]- ds.domain_left_edge[0])/dims[0]
    dy = (ds.domain_right_edge[1]- ds.domain_left_edge[1])/dims[1]
    dz = (ds.domain_right_edge[2]- ds.domain_left_edge[2])/dims[2]   
    xl = np.arange(ds.domain_left_edge[0], ds.domain_right_edge[0], dx , dtype='float64')
    yl = np.arange(ds.domain_left_edge[1], ds.domain_right_edge[1], dy , dtype='float64')
    zl = np.arange(ds.domain_left_edge[2], ds.domain_right_edge[2], dz , dtype='float64')
            
    idxl = 0 #np.argmin(abs(xl-domain_left_edge[0]))
    idyl = 0 #np.argmin(abs(yl-domain_left_edge[1]))
    idzl = 0 #np.argmin(abs(zl-domain_left_edge[2]))

    idxr = len(xl) #np.argmin(abs(xl-domain_right_edge[0]))
    idyr = len(yl) #np.argmin(abs(yl-domain_right_edge[1]))
    idzr = len(zl) #np.argmin(abs(zl-domain_right_edge[2]))
    
    #if (ds.domain_dimensions[2] == 1):
    #    idzr = idzr+1
        
    x = xl
    y = yl #[idyl:idyr]
    z = zl #[idzl:idzr]
    if(rank == 0):print('Time to initialize grids: {:1.3f} [s]'.format(time.time() - t1))

    nx = len(x)
    ny = len(y)
    nz = len(z)

    t1 = time.time()
    domlo = ds.domain_left_edge
    
    idxZlo_global, idxZhi_global = get_mpi_load(nz,rank)

    decomposed_domain = [domlo[0],domlo[1],z[idxZlo_global] * ds.units.code_length]
    ncells            = [      nx,      ny,            int(nz/NCPUs)]

    cube = ds.covering_grid(
        level,
        left_edge=decomposed_domain,
        dims=ncells,
        )

    if(rank == 0):print('Time in covering_grid: {:1.3f} [s]'.format(time.time() - t1))

    recvbuf = None
    if(rank == 0):
        recvbuf = np.empty([NCPUs,nx,ny,int(nz/NCPUs)])

    if (write_vtk):
        pointdata = {}
        for i in range(len(field_list)):
            # fld = np.array(cube[field_list[i]][idxl:idxr,idyl:idyr,idzl:idzr])
            # pointdata.update({field_list[i] : fld})

            t2 = time.time()
            fld = load_data(cube,field_list[i],[idxl,idyl,idzl],[idxr,idyr,idzr],0)
            if(rank == 0): print('Time load data: {:1.3f} [s]'.format(time.time()-t2))

            t2 = time.time()

            comm.Gather(fld, recvbuf, root=0)
            if(rank == 0): print('Time in communication: {:1.3f} [s]'.format(time.time()-t2))

            if(rank == 0):
                t2 = time.time()
                idxZlo, idxZhi = get_mpi_load(nz,rank)

                data_output = np.zeros((nx,ny,nz))
                data_output[:,:,idxZlo:idxZhi] = fld[:,:,idxZlo:idxZhi]

                for procid in range(1,NCPUs):
                    idxZlo_global, idxZhi_global = get_mpi_load(nz,procid)
                    data_output[:,:,idxZlo_global:idxZhi_global] = recvbuf[procid,:,:,idxZlo:idxZhi]

                print('Time to gather all data: {:1.3f} [s]'.format(time.time()-t2))

                # plt.imshow(cube[field_list[i]][int(dims[0]/2),:,:],origin="lower")
                # plt.imshow(data_output[int(dims[0]/2),:,:],origin="lower")
                # plt.show()

                pointdata[field_list[i]] = data_output


        if(rank == 0):
            t2 = time.time()
            gridToVTK(path+filename+"_max_level_"+str(level), x, y, z, pointData = pointdata)
            print('Time to write VTK: {:1.3f} [s]'.format(time.time()-t2))
    
    elif(write_h5):
        t1 = time.time()
        if(rank == 0):
            print("Data will be written in HDF5 format")
            # Create h5 file
            H5FILE = h5.File(filename+".h5",'w')
            t = H5FILE.create_dataset('time', data=float(ds.current_time), dtype='f8')
            H5_x = H5FILE.create_dataset('x', data=x, dtype='f8')
            H5_y = H5FILE.create_dataset('y', data=y, dtype='f8')
            H5_z = H5FILE.create_dataset('z', data=z, dtype='f8')

            print("Finished writing grid data. Moving on the scalar fields...")
        for i in range(len(field_list)):
            # print(field_list[i])
            # pointdata = {}

            t2 = time.time()
            # fld = np.array(cube[field_list[i]][idxl:idxr,idyl:idyr,idzl:idzr])

            fld = load_data(cube,field_list[i],[idxl,idyl,idzl],[idxr,idyr,idzr],0)

            comm.Gather(fld, recvbuf, root=0)
            if(rank == 0):
                data_output = np.zeros((nx,ny,nz))

                idxZlo, idxZhi = get_mpi_load(nz,rank)
                data_output[:,:,idxZlo:idxZhi] = recvbuf[rank,:,:,idxZlo:idxZhi]

                for procid in range(1,NCPUs):
                    idxZlo_global, idxZhi_global = get_mpi_load(nz,procid)
                    data_output[:,:,idxZlo_global:idxZhi_global] = recvbuf[procid,:,:,idxZlo:idxZhi]

                print('Time to load 1 scalar: {:1.3f} [s]'.format(time.time()-t2))
   
                t2 = time.time()
                H5_aux = H5FILE.create_dataset(field_list[i], (nx, ny, nz), dtype=precision, **hdf5plugin.Zfp(reversible=True))
                H5_aux[:,:,:] = data_output[:,:,:]
                print('Time write 1 scalar: {:1.3f} [s]'.format(time.time()-t2))

        if(rank == 0):print('Total time to write {} scalars: {:1.3f} [s]'.format(len(field_list),time.time()-t1))

        
if(rank == 0):print('Total time: {:1.3f} [s]'.format(time.time()-t_start))

