import yt
import cantera as ct
import numpy as np
import fnmatch
import os


def get_line(filename,scalars):
    print("Writing the following scalars to file: ",scalars)
    ray = ds.ortho_ray(2,(0.0,0.0e-3))
    f = ray.save_as_dataset(filename,fields=scalars)
    return ray

cwd = os.getcwd()
# cwd = 'results_4levels_soot_a_500_noClip/'

istart = 0
iskip  = 10
folders=[
        'humidification/',
        ]
for path in folders:

    list_all=fnmatch.filter(os.listdir(path), 'plt*')
    list_all.sort()    
    print("Found output   : {}".format(list_all))

    list_output = []
    for file in list_all:
        if(file.endswith('.png') == False):
            list_output.append(file)

    list_process=list_output[istart:len(list_output):iskip]
    #list_process=[list_process[-1]]
    list_process=['plt252459']
    print("To be processed: {} \n".format(list_process))

    for file in list_process:
        print(file)
        ds = yt.load("{}/{}".format(path,file))

        time = float(ds.current_time)

        if(time < 600.):
            continue

        names = [v[1] for v in ds.field_list]
        print(names)
        # slc.annotate_contour('c')

        #filename=path+'/line_data_t_{:1.2e}.h5'.format(time)
        filename=path+'/{}.h5'.format(file)
        oray = get_line(filename,names)


