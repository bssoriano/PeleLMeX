
import fnmatch
import shutil
import os.path

# inputs
nprocs=1
istart = 1 # Start index
iskip  = 2 # Skip index

path = os.getcwd()

folders = [path]
print("Folders found:",folders)

for path in folders:
    print("Reading files from: ", path)
    all_files=fnmatch.filter(os.listdir(path), 'chk*')
    all_files.sort()    

    #remove unwanted files
    list_output = []
    for filename in all_files:
        if(filename.endswith(".png") or filename.endswith(".tar.gz") or filename.endswith(".csv")):
            #or filename == "chk00000" or filename == "chk02000" or filename == "chk05000" or filename == "chk10000" or filename == 'chk12000'
            #or filename == "chk15000" or filename == "chk20000" or filename == "chk30000" or filename == "chk33300"):
    #        or filename == "chk15000" or filename == "chk20000"):  
            pass
        else:

            list_output.append(filename)

    print('All files found: ',len(list_output))
    list_output=list_output[istart:len(list_output)-1:iskip]
    # list_output=list_output[:-2]
    print("To be deleted   : {}".format(len(list_output)))

    print("Proceed? (yes/no)")
    proceed = input()

    if(proceed == 'yes'):
        print('Deleting files...')
        for case in list_output:
            shutil.rmtree(path+'/'+case)
            # os.remove(case)
    else:
        print("Answer was {}. Exiting...".format(proceed))
