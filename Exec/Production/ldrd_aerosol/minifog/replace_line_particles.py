import os 
import fnmatch
import shutil

def replace_line(filename):
    with open(filename+'/particles/Header',"r") as f:
        newline=[]
        for word in f.readlines():        
            newline.append(word.replace("One","Zero"))  ## Replace the keyword while you copy.
    

    shutil.copyfile(filename+'/particles/Header', filename+'/particles/Header_old')

    with open(filename+'/particles/Header',"w") as f:
        for word in newline:        
            f.write(word)  ## Replace the keyword while you copy.

# Select the time 
istart = 0 # Start index
iskip  = 1 # Skip index
path = os.getcwd()

list_output=[x for x in fnmatch.filter(os.listdir(path), 'plt*') if('.png' not in x and '.csv' not in x)]
list_output.sort()    
print("Found output   : {}".format(list_output))

list_process=list_output[istart:len(list_output):iskip]
print("To be processed: {} \n".format(list_process))


for file in list_process:
    replace_line(path+'/'+file)

