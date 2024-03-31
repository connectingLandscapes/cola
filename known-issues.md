# Known issues for ´cola´ environment


## Connectivity analysis issues in the installation.
 
 We need to sort the following steps:
 1. Install reticulate python package
 2. Install miniconda software
 3. Install cola environment
 4. Install cola environment packages
 5. Define connections between R and Python
 
 
## How to solve each step:
 
 1. Install reticulate python package
  
  Testing:
  Expected answer:
  
  Known issue:
  Solution:
  
 2. Install miniconda software
  
  Testing:
  Expected answer:
  
  Known issue: cola folder existing but not properly configurated
  Solution: delete folder and install environment again
  
  Known issue: 
  Solution:
  
 3. Install cola environment
  
  Testing:
  Expected answer:
  
  Known issue: cola available but without name
     Error -- no name of conda under "conda info --envs"
  Solution: 
   (base) C:\Users\Admin>conda activate C:\Users\Admin\AppData\Local\r-miniconda\envs\cola # activate unnamed env
   conda config --append envs_dirs C:\Users\Admin\AppData\Local\r-miniconda\envs ## add unamed envs
   https://stackoverflow.com/questions/57527131/conda-environment-has-no-name-visible-in-conda-env-list-how-do-i-activate-it-a
   
   
 4. Install cola environment packages
 
  Testing:
  Expected answer:
  
  Known issue:
  Solution:
  
 5. Define connections between R and Python
 
  Testing:
  Expected answer:
  
  Known issue:
  Solution:
