# Known issues for ´cola´ environment


## Connectivity analysis issues in the installation.
 
 We need to sort the following steps:
 1. Install reticulate python package
 2. Install miniconda software
 3. Install cola environment
 4. Install cola environment packages
 5. Define connections between R and Python
 
 
## How to solve each step:
 
 1. Install reticulate R package
  
  Testing: `library(reticulate)`  
  Expected answer: No error
  
  Known issue: Library not installed
  Solution: `install.package('reticulate')`. Be sure to have RTools and admin privilegies if some error arrise during the installation.
  
  
  https://github.com/r-lib/remotes/issues/641
  

 2. Install cola R package
  
  Testing: 
  Expected answer:
  
  Known issue: 
  devtools::install_github('gonzalezivan90/cola')
Using GitHub PAT from the git credential store.
Error: Failed to install 'unknown package' from GitHub:
  HTTP error 401.
  Bad credentials

  Rate limit remaining: 55/60
  Rate limit reset at: 
  Solution:  
  
  https://github.com/r-lib/devtools/issues/1566#issuecomment-320504796 
  https://github.com/scibrokes/setup-centOS7-DO/issues/3
  
  Known issue: 
  Solution:



  # Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))
  # origLibs <- installed.packages()
  sapply(.libPaths(), list.files, pattern = 'cola', recursive = FALSE)

  
 2. Install miniconda software
  
  Testing: `reticulate::conda_list()`
  Expected answer: Table with local paths to python conda versions
  `     name                                                  python`
  `1    base C:\\Users\\USER\\AppData\\Local\\r-miniconda/python.exe`
  
  Known issue: cola folder existing but not properly configurated
  Solution: delete folder and install environment again
  
  Known issue: 
  Solution:
  
 3. Install cola environment
  
  Testing: install_miniconda()
  Expected answer:
  2         cola         C:\\Users\\ig299\\AppData\\Local\\r-miniconda\\envs\\cola/python.exe

 

  
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
