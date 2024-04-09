# Known issues for ´cola´ environment
### Connectivity analysis issues in the installation.

As there's many components and libraries involved in this software, your particular machine might present particular configurations that require special attention. Here we recapitulate some of the known issues that you might face. We also provide some solutions for those issues. Feel free to provide your case and if possible, the solution.

We need to evaluate several steps in two stages: the python conda environment and the R Dashboard:

 
#### A. Setting up cola R package and Python conda environment
 
 1. Install ´cola´ R package
 2. Install ´reticulate´ R package
 3. Install miniconda / python software
 4. Install ´cola´ conda environment
 5. Install ´cola´ conda environment packages
 6. Define connections between R and Python
 

#### How to solve each step:

Running ´cola::setup_cola()´ should set up all the different steps and print which of them are successfully finished. Both ´cola::setup_cola()´ and ´cola::commonErrors()´ functions shown the finished steps. See the latest number in the console print and look for that step in the following titles.


######  1. Installing cola R package
  
  Testing: library(cola)
  Expected answer: No answer, just ´>´
  

  Known issue: You might see this message when try to install ´cola´ with the command ´devtools::install_github('connectingLandscapes/cola')´
  
  The message shown is:
  ´devtools::install_github('connectingLandscapes/cola')´
  ´Using GitHub PAT from the git credential store.´
  ´Error: Failed to install 'unknown package' from GitHub:´
  ´HTTP error 401.´
  ´Bad credentials´

  ´Rate limit remaining: 55/60´
  ´Rate limit reset at:´
  
  Solution:  We found this error in old R versions. We solved it installing a new version of R (4.3.3), rtools 4.3  (https://cran.r-project.org/bin/windows/Rtools/), and devtools 2.4.5
  
  We try this links with no success. But 
  https://stackoverflow.com/questions/70908295/failed-to-install-unknown-package-from-github
  https://github.com/r-lib/devtools/issues/1566#issuecomment-320504796 
  https://github.com/scibrokes/setup-centOS7-DO/issues/3
  https://github.com/r-lib/remotes/issues/641


  Error in loadNamespace(x) : there is no package called ‘devtools’
  
  
  --
  
  reticulate::install_miniconda()
reticulate::py_config()
Error in py_module_import(module, convert = convert) : 
  ModuleNotFoundError: No module named 'rpytools'. 
Solution: load reticulate library.  python - No module named 'rpytools'? - Stack Overflow
https://stackoverflow.com/questions/54791126/no-module-named-rpytools

--

reticulate::install_miniconda()
Error: Miniconda is already installed at path "C:/Users/Admin/AppData/Local/r-miniconda".
- Use `reticulate::install_miniconda(force = TRUE)` to overwrite the previous installation.

Error: Miniconda installed already. 
Solution:  force TRUE


Warning messages:
1: In file.copy(savedcopy, lib, recursive = TRUE) :
  problem copying C:\Users\gonza\AppData\Local\R\win-library\4.2\00LOCK\curl\libs\x64\curl.dll to C:\Users\gonza\AppData\Local\R\win-library\4.2\curl\libs\x64\curl.dll: Permission denied
2: In file.copy(savedcopy, lib, recursive = TRUE) :
  problem copying C:\Users\gonza\AppData\Local\R\win-library\4.2\00LOCK\vctrs\libs\x64\vctrs.dll to C:\Users\gonza\AppData\Local\R\win-library\4.2\vctrs\libs\x64\vctrs.dll: Permission denied



###### 2. Install reticulate R package
  
  *Testing:* `library(reticulate)`
  _Expected answer:_ No error
  
  *Known issue:* `Library not installed` after `library(reticulate)`. This indicates the package is not available either because it was not installed or it was a problem during the installation.
  _Solution:_ Run `install.package('reticulate')`. Be sure to having the last RTools (https://cran.r-project.org/bin/windows/Rtools/) and admin privileges if some error arises during the installation. Click on R/Rstudio icon > right click > run as administrator.
  
 
  
###### 3. Install miniconda software
  
  *Testing:* `reticulate::conda_list()`
  _Expected answer:_ Table with local paths to python conda versions
  `     name                                                  python`
  `1    base C:\\Users\\USER\\AppData\\Local\\r-miniconda/python.exe`
  
  Also verify the `base` python path and was installed properly.
  `file.exists( subset(reticulate::conda_list(), name == 'base')$python)`
  This should result in `TRUE`
  
  
  *Known issue:* cola folder existing but not properly configuration
  _Solution:_ Delete folder and install environment again. Go to the path `reticulate::miniconda_path()` and check is not empty. Some other paths, depending under which user you installed miniconda (`reticulate::install_miniconda()`) can include: 
  C:\Users\Admin\AppData\Local\r-miniconda\envs
  C:\Users\Admin\miniconda3\envs
    
  
  *Known issue:*
  _Solution:_
  

  *Testing:* install_miniconda()
  _Expected answer:_

 
  *Known issue:*
  _Solution:_
  
  *Known issue:* cola available but without name
  _Solution:_
  
´Known issue: cola available but without name´
     ´Error -- no name of conda under "conda info --envs"´
  ´Solution: ´
   ´(base) C:\Users\Admin>conda activate C:\Users\Admin\AppData\Local\r-miniconda\envs\cola # activate unnamed env´
   ´conda config --append envs_dirs C:\Users\Admin\AppData\Local\r-miniconda\envs ## add unamed envs´
   ´https://stackoverflow.com/questions/57527131/conda-environment-has-no-name-visible-in-conda-env-list-how-do-i-activate-it-a´


###### 4. Install cola environment packages
 
  *Testing:*
  _Expected answer:_
  
  *Known issue:*
  _Solution:_
  
    2         cola         C:\\Users\\ig299\\AppData\\Local\\r-miniconda\\envs\\cola/python.exe
    
    
    https://stackoverflow.com/questions/54437115/error-networkit-compilation-requires-cmake-288
    

  
  ---
  Error 2: 
Collecting networkit
  Downloading networkit-10.1.tar.gz (5.7 MB)
  Preparing metadata (setup.py): started
  Preparing metadata (setup.py): finished with status 'error'
  ERROR: Command errored out with exit status 1 ………………
from Cython.Build import cythonize
  ModuleNotFoundError: No module named 'Cython'

---

on - ERROR: NetworKit compilation requires cmake. #288 - Stack Overflow

Here is how to fix this error in the main use cases:
    - use 'pip install scikit-learn' rather than 'pip install sklearn'
    - replace 'sklearn' by 'scikit-learn' in your pip requirements files
      (requirements.txt, setup.py, setup.cfg, Pipfile, etc ...)
    - if the 'sklearn' package is used by one of your dependencies,
      it would be great if you take some time to track which package uses
      'sklearn' instead of 'scikit-learn' and report it to their issue tracker
    - as a last resort, set the environment variable
      SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True to avoid this error

--

Error 4: 

 ERROR: NetworKit compilation requires cmake.


###### 5. Define connections between R and Python
 
  *Testing:*
  _Expected answer:_
  
  *Known issue:*
  _Solution:_



  N:\Mi unidad\git\cola
Traceback (most recent call last):
  File "C:/Users/Admin/Documents/R/win-library/4.0/cola/python/welcome.py", line 17, in <module>
  from rasterio.crs import CRS
File "C:\Users\Admin\AppData\Local\R-MINI~1\envs\cola\lib\site-packages\rasterio\__init__.py", line 9, in <module>
  from rasterio._base import gdal_version
ImportError: DLL load failed while importing _base: No se encontrÃ³ el proceso especificado.



  +Step 5/5 Setting up local variables
Warning message:
In system(cmd2test, intern = TRUE) :
  running command 'C:\Users\Admin\AppData\Local\r-miniconda\envs\cola/python.exe C:/Users/Admin/Documents/R/win-library/4.0/cola/python/welcome.py' had status 1
  
  
  
  ---
  > cola::commonErrors()
 
 We found some errors. Running `cola::setup_cola()` should help you to configure the package.
 Please refer to https://github.com/gonzalezivan90/cola/blob/main/known-issues.md for more details
   === All dependencies and requirements installed. Look for futher details in the repository documentation ===
   ---
   

  Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))
  origLibs <- installed.packages()
  Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))
  origLibs <- installed.packages()
  sapply(.libPaths(), list.files, pattern = 'cola', recursive = FALSE)
  
  
#### B. Setting up cola R packages for the dashboard

  *Testing:*
  _Expected answer:_
  
  *Known issue:*
  _Solution:_
  

## Other usfeull comands
>
 "C:/Users/ig299/AppData/Local/r-miniconda/condabin/conda.bat" update --yes --name base conda

