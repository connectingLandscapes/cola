# Known issues for `cola` environment
### Connectivity analysis issues in the installation.

As there's many components and libraries involved in this software, your particular machine might present particular configurations that require special attention. Here we recapitulate some of the known issues that you might face. We also provide some solutions for those issues. *Feel free to provide your case and if possible, the solution.*

We need to evaluate several steps in two stages: the python conda environment and the R Dashboard:


#### A. Setting up cola R package and Python conda environment
 
 1. Install `cola` R package
 2. Install `reticulate` R package
 3. Install miniconda / python software
 4. Install `cola` conda environment
 5. Install `cola` conda environment packages
 6. Define connections between R and Python
 

#### How to solve each step:

Running `cola::setup_cola()` should set up all the different steps and print which of them are successfully finished. Both `cola::setup_cola()` and `cola::commonErrors()` functions shown the finished steps. See the latest number in the console print and look for that step in the following titles.


######  1. Installing cola R package
  During the instruction `devtools::install_github('connectingLandscapes/cola')` or after the installation, when `library(cola)`:
  
  *Testing:* `library(cola)`
  _Expected answer:_ No answer, just `>`
  

  *Known issue:* You might see this message when try to install `cola` with the command `devtools::install_github('connectingLandscapes/cola')`
  
  The message shown is:
  `devtools::install_github('connectingLandscapes/cola')`
  `Using GitHub PAT from the git credential store.`
  `Error: Failed to install 'unknown package' from GitHub:`
  `HTTP error 401.`
  `Bad credentials`

  `Rate limit remaining: 55/60`
  `Rate limit reset at:`
  
  _Solution:_  We found this error in old R versions. We solved it installing a new version of R (4.3.3), rtools 4.3  (https://cran.r-project.org/bin/windows/Rtools/), and devtools 2.4.5
  
  We try this links with no success. But 
  https://stackoverflow.com/questions/70908295/failed-to-install-unknown-package-from-github
  https://github.com/r-lib/devtools/issues/1566#issuecomment-320504796 
  https://github.com/scibrokes/setup-centOS7-DO/issues/3
  https://github.com/r-lib/remotes/issues/641



  *Known issue:* When installing the package,  `devtools::install_github('connectingLandscapes/cola')`, you might see:
  `Error in loadNamespace(x) : there is no package called ‘devtools’`
  
  _Solution:_ Install and load `devtools`. Use `install.packages('devtools')` for this. Consider using the latest version of R and also install RTools (excecutable, not an R library)
  
  
  
  *Known issue:* `Error in py_module_import(module, convert = convert) :` `ModuleNotFoundError: No module named 'rpytools'. `
  
  _Solution:_  Load   `reticulate` library: `library(reticulate)`
  https://stackoverflow.com/questions/54791126/no-module-named-rpytools



  *Known issue:* Getting the message:
  `problem copying C:\Users\...\AppData\Local\R\win-library\4.2\00LOCK\'some.file' to C:\Users\...\AppData\Local\R\win-library\4.2\'other.file': Permission denied`.
  The libraries and paths might change, but still having the `00LOCK` text.
  
  _Solution:_ Be sure to close all R sessions that might be using or writing the folder with the libraries. Then, try the command again. If the error persist, when closing all the R sessions, check if the library folder (found them by typing: `.libPaths()`) still have the lock file or folder, and remove it manually. 



###### 2. Install reticulate R package
  
  Installing `reticulate` R package should be done by the `cola::setup_cola()` function. Here some issues detected:
  
  *Testing:* `library(reticulate)`
  _Expected answer:_ No error
  
  *Known issue:* `Library not installed` after `library(reticulate)`. This indicates the package is not available either because it was not installed or it was a problem during the installation.
  _Solution:_ Run `install.package('reticulate')`. Be sure to having the last RTools (https://cran.r-project.org/bin/windows/Rtools/) and admin privileges if some error arises during the installation. Click on R/Rstudio icon > right click > run as administrator.
  
  

###### 3. Install miniconda software
  
  Installing `miniconda` should be done by the `cola::setup_cola()` function. There's several reported problems trying to install miniconda, so we encourage not to uninstall it once is installed. Here some issues detected:


  *Testing:* `reticulate::conda_list()`
  _Expected answer:_ Table with local paths to python (mini) conda versions
  `     name                                                  python`
  `1    base C:\\Users\\USER\\AppData\\Local\\r-miniconda/python.exe`
  
  Also verify the `base` python path and was installed properly.
  `file.exists( subset(reticulate::conda_list(), name == 'base')$python)`
  This should result in `TRUE`
  
  
  *Known issue:* *Known issue:*  `Error: Miniconda is already installed at path "C:/Users/Admin/AppData/Local/r-miniconda".`
  _Solution:_ Run `reticulate::install_miniconda(force = TRUE)`
  
  
  *Known issue:* Getting `Error: Miniconda installation failed [unknown reason]`. This might result from a broken installation, so some files and folders can exists on your machine, but not completed.
  `> reticulate::install_miniconda( )`
  ` running reticulate::install_miniconda( )`
  ` * Installing Miniconda -- please wait a moment ...`
  ` * Downloading "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" ...`
  ` trying URL 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe'`
  ` Content type 'application/octet-stream' length 74687656 bytes (71.2 MB)`
  ` downloaded 71.2 MB `
  `Error: Miniconda installation failed [unknown reason]`

  _Solution:_ There's many reasons. We solved the problem by:
  - Finding the uninstall executable of miniconda `Uninstall-Miniconda3.exe`, typically in the folders given by `reticulate::miniconda_path()` or `reticulate::conda_list()`. If those commands throw errors, try on these paths:
     + `C:/Users/Admin/miniconda3`
     + `C:/Users/Admin/AppData/Local/r-miniconda`
  - Use the .exe file to uninstall conda
  - Try `reticulate::install_miniconda( )` again
    https://github.com/rstudio/reticulate/issues/1297
    
    
  Also, this post suggests to:
    `remotes::install_github("hafen/rminiconda")`
    `rminiconda::install_miniconda(name='cola')`
    `py <- rminiconda::find_miniconda_python("cola")`
    `reticulate::use_python(py, required = TRUE)`
  
  
  
###### 4. Install cola conda environment
Installing `cola` conda envronment should be done by the `cola::setup_cola()` function. Here some issues detected:

  *Testing:* `reticulate::conda_list()`
  _Expected answer:_ Finding a `cola` conda environment in the resulting table
```  
          name                                                                       python
1         base                     C:\\Users\\USER\\AppData\\Local\\r-miniconda/python.exe
2         cola         C:\\Users\\USER\\AppData\\Local\\r-miniconda\\envs\\cola/python.exe
3 r-reticulate C:\\Users\\USER\\AppData\\Local\\r-miniconda\\envs\\r-reticulate/python.exe
```

and having a message from `cola::setup_cola()` similar to: `The python version is Python 3.9.19`

  *Known issue:* `cola` folder existing but not properly configurated. 
  _Solution:_ Delete folder and install environment again. Go to the path `reticulate::miniconda_path()` and check is not empty. Some other paths, depending under which user you installed miniconda (`reticulate::install_miniconda()`) can include: 
  - `C:\Users\USER\AppData\Local\r-miniconda\envs\cola`
  - `C:\Users\USER\miniconda3\envs\cola`


  
  *Known issue:* `cola` conda environment available but without name. Your conda manager installed the environment correctly but is not named. This might occur because the path where R install the environment is different than the manager conda path. Likely occur when changing users, permissions. 
  `(base) C:\Users\Admin>conda activate C:\Users\USER\AppData\Local\r-miniconda\envs\cola`
   `Error -- no name of conda under "conda info --envs"`
   
   After running `reticulate::conda_list()` or in the conda prompt `conda info --envs`:
  
```
base                  *  C:\ProgramData\Miniconda3
                         C:\Users\Admin\AppData\Local\R-MINI~1\envs\r-reticulate
                         C:\Users\Admin\AppData\Local\r-miniconda\envs\cola
```

  
  _Solution:_ In your conda console (outside R or R studio) and add the conda enviroment folder path to the conda manager paths.
   `(base) conda config --append envs_dirs C:\Users\Admin\AppData\Local\r-miniconda\envs`
   
   https://stackoverflow.com/questions/57527131/conda-environment-has-no-name-visible-in-conda-env-list-how-do-i-activate-it-a




###### 5. Install cola conda environment packages
Installing the `cola` conda environment packages should be done by the `cola::setup_cola()` function. Here some issues detected:

  *Testing:* Running `cola::setup_cola()` and `cola::commonErrors()` getting the following messages:
  _Expected answer:_ `All required modules installed!` or `4. All 'cola' conda environment packages installed`


  *Known issue:* Some Windows versions require the last VisualStudio version for Cython:
  `Collecting networkit`
  `Downloading networkit-10.1.tar.gz (5.7 MB)`
  `Preparing metadata (setup.py): started`
  `Preparing metadata (setup.py): finished with status 'error'`
  `ERROR: Command errored out with exit status 1`
  `from Cython.Build import cythonize`
  `ModuleNotFoundError: No module named 'Cython'`
  
  or similar messages:
  `ERROR: NetworKit compilation requires cmake.`
  `error: Microsoft Visual C++ 14.0 or greater is required`
  
  
  _Solution:_ Follow this (tutorial) [https://stackoverflow.com/questions/64261546/how-to-solve-error-microsoft-visual-c-14-0-or-greater-is-required-when-inst], install the `vs_BuildTools.exe` file checking the 'C++ build tools'


###### 6. Define connections between R and Python
These errors are common because each computer has a particular configuration. Here some of the solutions we found.
 
  *Testing:* Running `cola::setup_cola()` and `cola::commonErrors()` getting the following messages:
  _Expected answer:_ `All required modules installed!` or `=== All dependencies and requirements installed`
  
  *Testing:* Also ypu need to have this two internal variables defines: `Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))`
  _Expected answer:_
  `                                        COLA_MINICONDA_PATH                                                      COLA_SCRIPTS_PATH `
`"C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe"  "C:/Users/ig299/AppData/Local/Programs/R/R-4.3.3/library/cola/python"`
  

  *Known issue:* Some python libraries have not correctly defined internal paths to found other modules. Some errors include something like:
  `from rasterio._version import gdal_version, get_geos_version, get_proj_version"`
  `"ImportError: DLL load failed while importing _version`
  
  This occurs because even when the correct version of python is installed (step 4), the internal paths and libraries installed in the step 5 are not found.
  
  _Solution:_ Use the 'base' conda environment instead of 'cola'. For this, we need to: 
  1. Follow this instruction to activate conda by default in the command line: 
    As Admin in powershell: set-executionpolicy unrestricted
  2. Go to `%USERPROFILE%\Miniconda3\Scripts` or something like `C:\Users\Admin\miniconda3\Scripts`, using `cd` command: `cd C:\Users\Admin\miniconda3\Scripts`
  3. Run `conda --version` to check the version
  4. Run `conda init powershell` to activate conda
  5. Run again `cola::setup_cola(envname = 'base')`
  

Bug / issue reporting in cola R package
https://docs.google.com/forms/d/e/1FAIpQLSdFsM1e02biuauaWE4Svwtu5QMneKU7Ilfa8pAJHiRy3a-KGw/viewform?usp=sf_link

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

