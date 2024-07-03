# Known issues for `cola` installation

As there's many components and libraries involved in this software, your particular machine might present particular configurations that require special attention. Here we recapitulate some of the known issues that you might face. We also provide some solutions for those issues. **Feel free to provide your case and if possible, the solution.**

We need to evaluate several steps in two stages: ***A) the python conda environment*** and ***B) the R Dashboard***


#### **A. Setting up cola R package and Python conda environment**
 
 1. Install `cola` R package
 2. Install `reticulate` R package
 3. Install miniconda / python software
 4. Install `cola` conda environment
 5. Install `cola` conda environment packages
 6. Define connections between R and Python
 

##### How to solve each step:

Running `cola::setup_cola()` should set up all the different steps and print which of them are successfully finished. Both `cola::setup_cola()` and `cola::commonErrors()` functions shown the finished steps. See the latest number in the console print and look for that step in the following titles.




-------------
-------------

####  **1. Installing cola R package**
  During the instruction `devtools::install_github('connectingLandscapes/cola')` or after the installation, when `library(cola)`:
  
  **Testing:** `library(cola)`
  
  ***Expected answer:*** No answer, just `>`
  
-------------


  **Known issue:** You might see this message when try to install `cola` with the command `devtools::install_github('connectingLandscapes/cola')`
  
  The message shown is:
  ```
  devtools::install_github('connectingLandscapes/cola')
  Using GitHub PAT from the git credential store.
  Error: Failed to install 'unknown package' from GitHub:
  HTTP error 401.
  Bad credentials

  Rate limit remaining: 55/60
  Rate limit reset at:
  ```
  
  ***Solution:***  We found this error in old R versions. We solved it installing a new version of R (4.3.3), [rtools 4.3]  (https://cran.r-project.org/bin/windows/Rtools/), and devtools 2.4.5
  
  
  Here more resources about the issue.  
  https://stackoverflow.com/questions/70908295/failed-to-install-unknown-package-from-github
  
  https://github.com/r-lib/devtools/issues/1566#issuecomment-320504796 
  
  https://github.com/scibrokes/setup-centOS7-DO/issues/3
  
  https://github.com/r-lib/remotes/issues/641


-------------


  **Known issue:** When installing the package,  `devtools::install_github('connectingLandscapes/cola')`, you might see:
  `Error in loadNamespace(x) : there is no package called ‘devtools’`
  
 ***Solution:*** Install and load `devtools`. Use `install.packages('devtools')` for this. Consider using the latest version of R and also install RTools (excecutable, not an R library)


-------------


  **Known issue:** `Error in py_module_import(module, convert = convert) :` `ModuleNotFoundError: No module named 'rpytools'. `
  
 ***Solution:***  Load   `reticulate` library: `library(reticulate)`. 
 [Here a reference](https://stackoverflow.com/questions/54791126/no-module-named-rpytools)


-------------


  **Known issue:** Getting the message:
  `problem copying C:\Users\...\AppData\Local\R\win-library\4.2\00LOCK\'some.file' to C:\Users\...\AppData\Local\R\win-library\4.2\'other.file': Permission denied`.
  The libraries and paths might change, but still having the `00LOCK` text.
  
 ***Solution:*** Be sure to close all R sessions that might be using or writing the folder with the libraries. Then, try the command again. If the error persist, when closing all the R sessions, check if the library folder (found them by typing: `.libPaths()`) still have the lock file or folder, and remove it manually. 


-------------


 **Known issue:** Getting the message: `there is no package called ‘cola’`. You probably tried to install cola already, so let's check is there's a broken installation.
  
 ***Solution:*** Try `remove.packages('cola')` to remove previous installations, or checBe sure to close all R sessions that might be using or writing the folder with the libraries. Then, try the command again. If the error persist, when closing all the R sessions, check if the library folder (found them by typing: `.libPaths()` or `sapply(.libPaths(), list.files, pattern = 'cola', recursive = FALSE)`) still have the lock file or folder, and remove it manually.
 
    -------------


 **Known issue:** Can't install the package, probbaly for a interruption in the internet connection or was interrupted by a process in your system. You can get the message: 
 
 ```
 Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
  download from 'https://api.github.com/repos/connectingLandscapes/cola/tarball/HEAD' failed
  ```
  
  
 ***Solution:*** Try again `devtools::install_github('connectingLandscapes/cola'`. If the error persists, restart your R session. Here the log of the console of the reported case and solution: 
 ```
> devtools::install_github('connectingLandscapes/cola')  # <----- HERE first try 
Downloading GitHub repo connectingLandscapes/cola@HEAD
Error in utils::download.file(url, path, method = method, quiet = quiet,  : 
  download from 'https://api.github.com/repos/connectingLandscapes/cola/tarball/HEAD' failed
> devtools::install_github('connectingLandscapes/cola') # <----- HERE second try, working, with no changes
Downloading GitHub repo connectingLandscapes/cola@HEAD
── R CMD build ───────────────────────────────────────────────────────────────────────────────────────
✔  checking for ... 

```

-------------
-------------

####  **2. Install `reticulate` R package**
  
  Installing `reticulate` R package should be done by the `cola::setup_cola()` function. Here some issues detected:
  
  **Testing:** `library(reticulate)`
 
 ***Expected answer:*** No error
 

-------------



  **Known issue:** `Library not installed` after `library(reticulate)`. This indicates the package is not available either because it was not installed or it was a problem during the installation.
 
 ***Solution:*** Run `install.package('reticulate')`. Be sure to having the last [RTools](https://cran.r-project.org/bin/windows/Rtools/) and admin privileges if some error arises during the installation. Click on R/Rstudio icon > right click > run as administrator.
  
  



####  **3. Install miniconda software**
  
  Installing `miniconda` should be done by the `cola::setup_cola()` function. There's several reported problems trying to install miniconda, so we encourage not to uninstall it once is installed. Here some issues detected:


  **Testing:** `reticulate::conda_list()`
 
 ***Expected answer:*** Table with local paths to python (mini) conda versions
 
 ```  
         name                                                  python
         
  1    base C:\\Users\\USER\\AppData\\Local\\r-miniconda/python.exe
  ```  
  
  Also verify the `base` python path and was installed properly.
  `file.exists( subset(reticulate::conda_list(), name == 'base')$python)`
  This should result in `TRUE`
  

-------------



  
  **Known issue:**  `Error: Miniconda is already installed at path "C:/Users/Admin/AppData/Local/r-miniconda".`
 
 ***Solution:*** Run `reticulate::install_miniconda(force = TRUE)`
  

-------------



  
  **Known issue:** Getting `Error: Miniconda installation failed [unknown reason]`. This might result from a broken installation, so some files and folders can exists on your machine, but not completed.
  
  ```
  > reticulate::install_miniconda( )
  
   running reticulate::install_miniconda( )
   ** Installing Miniconda -- please wait a moment ...
   ** Downloading "https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe" ...
   trying URL 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe'
   Content type 'application/octet-stream' length 74687656 bytes (71.2 MB)
   downloaded 71.2 MB 
  Error: Miniconda installation failed [unknown reason]
  ```


 ***Solution:*** There's many reasons. We solved the problem by:
  - Finding the uninstall executable of miniconda `Uninstall-Miniconda3.exe`, typically in the folders given by `reticulate::miniconda_path()` or `reticulate::conda_list()`. If those commands throw errors, try on these paths:
     + `C:/Users/Admin/miniconda3`
     + `C:/Users/Admin/AppData/Local/r-miniconda`
  - Use the .exe file to uninstall conda
  - Try `reticulate::install_miniconda( )` again. Here a [reference](https://github.com/rstudio/reticulate/issues/1297)
    
  Also, this post suggests to:
    `remotes::install_github("hafen/rminiconda")`
    `rminiconda::install_miniconda(name='cola')`
    `py <- rminiconda::find_miniconda_python("cola")`
    `reticulate::use_python(py, required = TRUE)`
  
  
  
-------------


  **Known issue:**  `Error: Miniconda is already installed at path "C:/Users/Admin/AppData/Local/r-miniconda".`
 
  ***Solution:*** Run `reticulate::install_miniconda(force = TRUE)`
  

-------------

  


  **Known issue:**  Get the message `Error: Unable to find conda binary. Is Anaconda installed?`
 
  ***Solution:*** Find where `conda` or `conda.bat` lives. Add that path to the R command `options(reticulate.conda_binary = "C:/path/to/conda.bat")`, or the equivalent.
  The paths might include: `C:/Users/Admin/AppData/Local/r-miniconda/condabin/conda.bat` or `C:/Users/Admin/miniconda3/condabin/conda.bat` among others.
  

-------------
-------------
  
  
####  **4. Install cola conda environment**
Installing `cola` conda environment should be done by the `cola::setup_cola()` function. Here some issues detected:

  **Testing:** `reticulate::conda_list()`
  
 ***Expected answer:*** Finding a `cola` conda environment in the resulting table
```  
          name                                                                       python
1         base                     C:\\Users\\USER\\AppData\\Local\\r-miniconda/python.exe
2         cola         C:\\Users\\USER\\AppData\\Local\\r-miniconda\\envs\\cola/python.exe
3 r-reticulate C:\\Users\\USER\\AppData\\Local\\r-miniconda\\envs\\r-reticulate/python.exe
```
and having a message from `cola::setup_cola()` similar to: `The python version is Python 3.9.19`



-------------



  **Known issue:** `cola` folder existing but not properly configurated. 
  
 ***Solution:*** Delete folder and install environment again. Go to the path `reticulate::miniconda_path()` and check is not empty. Some other paths, depending under which user you installed miniconda (`reticulate::install_miniconda()`) can include: 
  - `C:\Users\USER\AppData\Local\r-miniconda\envs\cola`
  - `C:\Users\USER\miniconda3\envs\cola`
  
  [Here a reference](https://stackoverflow.com/questions/60974507/error-unable-to-find-conda-binary-is-anaconda-installed)



-------------



  **Known issue:** `cola` environment can't be installed
  
 ***Solution:*** Install conda environment manually:
 `conda_install(envname = 'cola', packages = c('gdal', 'h5py', 'numexpr', 'rasterio', 'pytables', 'pandas',  'cython', 'numba',                   'networkit', 'fiona', 'shapely', 'geopandas', 'kdepy','scikit-image'))` 


-------------



  
  **Known issue:** `cola` conda environment available but without name. Your conda manager installed the environment correctly but is not named. This might occur because the path where R install the environment is different than the manager conda path. Likely occur when changing users, permissions. 
  `(base) C:\Users\Admin>conda activate C:\Users\USER\AppData\Local\r-miniconda\envs\cola`
   `Error -- no name of conda under "conda info --envs"`
   
   After running `reticulate::conda_list()` or in the conda prompt `conda info --envs`:
  
```
base                  **  C:\ProgramData\Miniconda3
                         C:\Users\Admin\AppData\Local\R-MINI~1\envs\r-reticulate
                         C:\Users\Admin\AppData\Local\r-miniconda\envs\cola
```


  
 ***Solution:*** In your conda console (outside R or R studio) add the conda environment folder path to the conda manager paths.
 `(base) conda config --append envs_dirs C:\Users\Admin\AppData\Local\r-miniconda\envs`, with the equivalent path in your computer.
   
   Here a [reference](https://stackoverflow.com/questions/57527131/conda-environment-has-no-name-visible-in-conda-env-list-how-do-i-activate-it-a)


-------------
-------------

####  **5. Install cola conda environment packages**
Installing the `cola` conda environment packages should be done by the `cola::setup_cola()` function. Here some issues detected:

  **Testing:** Running `cola::setup_cola()` and `cola::commonErrors()` getting the following messages:
  
 ***Expected answer:*** `All required modules installed!` or `4. All 'cola' conda environment packages installed`



-------------


  **Known issue:** Some Windows versions require the last VisualStudio version for Cython:
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
  
  
 ***Solution:*** Follow this [tutorial](https://stackoverflow.com/questions/64261546/how-to-solve-error-microsoft-visual-c-14-0-or-greater-is-required-when-inst), install the `vs_BuildTools.exe` file checking the 'C++ build tools'


-------------
-------------

####  **6. Define connections between R and Python**
These errors are common because each computer has a particular configuration. Here some of the solutions we found.
 
  **Testing:** Running `cola::setup_cola()` and `cola::commonErrors()` getting the following messages:
  
 ***Expected answer:*** `All required modules installed!` or `=== All dependencies and requirements installed`
  

  **Testing:** Also you need to have this two internal variables defines: `Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))`
  
 ***Expected answer:***
 
 ```
 COLA_MINICONDA_PATH
"C:\Users\USER\AppData\Local\r-miniconda\envs\cola\python.exe"  
COLA_SCRIPTS_PATH 
"C:/Users/ig299/AppData/Local/Programs/R/R-4.3.3/library/cola/python"
```  
  

-------------


  **Known issue:** Some python libraries have not correctly defined internal paths to found other modules. Some errors include something like:
  `from rasterio._version import gdal_version, get_geos_version, get_proj_version"`
  `"ImportError: DLL load failed while importing***version`
  
  This occurs because even when the correct version of python is installed (step 4), the internal paths and libraries installed in the step 5 are not found.
  
 ***Solution:*** Use the 'base' conda environment instead of 'cola'. For this, we need to: 
  1. Follow these [instructions](https://gist.github.com/martinsotir/2bd2e16332dff71e0fa5be3ed3468a6c) to activate conda by default in the command line: 
    As Admin in powershell: `set-executionpolicy unrestricted`
  2. Go to `%USERPROFILE%\Miniconda3\Scripts` or something like `C:\Users\Admin\miniconda3\Scripts`, using `cd` command: `cd C:\Users\Admin\miniconda3\Scripts`
  3. Run `conda --version` to check the version
  4. Run `conda init powershell` to activate conda
  5. Run again `cola::setup_cola(envname = 'base')`
  
  
  -------------


  **Known issue:** The `cola::cola_setup () `log doesn't print `=== Ready to connect landscapes! ===`. We can't guess the appropriate way to run `cola`.
    This occurs because even when the correct version of python is installed (step 4), the internal paths and libraries installed in the step are not finished. In a new R session, the following variables shouldn't be empty:
  ```
  Sys.getenv("COLA_SCRIPTS_PATH")
  [1] "C:/Users/Admin/Documents/R/win-library/4.0/cola/python"
  
  Sys.getenv("COLA_PYTHON_PATH")
    Case A (most of the cases):
  [1] "C:/Users/Admin/AppData/Local/r-miniconda/envs/cola/python.exe"
    Case B:
  [1] "C:/Users/Admin/AppData/Local/r-miniconda/condabin/conda.bat run --cwd C:/Users/Admin/Documents/R/win-library/4.0/cola/python -n cola python "
  ```
  
 ***Solution:*** Set up the manually the environmental variables. This requires that previous steps are completes: conda environment installed with all the libraries. We need R connect Python properly, so configure R internally.
 
      
    1. Find the `welcome.py` file and it's path: 
    ```
      (welcomepy <- system.file("python/welcome.py", package = "cola"))
      (cola_scripts_path <- dirname(welcomepy))
    ```
    
    2. Find the `python.exe` file
    
    ```
      (pyCola <- subset(reticulate::conda_list(), name == 'cola')$python )
    ```
    
    
    3. Run the python file in R. Try this:
    ```
      ## Check the command line
      cat(cmd <- paste( pyCola, welcomepy) )
      
      system(cmd, intern = TRUE)
      # Some errors can include "ModuleNotFoundError: No module named '_gdal'"
    ```
    
    The result must include  `=== Ready to connect landscapes! ===`.
    
    In some computers this might no work. If is your case, use the following (adatp the paths to your case);
    
    ```
      
      ## Option 1;
      (pyCola <- 'conda run -n cola python ')
      cat(cmd <- paste0(pycola, welcomepy) )
      system(cmd, intern = TRUE)
      # One possible error is "'conda' not found"
      
      ## Option 2 (keep the last space at the end):
      (pyCola <- paste0(reticulate::conda_binary(), ' run --cwd ', cola_scripts_path, ' -n ', envName,' python '))
      cat(cmd <- paste0(pyCola, welcomepy)))

      system(cmd, intern = TRUE)
    ``` 
    
    4. If the previous command worked, saved the paths manually using `pyCola` as `COLA_PYTHON_PATH` and `cola_scripts_path` as  `COLA_SCRIPTS_PATH` depending on the way worked for you previously.
    
    ``` 
    COLA_PYTHON_PATH="C:/Users/Admin/AppData/Local/r-miniconda/envs/cola/python.exe"
    COLA_SCRIPTS_PATH="C:/Users/Admin/Documents/R/win-library/4.0/cola/python"
    ``` 
    
    Save those paths permanently in your system using one of the following options:
    - System environment variables following this [instructions](
      https://stackoverflow.com/questions/69550830/init-fs-encoding-failed-to-get-the-python-codec-of-the-filesystem-encoding)
      
    - Edit the `.Renviron` file located in `Sys.getenv("HOME")` and add the lines mentioned before. In my computer the file is lcoated in `C:\Users\Admin\Documents\.Renviron` at the file content looks like:
     ``` 
     COLA_PYTHON_PATH="C:/Users/Admin/AppData/Local/r-miniconda/condabin/conda.bat run --cwd C:/Users/Admin/Documents/R/win-library/4.0/cola/python -n cola python "
     COLA_SCRIPTS_PATH="C:/Users/Admin/Documents/R/win-library/4.0/cola/python"
    ``` 
    


-------------
-------------

####  **7. Installing and running packages for the DSS dashboard**
The installation of the required packages for `cola` dashboard  should be done by the `cola::setup_cola_dss()` function. Here some issues detected:

  **Testing:** Running `cola::setup_cola_dss()` and getting the following messages:
  
 ***Expected answer:*** `=== All libraries required for COLA's DSS installed ===`


-------------



  **Known issue:** Some packages versions requires to be updated. The package migth exists in your machine, but is not usable. Some of these errors might include:
  
  
  `Warning: Error in : 'addRasterImage()' for SpatRaster objects requires {terra} 1.6-3 or higher`
  `Iinherits(x, "RasterLayer") is not TRUE`
  
  
 ***Solution:*** Be sure to use the last version of R. And install manually the package the log points out, like `install.package('terra')` 


-------------
-------------

####  **Check your local paths**
 You can also try this command to see what's your cola path: 
 
  `Sys.getenv(c('COLA_MINICONDA_PATH', 'COLA_SCRIPTS_PATH'))`

  Those paths must point to two different routes containing the python and scripts paths


#### **B. Setting up cola R packages for the dashboard**
The errors in this section are usually related with some libraries which require manual installations

  **Testing:** Running `cola::setup_cola_dss()`
  
 ***Expected answer:*** `=== All libraries required for COLA's DSS installed ===`
  
  

-------------



  **Known issue:** `Error in library(tidyverse) : there is no package called ‘tidyverse’`
  
 ***Solution:*** `install.packages('tidyverse')`. Adapt to the corresponding library. First try saying _NO_ to install dependencies that requires compilation. If error persists hit _YES_
  
  
-------------



  **Known issue:**  R tries to install a package that requires another dependeny. Install that dependency first.
  `Error in loadNamespace(j <- i[[1L]], c(lib.loc, .libPaths()), versionCheck = vI[[j]]) :`
  `namespace 'fastmap' 1.1.0 is already loaded, but >= 1.1.1 is required`
  
 ***Solution:*** Restart R. Then install the dependency mentioned in the message: `install.packages('fastmap')`. 
    
-------------
    
  **Known issue:** Here your contribution!
  
 ***Solution:*** Here your contribution!
  

-------------
-------------
  
#### **Submit your issue**
  
Consider report an issue on GitHub, or submit your case [here](https://docs.google.com/forms/d/e/1FAIpQLSdFsM1e02biuauaWE4Svwtu5QMneKU7Ilfa8pAJHiRy3a-KGw/viewform?usp=sf_link)



#### **Other usfeull comands**

Command for installing conda in Windows:

***Install cola conda environment*** 

```
conda create -n cola gdal h5py numexpr rasterio pytables pandas cython numba networkit fiona shapely geopandas scikit-image -c conda-forge
```


***List available conda environments***
```
conda info --envs
```

***Remove conda environments***
```
# # conda remove --name cola --all
```


