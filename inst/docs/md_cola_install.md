# Installing Connecting landscapes `cola` 

Non-expert code users can follow [this tutorial](https://docs.google.com/presentation/d/1d2TgZGqDut8_HRP-xZnRnmA3W0OzP3-U/edit#slide=id.p1)

Find on this page the details and instructions for **a)** installation and **b)** customization *cola* installation parameters.


This code contains the quick-and-dirt way to install `cola`.  For the full procedure from scratch check the details in the next session.

```
## Windows users must have git or github installed.

if(!require(devtools)){install.packages('devtools')}
library(devtools) # Library for installing packages from github
devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never') ## Type "3": None ## Installs cola

## If error: Go to https://github.com/connectingLandscapes/cola >> Green button >>
## Download zip >> D:/path/to/cola-main.zip. Then run the following commands:
# install.packages("D:/path/to/cola-main.zip", repos = NULL, type = "win.binary")
## or
# devtools::install_local('D:/path/to/cola-main.zip')

library(cola) # Load cola
cola::setup_cola() # Setup cola. Run this line until all problems are solved.
# Expect 2 pop up windows for accepting to install Miniconda and cola environment

cola::diagnose_cola() # check installation problems -- also suggests solutions
cola::setup_cola() # Run again until solving all installation issue. Required to be
#   finished in order to set up the package properly. 

file.edit(file.path(Sys.getenv("HOME"), ".Renviron")) # Edit parameters
.rs.restartR() # Restart RStudio
Sys.getenv(c('COLA_PYTHON_PATH', 'COLA_SCRIPTS_PATH')) # Validate DSS parameters
Sys.getenv(c('COLA_DATA_PATH', 'COLA_SCRIPTS_PATH', 'COLA_DSS_UPL_MB', 'COLA_VIZ_THRES_PIX', 'COLA_VIZ_RES_NCOL', 'COLA_VIZ_RES_NROW', 'COLA_NCORES' ))

## Install the front end. Windows users must have Rtools and git installed.
cola::setup_cola_dss() # Setup the DSS libraries. Takes a while. Run until getting a success message
.rs.restartR() # Restart R
library(cola) # Load Cola
cola::cola_dss() # Launch the DSS

```




### A) Installation 

This package integrates R and Python modules. 

It's required to install several components (once). The structure of this software is:

- R as the base. Use the latest available [here](https://cran.r-project.org/bin/windows/base/)

- Git as R package mange; Do you have  [git](https://git-scm.com/downloads) or [GitHub ](https://desktop.github.com/download/)


- Only for windows: Have you [Rtools](https://cran.r-project.org/bin/windows/Rtools/) already? `devtools::find_rtools()` must be TRUE

- Python as the engine for the main scripts

- Miniconda (conda) environment as the package containing all Python dependencies

- R shiny is for the dashboard decision support system.


Some of the different computers might have particular conditions or requirements for installing all these components, so we made a section where you can find known issues an their potential solution. If you found a new one please share it with us so potential new users can see it.



  ***If you have any issues installing and setting up the software, please check the [KNOWN-ISSUES PAGE](https://github.com/connectingLandscapes/cola/blob/main/inst/docs/md_known_issues.md) (open in a new tab if launched in R/Rstudio/Shiny) for potential solutions*** 
   

------------------
------------------

#### **1.  Install cola R package.**

Consider use the **option 3 (None)** for installing new packages at the first try. If an error arises, update **all of them (option 1)**.

```
if (!require(devtools)){
   install.packages('devtools')
}
devtools::install_github('connectingLandscapes/cola') ## option 3: None
```

The installation log will shown in console (type 3 for fast installation):

```
Downloading GitHub repo connectingLandscapes/cola@HEAD
These packages have more recent versions available.
It is recommended to update all of them.
Which would you like to update?

 1: All                                    
 2: CRAN packages only                     
 3: None                              

 ... (more text)

-- R CMD build -----------------------------------------------------------------------------------------------------------
v  checking for file 'C:\tempR\RtmpM1f5us\remotes11c8f5a7145\connectingLandscapes-cola-5abb74a/DESCRIPTION' (378ms)
-  preparing 'cola': (1.8s)
v  checking DESCRIPTION meta-information ...
-  checking for LF line-endings in source and make files and shell scripts
-  checking for empty or unneeded directories
     NB: this package now depends on R (>= 3.5.0)
     WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects: 'cola/data/crs_df.rda'  'cola/data/perfomance-uper_per.rda'  WARNING: Added dependency on R >= 3.5.0 because serialized objects in  serialize/load version 3 cannot be read in older versions of R.  File(s) containing such objects:  'cola/./showcase/showcase.RData'
-  building 'cola_0.1.0.tar.gz'
   
Installing package into ‘C:/Users/Admin/Documents/R/win-library/4.0’
(as ‘lib’ is unspecified)
* installing *source* package 'cola' ...
** using staged installation
** R
** data
*** moving datasets to lazyload DB
** inst
** byte-compile and prepare package for lazy loading
** help
No man pages found in package  'cola' 
*** installing help indices
** building package indices
** testing if installed package can be loaded from temporary location
*** arch - i386
*** arch - x64
** testing if installed package can be loaded from final location
*** arch - i386
*** arch - x64
** testing if installed package keeps a record of temporary installation path
* DONE (cola)
```
  
  Here, the cola R package is installed. Now we need to set up A) the `conda` environment and B) the shiny the dashboard


------------------


#####  **2. Setting up cola requirements:**

This might take several minutes and will show the following message once you *get all the steps installed*. But during the first time you will have a longer log. See next section.

```
library(cola)
cola::setup_cola()
```

A final success installation should print:
```
  +Step 1/5: Installing & checking reticulate R package
    `reticulate` installed already!
  +Step 2/5 Installing & checking miniconda
    miniconda found at C:/Users/Admin/AppData/Local/r-miniconda!
  +Step 3/5 Installing & checking conda environment
    `cola` conda environment installed in C:\Users\Admin\AppData\Local\r-miniconda\envs\cola/python.exe
    `cola` conda environment named correctly!
    The python version is Python 3.9.19
  +Step 4/5 Installing & checking conda modules
    All required conda modules installed!
  +Step 5/5 Setting up local variables
    === Ready to connect landscapes! ===
```

----------------------

**If there's the first time installing cola you will see the following messages during the different steps:**


```
   +Step 1/5: Installing & checking reticulate R package
```


If `reticulate` R package doesn't exists, you need to install it. A prompt asking if you want to install packages that requires compilation will appear. If you say "yes" this might take longer. You can also say "No", but it might require update the packages at some point, depending in your current.

![Fig. 1a. Compilation message](./01_question_compilation.JPG "Install reticulate")

You might need to restart your R session as well. Say "YES".

![Fig. 1b. Compilation message](./02_restart_session.JPG "Restart R session")


Also, you might be asked to install a reticulate conda environment. Say "YES".

![Fig. 1c. Compilation message](./03_instal_conda_reticulate.png "Install conda r-reticulate")

A correct `reticulate` R package installation will show:

```
  +Step 1/5: Installing & checking reticulate R package
    `reticulate` installed already!
```

----------------------


```
  +Step 2/5 Installing & checking miniconda
```


If miniconda is not installed, you will be asked to accept the installation. Say "yes"


![Fig. 2. Miniconda installation](./miniconda-02-execute.png "Install miniconda")


A correct miniconda installation will show:

```
  +Step 2/5 Installing & checking miniconda
    miniconda found at C:/Users/Admin/AppData/Local/r-miniconda!
```

----------------------

```

+Step 3/5 Installing & checking conda environment

+Step 4/5 Installing & checking conda modules

```


This steps creates a conda environment, named `cola` by default. You will ask to create the environment. Say "yes"

![Fig. 2. Install cola conda environment](./question_cola.JPG "Install miniconda")



It will try several ways in the case some options doesn't works in your machine:

    1. Install the environment using a YML file. If fails,
    
    2. Creating the environment with the required libraries. If fails,
  
    3. Create the environment, and then installing each library independently
  
  
  
This step will take some minutes since requires download spatial and specialized packages.


A correct installation should shown a message like:

```
  +Step 3/5 Installing & checking conda environment
  
    `cola` conda environment installed in C:\Users\Admin\AppData\Local\r-miniconda\envs\cola/python.exe
    
    `cola` conda environment named correctly!
    
    The python version is Python 3.9.19
    
  +Step 4/5 Installing & checking conda modules
  
    All required conda modules installed!
    
```

----------------------

```
  +Step 5/5 Setting up local variables
```

During this step R will define the ways to interact with the python version and environment. The installation script will find the path to the scripts and identify how your system should call cola python. At the same time, here we load all the required modules to test if were installed correctly, and run a sample analysis to ensure 

If all the test passed, this will be the shown:

```
  +Step 5/5 Setting up local variables
    === Ready to connect landscapes! ===
```



#####  **3. Setting up cola dashboard:**

Now let's install all the native R libraries to deply the dashboard.

```
cola::setup_cola_dss()
```

If all libraries were installed, the following 

```
    === All libraries required for COLA's DSS installed ===
```


#####  **4. Load cola and run some functions:**

Let's run one of the functions. Habitat suitability to habitat resistance: s2res

```
library(cola)
input_tif <- system.file(package = 'cola', 'sampledata/sampleTif.tif')
out_shp <- 'outputPoints.shp'

simulated_points <- points_py(intif = input_tif, outshp = out_shp,
   param3 = 0, param4 =  100, param5 = 50)
  
if(file.exists(out_tif)){
  plot(sf::read_sf(simulated_points$file))  
}
```


#####  **5. Uninstall cola **

```
# Not run
# reticulate::conda_remove('cola')
# remove.packages( "cola" )


```



### .    B) Customize `cola` parameters

There's some parameters we can customize for `cola` performance, in each of the following parameters:

 - Temporal directory `COLA_DATA_PATH`: where a folder will be created for each DSS session, and will be used to save the results. By default we use R temporal folder `tempdir()`. Keep the last slash or backshlash.
 - Number of cores `COLA_NCORES`: To be used for running the programs. Default is `1`.
 - Maximum size of the layer `COLA_DSS_UPL_MB`: Size to load in the shiniy dashboard in mega bytes (MB) units. Default is `250`.
 - Maximum resolution for the layer representation `COLA_VIZ_THREs_PIX`, `COLA_VIZ_RES_NCOL`,
  `COLA_VIZ_RES_NROW`: The original layer created or uploaded will be resampled in the front end leaflet visor. This value is the upper threshold fo pixel numbers in the layer. Layers with number of pixels above this number will be resampled to `c(1000, 1000)` pixels. The default value is `1000000` for the threshold, and `1000` for the resampling size in columns and rows. Increase this number to allow your machine load bigger original files in the front end. Decrease this number to avoid loading big files in your front end. 
 
 **How to edit those variables? Don't edit this variables unless you are an advanced user**
 
 Find your R configuration file `.Renviron` file located in `Sys.getenv("R_USER")`.
 Edit with a text editor this file and edit with caution the values. Keep the proper file separator (\ or /), quotes, decimal separators, etc.

Open the file;
```
file.edit(file.path(Sys.getenv("HOME"), ".Renviron"))

```

Populate the values as you consider:
 
 ```
COLA_PYTHON_PATH="C:\Users\gonza\AppData\Local\r-miniconda\envs\cola/python.exe"
COLA_SCRIPTS_PATH="C:/Users/gonza/AppData/Local/R/win-library/4.2/cola/python"
COLA_DATA_PATH="C:/temp/cola"
COLA_NCORES=1
COLA_DSS_UPL_MB=250
COLA_VIZ_THRES_PIX=1000000
COLA_VIZ_RES_NCOL=1000
COLA_VIZ_RES_NROW=1000
 ```
 
 You should also see `COLA_PYTHON_PATH` and `COLA_SCRIPTS_PATH`, but ***DON'T*** edit them manually.
 
 This file is reset to the default values every time you run `cola::setup_cola()`.
 
