# Connectiviy analysis `cola` 


### Installation 

This package integrates R and Python modules. 

It's required to install several components (once). The structure of this software is:
- R as the base. Use the latest available [here](https://cran.r-project.org/bin/windows/base/)
- Have you [Rtools](https://cran.r-project.org/bin/windows/Rtools/) already? `devtools::find_rtools()` must be TRUE
- Python as the engine for the main scripts
- Miniconda (conda) environment as the package containing all python dependencies
- R shiny for the dashboard Decision support system.

Some of the different computers might have particular conditions or requirements for installing all these components, so we made a section where you can find known issues an their potential solution. If you found a new one please share it with us so potential new users can see it.


#### **1.  Install cola R package.**

Consider use the option 3 (None) for installing new packages at the first try. If an error arises, update all of them (option 1)
```{r}
if (!require(devtools)){
   install.packages('devtools')
}

devtools::install_github('connectingLandscapes/cola') ## option 3: None
```

The installation log Will shown in console:

```
Downloading GitHub repo connectingLandscapes/cola@HEAD
These packages have more recent versions available.
It is recommended to update all of them.
Which would you like to update?

 1: All                                    
 2: CRAN packages only                     
 3: None                              

...

** testing if installed package can be loaded from final location
*** arch - i386
*** arch - x64
** testing if installed package keeps a record of temporary installation path
* DONE (cola)
```
  
  
2. Setting up cola requirements:
This might take several minutes and 
```{r}
library(cola)
cola::setup_cola()

## equivalent to 
# cola::setup_cola(envName = 'cola', libs2Install = c('gdal', 'h5py', 'numexpr', 'rasterio', 'pytables', 'pandas', 'cython', 'numba', 'networkit', 'fiona', 'shapely', 'geopandas', 'scikit-image'), nSteps = 5)


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

3. Setting up cola dashboard:
```{r}
cola::setup_cola_dss()

    === All libraries required for COLA's DSS installed ===
>
```

4. Run function:
```{r}
cola::setup_cola_dss()

```

5. Run DSS:
```{r}
cola::cola_dss()

    === All libraries required for COLA's DSS installed ===
>
```

6. Uninstall cola:
```{r}
# remove.packages( "cola" )
```





Tests:
- Windows
- Linux

