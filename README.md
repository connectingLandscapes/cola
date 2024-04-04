#  Connectiviy analysis - cola

This package integrates R and Python modules. 



## Installing the package
It's required to install several components (once). The srtucture of this software is:
- R as the base
- Python as the engine
- Miniconda (conda) environment as the package containing all python dependencies
- R shiny for the dashboard Decision support system.

Some of the different computers might have particular conditions or requirements for installing all these components, so we made a section where you can find known issues an their potential solution. If you found a new one please share it with us so potential new users can see it.

1.  Install cola
```{r}
if (!require(devtools)){
   install.packages('devtools')
}

devtools::install_github('connectingLandscapes/cola') ## option 3: None
```
  
2. Setting up cola requirements:
```{r}
library(cola)
cola::setup_cola()

## equivalent to 
# cola::setup_cola(envName = 'cola', libs2Install = c('gdal', 'h5py', 'numexpr', 'rasterio', 'pytables', 'pandas', 'cython', 'numba', 'networkit', 'fiona', 'shapely', 'geopandas', 'scikit-image'), nSteps = 5)
```

3. Setting up cola dashboard:
```{r}
cola::setup_cola_dss()
```

4. Run function:
```{r}
cola::setup_cola_dss()
```

5. Run DSS:
```{r}
cola::cola_dss()
```

6. Uninstall cola:
```{r}
# remove.packages( "cola" )
```




Tests:
- Windows
- Linux

