library(rgee)
library(stars)
library(sf)
#install.packages('googleCloudStorageR')


library(reticulate)

miniconda_path()
py_config()

# 'C:/Users/gonza/AppData/Local/r-miniconda'
# python.exe
# r-reticulate


minic_path <- miniconda_path()
setwd(minic_path)
  
if ( dir.exists(minic_path) ) {
  
  mini_act <- system('conda activate base', intern = TRUE) #
  mini_ini <- system('conda init cmd.exe', intern = TRUE) #
  mini_act2 <- system('conda activate base', intern = TRUE) #
  
  system('Scripts/conda-script.py activate base', intern = TRUE)
  system('Scripts/activate.bat activate base', intern = TRUE)
  
  (mini_init3 <- system('Scripts/activate.bat base', intern = TRUE))
  
  (mini_init3 <- system(paste0('Scripts/conda.exe create -n cola networkit networkx ',
                               'h5py pip ipykernel ipywidgets pandas geopandas ',
                               '-c conda-forge --yes'), intern = TRUE))
  
  # https://rdrr.io/github/r-spatial/rgee/src/R/ee_install.R#sym-ee_install_set_pyenv_env
  ### TASK --- run the file || conda_01_install_RUN-AS_ADMIN.bat | Proceed ([y]/n) ? 
  # call %USERPROFILE%\Anaconda3\Scripts\activate.bat %USERPROFILE%\Anaconda3
  # ::call C:\ProgramData\Miniconda3\Scripts\activate.bat C:\ProgramData\Miniconda3
  # call conda activate base
  # call conda create -n cdujlab networkit networkx h5py pip ipykernel ipywidgets pandas geopandas jupyterlab nb_conda_kernels jupyter_contrib_nbextensions ipysheet -c conda-forge
  # call conda activate cdujlab
  # call conda info --envs
  # 
  
  
  
  
  if (file.exists(file.path(minic_path, 'python.exe')))
  
  
  
}


ins <- tryCatch(install_miniconda(path = miniconda_path(), update = FALSE, force = FALSE), error = function (e) e)


if (class('error') %in% ins & grep('Miniconda is already installed')){
  
} 