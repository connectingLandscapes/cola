# `cola` installation test


Consider report any bug via issue on github, or using this [slides](https://docs.google.com/presentation/d/1VXRMBr9OgHOxRzSsFQxm0nx2jlxwIWMg9Vtj_O4yUU0/edit#slide=id.g2c93da8899f_0_8), [document](https://docs.google.com/document/d/11olc9-V-T748g6WstRyv35wZtgGZnBvYZ0F85AyUF-k/edit), or [form](https://docs.google.com/forms/d/e/1FAIpQLSdFsM1e02biuauaWE4Svwtu5QMneKU7Ilfa8pAJHiRy3a-KGw/viewform)


Here is the repo with:

-Intro: https://github.com/connectingLandscapes/cola

-Installing: https://github.com/connectingLandscapes/cola/blob/main/inst/docs/md_cola_install.md

-Installing known-issues: https://github.com/connectingLandscapes/cola/blob/main/inst/docs/md_known_issues.md

- Functions documentation: https://github.com/connectingLandscapes/cola/blob/main/inst/docs/md_colafun.md 


### Run the following lines

```
# start a fresh R session
.rs.restartR() 

# remove old versions
remove.packages('cola') 

# start a fresh R session
.rs.restartR()

# Get the package. Set option 3: None
devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never') 

# load the library
library(cola) 

# Run errors check. Should arise some since the package is not configured.
cola::diagnose_cola() 

# set up everything, and install miniconda. Takes severla minutes -----
cola::setup_cola() 

# Hopefully all good
cola::diagnose_cola() 

## Check internal paths. Must not be empty
Sys.getenv(c('COLA_PYTHON_PATH', 'COLA_SCRIPTS_PATH')) 

# Install dashboard libraries. Takes severla minutes -----
cola::setup_cola_dss() 

# run the dashboard
cola::cola_dss() 

```
