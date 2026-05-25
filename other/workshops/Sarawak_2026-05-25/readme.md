# Workshop
## Kasane, Botswana. April 1st-3rd 2025

##### *Scan this code*
![QRcode](https://raw.githubusercontent.com/connectingLandscapes/cola/refs/heads/main/workshops/KAZA_2025-04-01/qr-kaza.png)


# Agenda

1. Intro by David McDonald. 
2. Concepts and workflow by Dawn Burnham. 
3. CoLa DSS by Patrick Jantz and Ivan Gonzalez. [Slides](https://docs.google.com/presentation/d/1ZSHY7Xd6p3OlisT7DTARHNSgWUdjwYd4/edit?slide=id.g3e9ee1163d6_0_69#slide=id.g3e9ee1163d6_0_69)
4. Study case by Andrew Hearn.


# Practical session

 1. Access servers here, wait 10 seconds, open in anew tab, and copy session ID:
    - [Server A]( http://34.55.156.96:3838/connecting-landscapes/)
    - [Server B](http://34.45.121.125:3838/connecting-landscapes/)  
    
 3. Download some sample data
    **These are sample data. They are incomplete and only to be used for training purposes. Please do not share the datasets beyond the people in this workshop.**
    - [SDM 1km](https://rcdata.nau.edu/above_gedi/sarawak/sdm1km/)
    - [ZIP file of 39 species, 85MB](https://rcdata.nau.edu/above_gedi/sarawak/sdm1km/__sdm1km.zip)
    - [SDM 90m](https://rcdata.nau.edu/above_gedi/sarawak/proj/)
    
 4. Fill out this survey to tell us about yourselves and your experience (2 min):
 (https://forms.gle/jW98Aa5X7Gp5FnnU6). And an another survey about your the DSS perspective (https://docs.google.com/forms/d/e/1FAIpQLSdE9QMHnBv3FhXy8zJgOrvvx39ltetJf7-mtIbv4kZQuElubg/viewform?usp=sf_link)


# Installing
Follow the steps presented [here](https://github.com/connectingLandscapes/cola/blob/main/inst/docs/md_cola_install.md)
 
Main steps (follow this order):
1. Install Git
2. Install R
3. Install Rtools
4. Install Rstudio (optional)

Then, in R run (use your house WIFI):
```
## Install package manager
if(!require(devtools)){install.packages('devtools')}

## Install or update CoLa. Take 2 mins
devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never') # Installing CoLa R package

## Set up CoLa. Takes 30 mins the first time you run it
## Run this function until all problems are solved
cola::setup_cola(ask = FALSE, dss = TRUE) 

# Run the DSS
cola::cola_dss()
```
