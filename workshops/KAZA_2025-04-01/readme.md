# Workshop
## Kasane, Botswana. April 1st-3rd 2025

##### *Scan this code*
![QRcode](https://raw.githubusercontent.com/connectingLandscapes/cola/refs/heads/main/workshops/KAZA_2025-04-01/qr-kaza.png)

![banner](https://raw.githubusercontent.com/connectingLandscapes/cola/refs/heads/main/workshops/KAZA_2025-04-01/ebanner.png)


# Agenda

1. Intro by David McDonald. [Slides](https://drive.google.com/open?id=1oAsKZAhw81zyzPKbyS0ehC6fG1ZYir1y&usp=drive_fs)
2. Concepts and workflow by Dawn Burnham. [Slides](https://docs.google.com/presentation/d/1nvz7o9J4YcKl4p9okQ_1BGo8AiVj5awr?rtpof=true&usp=drive_fs)
3. CoLa DSS. [Slides](https://drive.google.com/open?id=1o5YLn2k49cGJCzKHz5ncWYrKj5D4TUuZ&usp=drive_fs)
4. CoLa tutorial by Ivan Gonzalez. [Slides](https://docs.google.com/presentation/d/18iNtXGxe_NAlaNdxGC9xb_OBJrwRIzXI/edit?usp=sharing&ouid=103068293807996405041&rtpof=true&sd=true)


# Practical session

 1. Get to this link:
    - [Server A](http://34.57.99.166:3838/connecting-landscapes/). Wait for 10 seconds to load. Open on a new tab.
    
    - [Server B](http://34.57.191.163:3838/connecting-landscapes/)  
    
 3. Download some sample data> (https://github.com/connectingLandscapes/cola/tree/main/workshops/KAZA_2025-04-01/data)
    **These are sample data. They are incomplete and only to be used for training purposes. Please do not share the datasets beyond the people in this workshop.**
 4. Fill out this survey to tell us about yourselves and your experience (2 min):
 (https://forms.gle/jW98Aa5X7Gp5FnnU6). And an another survey about your the DSS perspective (https://docs.google.com/forms/d/e/1FAIpQLSdE9QMHnBv3FhXy8zJgOrvvx39ltetJf7-mtIbv4kZQuElubg/viewform?usp=sf_link)
 5. Follow the next steps: [Tutorial or slides here](https://docs.google.com/presentation/d/18iNtXGxe_NAlaNdxGC9xb_OBJrwRIzXI/edit?usp=sharing&ouid=103068293807996405041&rtpof=true&sd=true)


# Installing
Follow the steps presented [here](https://github.com/connectingLandscapes/cola/blob/main/inst/docs/md_cola_install.md)
 
Main steps (follow this order):
1. Install Git
2. Install R
3. Install Rtools
4. Install Rstudio (optional)

Then, in R run (use your room's WIFI):
```
if(!require(devtools)){install.packages('devtools')}
devtools::install_github('connectingLandscapes/cola', dependencies = NA, upgrade = 'never') # Installing CoLa R package
cola::setup_cola(ask = FALSE, dss = TRUE) # Setup all CoLa components. Run this line until all problems are solved.
```
