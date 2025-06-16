# Connectiviy analysis software: Connecting landscapes - `cola` 

We performed some tests comparing the resulting layers and the required time for processing them.

You can download the files from this [link](https://drive.google.com/open?id=1Y5dPEHrp4hdEtGLdKuUeoc5qFAgehacm&usp=drive_fs)

Download the following files:

- Resistance layer in format TIF, ASCII or RSG
- Source points in format SHP, TXT, or XY


For executing the different software we use the following scripts and configuration files. X refers to the size of the raster, Y to the number of points:
- UNICOR (Python): execute the file 'unicorX_Y.rip' and adapt the file paths.
- CircuitScape (Julia): execute the file 'juliaSizeX_Ypts.ini' (config file) and 'scriptJuliaX_Y.jl'. 
- CoLa (R or Python): execute the file 'cola_py_cmd.txt' and adapt the file paths
