======
README
======

----------------- 
CDPOP 1.3 release
-----------------
  
Welcome to the CDPOP v1.3 release! CDPOP is a program to simulate mutation, gene flow, genetic drift, and selection in complex landscapes for a wide range of biological and evolutionary scenarios. This release includes installation instructions, version notes, some examples, and technical documentation. 
  
Program Contributors: Erin Landguth, Brian Hand, Joe Glassy, Mike Jacobi, Tyler Julian, Allen Warren, Brenna Forester, Sam Cushman, Andrew Eckert, Andrew Shirk, Amy Whipple, Mitra Mennon
Link: http://github.com/ComputationalEcologyLab/CDPOP/
Version: 1.3.xx 
Python: 3.x
Release Date: 2022.03.09
README Update: 2022.03.09 (ell)
  
--------
Contents
--------
  
Included in this release are the following:

src -> CDPOP source files
doc -> README.txt, user manual, and history
data -> test example files
  
---------------------------------------
Requirements and Pre-requisite Software
---------------------------------------

Baseline Requirements. CDPOP requires the Python3 interpreter, NumPy package, and SciPy package. 

-----------------------
CDPOP Installation
----------------------- 

Linux or Windows: Unpack the CDPOP Archive. Navigate to the directory on your PC where you wish to install CDPOP, and unpack the supplied zip archive file using a free archive tool like 7Zip (7z.exe), Pkunzip, Unzip, or an equivalent. Seven-Zip (7Z.exe) is highly recommended since it can handle all common formats on Windows, MAC OS X and Linux. On Windows, it is best to setup a project specific modeling subdirectory to perform your modeling outside of any folder that has spaces in its name (like "My Documents").

-----------------
Example CDPOP Run
-----------------

The example run is for 16-points representing individuals with a Euclidean distance cost distance matrix. To run the following example, follow these steps:

1. Double check that the 3 directories provided in the archive are in the same directory. 

2. The included file inputvars.csv in the data directory specifies the parameters that can be changed and used in a sample CDPOP run. Open inputvars.csv in your editor of choice. A spreadsheet program like Microsoft Excel, allows for easy editing of the tabular values.

3. There will be 3 lines of information in inputvars.csv: a header line and 2 lines of information corresponding to 2 separate CDPOP runs (batch process). See the user_manual that contains a breakdown for each column header and the parameters that can be changed. The Input listed is for the first row in the file. Make sure you save inputvars in the same format – a comma delimited file. Select ‘Yes’ or ‘OK’ for any Excel questions about saving in this format.

5. Start the program: For example, if you use python from the command line, then open a terminal window and change your shell directory to the CDPOP src home directory (i.e., > cd C:\"homedirectorylocation"\src). 

6. Run the program: There are a number of ways to run this program. If you are using a command shell you can run the program by typing “python CDPOP.py C:/"homedirectorylocation"/data inputvars.csv output_test”. Note that there are 5 arguments here: "python" starts python, "CDPOP.py" runs CDPOP program, "C:/"homedirectorylocation"/data" is the directory location of the input test files, "inputvars.csv" is the parameter file, and "output_test" is the name of the directory that will be created with output in the directory specified by the second third argument.

7. Check for successful model run completion: The program will provide step-by-step output in the Shell window. Once completed, a simulation time will be printed out and folders batchrun0mcrun0, batchrun0mcrun1, batchrun0mcrun2, batchrun0mcrun3, batchrun0mcrun4, and batchrun1mcrun0 will be created in your CDPOP home directory to store output from the separate batch and/or Monte-Carlo runs. These folders are located in the data folder specified in (6). The output folder will have a unique date/time stamp after the name in case you want to run multiple CDPOP runs in this same directory. The program will also provide a log file with program steps in your CDPOP home directory.

Happy Simulations!

Erin.

Contact Information
Erin Landguth
Computational Ecology Laboratory
The University of Montana
32 Campus Drive
Missoula MT, 59812-1002
Erin.landguth@mso.umt.edu
