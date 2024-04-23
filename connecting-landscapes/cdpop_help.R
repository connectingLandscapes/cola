setwd('/home/shiny/connecting-landscapes/performance-tests/cdpop/')
#list.files()
#dir.create('/home/shiny/connecting-landscapes/lib/CDPOP/data/fromCola/')


inputvars_intros <- read.csv("0_inputvars.csv")
R_vars <- inputvars_intros[1, ]
t(R_vars)

R_vars$xyfilename <- 'xy_from_cola' # with no .csv 
R_vars$agefilename <- 'AgevarsA.csv' # with .csv ending
R_vars$mcruns <- 10
R_vars$looptime <- 100
R_vars$output_years <- 100
R_vars$gridformat <- 'cdpop'
R_vars$output_unicor <- 'N'
R_vars$cdclimgentime <- 0
R_vars$matecdmat <- ''
R_vars$dispcdmat <- ''


xyfile <- read.csv("~/connecting-landscapes/lib/CDPOP/data/xyfiles/xyED16.csv")
xy2 <- 0


# python CDPOP.py %userprofile%\dockerdata\CDPOP\data inputvars.csv outAnac1
## Works
# system(
#   paste0('/home/shiny/anaconda3/bin/python /home/shiny/connecting-landscapes/lib/CDPOP/src/CDPOP.py ',
#          '/home/shiny/connecting-landscapes/lib/CDPOP/data ',
#          'inputvars.csv outTry1'))



setwd('/home/shiny/connecting-landscapes/performance-tests/cdpop/cola')
list.files()
vars <- read.csv("inputvars.csv")
vars$

system(
  paste0('/home/shiny/anaconda3/bin/python /home/shiny/connecting-landscapes/lib/CDPOP/src/CDPOP.py ',
         '/home/shiny/connecting-landscapes/performance-tests/cdpop/cola ',
         'inputvars.csv', ' ',
         'outCDPOP'))
