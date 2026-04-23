x <- read.delim('N:/My Drive/git/cola/inst/ee/sat_ts_fusion/models/models.py', 
                sep = ';', check.names = FALSE, header = FALSE, colClasses = 'character')
x <- x[ , 1]
# x <- y
head(x)
nchar(x)
table(nchar(x) > 80)
table(nchar(x) > 160)
x[nchar(x) > 160]

sum(nc1 > 80) # 94
x[130]
#while ( sum(nc1 > 80) != 0  ){
  #print (1111111111111111111)
for (r in 1:20){
  x1 <- substr(x = x, 0, 1)
  nc1 <- nchar(x)
  y <- NULL
  length(x)
  for(i in 1:length(x)){ # 844 910
    # i = 131;  i = 6
    if (x1[i] == '#' & nc1[i] > 80) {
      print('fix')
      print(i)
      nc <- nchar(x[i])
      xa <- substr(x = x[i], 0, 80)
      xb <- substr(x = x[i], 81, nc)
      y <- c(y, xa, paste0('# ', xb))
    } else {
      y <- c(y, x[i])
    }
  }
  sum(nchar(y) > 80) # 76 67
  length(x) # 844 910
  length(y) # 910 958
  x <- y
  x1 <- substr(x = x, 0, 1)
  nc1 <- nchar(x)
  sum(nc1 > 80) # 67
  
  if( sum(nc1 > 80) != 0){
  }
  print(22222222222222222222)
}

length(x)
table(nchar(x) > 80)
x[nchar(x) > 80]
table(nchar(y) > 80)

writeLines(y, 'N:/My Drive/git/cola/inst/ee/sat_ts_fusion/models/models2.py')
