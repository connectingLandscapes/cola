# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

2+2

import pandas as pd
import numpy as np

def fpd(expr, dic):
    """
    Esta función calcula una expresión matemática
    :param expr: La función parte derecha
    :param dic: Diccionario con valores a incluir en la ecuación y con nombres correspondientes en el argumento 'expr'
    :return: Retorna el valor de la expresión 
    """
    for key in dic:
        exec(key + ' = ' + str(dic[key]))
        
    #exec('z = '+ expr)
    z = eval(expr)
    return z
#
#
def exp2(z):
    return z**2
#
dic =  {'a': 1,
         'b': 2.0,
         'c': 3,
         'd': 4,
         'e': 5,
         'li': [1.0, 2, 3, 4, 4, 5, 6]}
         
       
fpd('b/sum(li)', dic)
fpd(expr = 'exp2(b)/sum(li)', dic)


# Set input and output files
root = 'C:/Dropbox/Msc/02_ModHidro/Mod_PH/Metricas/'
dataFile = root + 'Obs-Sim.xlsx'
df = pd.read_excel(dataFile, sheet_name= 'data', index_col = 'Fecha')
obs = df['qobs']

dic2 = {'obs': np.ndarray.tolist(obs.values)}
fpd('np.mean(obs)', dic2)

#################
import ee
ee.Authenticate( )
ee.Initialize(project = eeproject)


print('EE initialized')

tasks = ee.batch.Task.list()
for t in tasks:
    print(t.status())
ee.data.listOperations()

##

predictors = predictors.clip(aoi)
   
# =============================================================================
# Generate 5,000 random points
data_cor = predictors.sample(scale=grainSize, numPixels=5000, geometries=True)
# Extract predictor variable values
pvals = predictors.sampleRegions(collection=data_cor, scale=grainSize)
# =============================================================================
