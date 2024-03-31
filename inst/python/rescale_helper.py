# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 09:38:29 2023
Script to graph suitaiblity against resistance
using different parameter sets.
@author: pj276
"""

#%%
# Imports
import pandas as pd
import sys
import cola_functions as cf
import numpy as np

def main() -> None:
    #%%
    # INPUTS
    # Suitability minimum value
    sMin = sys.argv[1]
    
    # Suitability maximum value
    sMax = sys.argv[2]

    # Maximum resistance value (min is set to 1)
    rMax = sys.argv[3]
    
    # Shape parameter
    shpParam = sys.argv[4]
    
    # No data value of suitability raster
    ndVal = sys.argv[5]
    
    # Convert suitability min to float or integer
    try:
        if cf.is_float(sMin):
            sMin = float(sMin)
        elif sMin.isdigit():
            sMin = int(sMin)
        else:
            float(sMin)
            int(sMin)
    except ValueError:
        print('Suitability min must be either a float or integer')

    # Convert suitability max to float or integer
    try:
        if cf.is_float(sMax):
            sMax = float(sMax)
        elif sMax.isdigit():
            sMax = int(sMax)
        else:
            float(sMax)
            int(sMax)
    except ValueError:
        print('Suitability max must be either a float or integer')

    # Convert resistance max to float or integer
    try:
        if cf.is_float(rMax):
            rMax = float(rMax)
        elif rMax.isdigit():
            rMax = int(rMax)
        else:
            float(rMax)
            int(rMax)
    except ValueError:
        print('Resistance max must be either a float or integer')

    # Convert shape param to float or integer
    try:
        if cf.is_float(shpParam):
            shpParam = float(shpParam)
        elif shpParam.lstrip('-').isdigit():
            shpParam = int(shpParam)
        else:
            float(shpParam)
            int(shpParam)
    except ValueError:
        print('Shape param must be either a float or integer')

    # Convert no data value to float or integer
    try:
        if cf.is_float(ndVal):
            ndVal = float(ndVal)
        elif ndVal.lstrip('-').isdigit():
            ndVal = int(ndVal)
        elif ndVal == 'nan':
            ndVal = -9999
        else:
            float(ndVal)
            int(ndVal)
    except ValueError:
        print('No data value must be either a float, integer, or nan string')

    #%%
    # Habitat suitability transform plot
    exampleValues = np.linspace(sMin,sMax,10)
    # Rescale to 0-1 range, then to 1 - max resistance range
    resVec = cf.vecrescale(exampleValues, sMin, sMax, rMax, ndVal, -9999, True, shpParam)
    dd = pd.DataFrame(np.array([exampleValues,resVec]).transpose(),columns=['suitability','resistance'])
    ax = dd.plot.line(x='suitability',y='resistance',lw=3, fontsize=16,grid=True)
    #ax.set(xlabel='Suitability',ylabel='Resistance')
    ax.set_xlabel('Suitability', fontsize=16)
    ax.set_ylabel('Resistance', fontsize=16)
    ax.get_legend().remove()

if __name__ == "__main__":
    main()