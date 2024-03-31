# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 19:40:33 2023

@author: pj276
"""
# Start by array slicing the hdf file by the number of
# nodes that will fit in memory. For each slice,
# find the pairwise combinations from reOrder
# that can be calculated with the set of nodes
# from that slice.
# After iterating through all slices, there will
# still be many pairwise combinations left. 
# Find the nodes that are part of these pairwise
# combinations, break them into batches that will
# fit in memory, then iterate through them.
# Keep doing this until all combinations have been 
# calculated.

# if group size is > batch limit, break group into
# batches


# Sort reOrder
atest = reOrder
atest = atest[np.lexsort((atest[:,1],atest[:,0]))]

# pandas
import pandas as pd
grp = pd.DataFrame(atest)
grp = grp.rename(columns={0: "s1", 1: "s2"})
grp1 = grp.groupby('s1').filter(lambda x: len(x) >= 300).groupby('s1')

grp2 = grp.groupby('s1').filter(lambda x: len(x) < 300).groupby('s1')

size = grp2.size()
group_names = []
for name,group in grp2: group_names.append(name)
size = pd.DataFrame({'size': size, 'gname': group_names})

aa = batchesCalc4(grp2, 300)

gbatch = []
for gid in aa[100]:
    gbatch.append(grp2.get_group(gid))
gbatch = pd.concat(gbatch)
gbatch = gbatch.to_numpy()

len(np.unique(gbatch))

# Merge groups up to the batch limit
# e.g. if 300 is the batch limit
# merge groups until the number of nodes
# in them reaches the limit
def batchesCalc(gSize, maxL):
    gList = []
    smallgList = []
    i = 0
    for index, row in gSize.iterrows():
        if (len(smallgList) == 0) and row['size'] > maxL:
            gList.append(row['gname'])
            smallgList = 0
            i = 0
        elif (i < maxL) and ((i + row['size']) < maxL):
            smallgList.append(row['gname'])
            i += row['size']
        elif (i < maxL) and (i + row['size'] > maxL):
            gList.append(smallgList)
            smallgList = []
            gList.append([row['gname']])
            i = 0
    return gList

# Merge groups up to the batch limit
# e.g. if 300 is the batch limit
# merge groups until the number of nodes
# in them reaches the limit
def batchesCalc2(grup, maxL):
    gList = []
    smallgList = []
    sid = []
    for name,group in grp2:
        sz = len(np.unique(group.to_numpy()))
        while sz < maxL:
            smallgList.append(name)
            sid = sid + np.unique(group.to_numpy()).tolist()
            sz += len(np.unique(sid))
        sid = []
        gList.append(smallgList)
        smallgList = []
    return gList
        
        if ((len(smallgList) == 0)) and (len(group) > maxL):
            gList.append(name)
            smallgList = 0
            i = 0
        elif (i < maxL) and ((i + sz) < maxL):
            smallgList.append(name)
            sid = sid + np.unique(group.to_numpy())
            sz = len(np.unique(sid))
            i += sz
        elif (i < maxL) and ((i + sz) > maxL):
            gList.append(smallgList)
            smallgList = []
            gList.append([name])
            i = 0
    return gList 
                
    for index, row in gSize.iterrows():
        if (len(smallgList) == 0) and row['size'] > maxL:
            gList.append(row['gname'])
            smallgList = 0
            i = 0
        elif (i < maxL) and ((i + row['size']) < maxL):
            smallgList.append(row['gname'])
            i += row['size']
        elif (i < maxL) and (i + row['size'] > maxL):
            gList.append(smallgList)
            smallgList = []
            gList.append([row['gname']])
            i = 0
    #gList.append(otra)
    return gList


# in them reaches the limit
def batchesCalc3(grup, maxL):
    gList = []
    smallgList = []
    uniqueNodes = []
    sz = 0
    for name,group in grp2:
        #sz = len(np.unique(group.to_numpy()).tolist())
        nuns = uniqueNodes + np.unique(group.to_numpy()).tolist()
        nuns = list(set(nuns))
        szi = len(nuns)
        if (sz < maxL) and ((sz + szi) < maxL):
            smallgList.append(name)
            uniqueNodes = uniqueNodes + nuns
            sz += szi
        elif (sz < maxL) and (sz + szi > maxL):
            gList.append(smallgList)
            smallgList = []
            gList.append([name])
            sz = 0
            nuns= []
        else:
            gList.append(smallgList)
            smallgList = []
            sz = 0
            nuns= []
    return gList

def batchesCalc3(grup, maxL):
    bigList = []
    smallList = []
    uniqueNodes = []
    ncount = 0
    for name,group in grp2:
        # Increment
        uniqueNodes = uniqueNodes + np.unique(group.to_numpy()).tolist()
        nuns = list(set(uniqueNodes)) # unique node ids
        szi = len(uniqueNodes) # length of unique node ids
        # Check against thresholds
        if (ncount < maxL) and ((ncount + szi) < maxL):
            smallList.append(name)
            uniqueNodes = list(set(uniqueNodes + nuns))
            ncount = szi
        elif (ncount < maxL) and ((ncount + szi) > maxL):
            bigList.append(smallList)
            smallList = []
            bigList.append([name])
            ncount = 0
            nuns= []
            uniqueNodes = []
#        else:
#            gList.append(smallgList)
#            smallgList = []
#            sz = 0
#            nuns= []
    return bigList

def batchesCalc4(grup, maxL):
    bigList = []
    uniqueNodes = []
    ncount = 0        
    for name,group in grp2:
        while ncount < maxL:
            uniqueNodes = uniqueNodes + np.unique(group.to_numpy()).tolist()
            nuns = list(set(uniqueNodes)) # unique node ids
            ncount = len(uniqueNodes)
        bigList.append(nuns)
        ncount = 0
        uniqueNodes = []
        nuns = []
    return bigList
        
        # Increment
        uniqueNodes = uniqueNodes + np.unique(group.to_numpy()).tolist()
        nuns = list(set(uniqueNodes)) # unique node ids
        szi = len(uniqueNodes) # length of unique node ids
        # Check against thresholds
        if (ncount < maxL) and ((ncount + szi) < maxL):
            smallList.append(name)
            uniqueNodes = list(set(uniqueNodes + nuns))
            ncount = szi
        elif (ncount < maxL) and ((ncount + szi) > maxL):
            bigList.append(smallList)
            smallList = []
            bigList.append([name])
            ncount = 0
            nuns= []
            uniqueNodes = []
#        else:
#            gList.append(smallgList)
#            smallgList = []
#            sz = 0
#            nuns= []
    return bigList





ul1 = np.unique(grp2.get_group(0).to_numpy()).tolist()

ul1 = list(set(ul1 + np.unique(grp2.get_group(6).to_numpy()).tolist()))

        
        while sz < maxL:
            smallgList.append(name)
            sid = sid + np.unique(group.to_numpy()).tolist()
            sz += len(np.unique(sid))
        sid = []
        gList.append(smallgList)
        smallgList = []
    return gList
            

       
    for j in gSize:
        
    
        i += 1
        
        tb1 = np.array_split(reO, i)
        blengths = [len(np.unique(j)) for j in tb1]
        maxblengths = np.max(blengths)     
    return i

#grp = grp.size().sort_values(ascending=False)
for name, group in grp:
    print(name)

grp.size


from itertools import combinations



nrows, ncols = reOrder.shape
ddtype={'names':['f{}'.format(i) for i in range(ncols)],
       'formats':ncols * [reOrder.dtype]}

reOrderTicker = []
#----
# Move indices over by 30 each time
x = np.arange(len(sources)+1)
y = np.array_split(x, nCBatches)
for s in range(5):
    s = y[s]
    s = s+225
    if len(s[s>len(sources)]) > 0:
        s = np.delete(s, np.where(s>len(sources)))
    #----

    pncx = list(combinations(s, 2))
    pncx = [np.array(i) for i in pncx]
    pncx = np.vstack(pncx)
    nmx = np.intersect1d(reOrder.view(ddtype), pncx.view(ddtype), return_indices=True)
    reOrderTicker = np.unique(reOrderTicker + nmx[1].tolist()).tolist()
    print(len(reOrderTicker)/len(reOrder)*100)

zz = np.delete(reOrder, reOrderTicker, axis=0)

ac1 = np.unique(zz[:,0], return_counts=True)
gg = np.vstack((ac1[0],ac1[1])).T
atest = gg[gg[:,1].argsort()[::-1]]

ac2 = np.unique(zz[:,1], return_counts=True)
gg2 = np.vstack((ac2[0],ac2[1])).T
atest2 = gg2[gg2[:,1].argsort()[::-1]]

s1 = atest[299:450,0]
s2 = atest2[299:450,0]
s = np.hstack((s1,s2)).astype("int")




# Shift indices over by 1/2 batch size
# This becomes like a new sourceBatches but with indices
# shifted to pick up new combinations of sources.
x = np.arange(len(sources)+1)
y = np.array_split(x, nCBatches)
y = [(a+(len(sBatches[0])/2)).astype('int') for a in y]
# Don't think I need to do this because this range would
# have already been included in the first round
# Add indices for 0 to 1/2 batch size
#y.insert(0, np.arange(y[0][0])) 
# Delete indices off the end that are > length sources
yy = y[-1][y[-1] > len(sources)-1]
y[-1] = np.delete(y[-1], np.nonzero(np.in1d(y[-1],yy))[0])

# GEt random combinations of half batch size
import random, copy
x = np.arange(len(sources)+1)
y = np.array_split(x, nCBatches*2)
ry = np.array_split(x, nCBatches*2)
random.shuffle(ry)
aa = [a for a in zip(y, ry)]
s = np.unique(np.hstack((aa[3][0], aa[3][1])))
s = s+121
if len(s[s>len(sources)]) > 0:
    s = np.delete(s, np.where(s>len(sources)))

# Get disjunct indices of sources
x = np.arange(len(sources)+1)
y = np.array_split(x, nCBatches*2)
aa = y[::4] # increment the last number to get every nth element
aa = [a for a in zip(aa[0::2], aa[1::2])] # get adjacent elements
bb = y[1::4]
bb = [b for b in zip(bb[0::2], bb[1::2])]
cc = y[2::4]
cc = [c for c in zip(cc[0::2], cc[1::2])]
dd = y[3::4]
dd = [d for d in zip(dd[0::2], dd[1::2])]

ee = aa + bb + cc + dd
s = np.hstack((ee[2][0], ee[2][1]))

# For example
# aa[0][0] is an array of the first range of numbers
# aa[0][1]] is an array of the second range of numbers
# Putting these two ranges togethe gets a disjunct set of indices

# Putting these two ranges together gets another disjunct set of  indices
# and so on
# aa[1][0] is an array of the first range of numbers
# aa[1][1]] is an array of the second range of numbers

# Same for bb...

# cc concatenates aa and bb


# The reason we do the first array slicing is because
# slicing contiguously is faster than slicing by a list.

ac1 = np.unique(reOrder[:,0], return_counts=True)
gg = np.vstack((ac1[0],ac1[1])).T
atest = gg[gg[:,1].argsort()[::-1]]


ac2 = np.unique(reOrder[:,1], return_counts=True)
gg2 = np.vstack((ac2[0],ac2[1])).T
atest2 = gg2[gg2[:,1].argsort()[::-1]]

s1 = atest[750:900,0]
s2 = atest2[750:900,0]
s = np.hstack((s1,s2)).astype("int")





# Get distance array
tic1 = time.perf_counter()
f = h5py.File(dahdf, mode='r')
darr = f['dmat'][0:500,:]
f.close()
toc1 = time.perf_counter()
print("Calculating batch " + str(b+1) + " took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)


f['dmat'].chunks
#(2, 60000 )

# LZF
afile = r'C:\Users\pj276\Projects\UNICOR_arch\unicor\kaza\t1.h5'
f = h5py.File(afile, mode='a')
f.create_dataset("dmat", (len(sources), len(nodeids)),
                 compression="lzf", shuffle=False,
                 chunks=True, dtype='float32')
f.close()
        
with h5py.File(afile, 'a', ) as f:
    f['dmat'][62:124,:] = ccArr

# Get distance array
tic1 = time.perf_counter()
f = h5py.File(afile, mode='r')
darr = f['dmat'][0:124,:]
f.close()
toc1 = time.perf_counter()
print("Reading took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)

# GZIP
afile = r'C:\Users\pj276\Projects\UNICOR_arch\unicor\kaza\t22.h5'
f = h5py.File(afile, mode='a')
f.create_dataset("dmat", (len(sources), len(nodeids)),
                 compression="gzip", shuffle=False,
                 chunks=True, dtype='float32')
f.close()
        
with h5py.File(afile, 'a', ) as f:
    f['dmat'][0:62,:] = ccArr

# Get distance array
tic1 = time.perf_counter()
f = h5py.File(afile, mode='r')
darr = f['dmat'][0:62,:]
f.close()
toc1 = time.perf_counter()
print("Reading took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)
           
def batchesCalc8(reO, maxL):
    alist = []
    cc = np.array([]).astype('int')
    for index, row in df.iterrows():
        cc = np.hstack((cc,row))
        if len(np.unique(cc)) >= maxL:
            alist.append(np.unique(cc))
            cc = np.array([]).astype('int')
    return alist

grp1 = grp1.apply(lambda x: x)
grp2 = grp2.apply(lambda x: x)

tic1 = time.perf_counter()
atest = batchesCalc8(grp2, 308)
toc1 = time.perf_counter()
print("Calculating batches took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)




# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:06:28 2023

@author: pj276
"""

import tables
dahdf = r"C:/Users/pj276/Projects/UNICOR_arch/unicor/kaza/lcc250_1kpts_1_dahdf.h5"
h5file = tables.open_file(f, 'r')
#table = h5file.root.dmat
#data = table.read([0:10,:])
aa = h5file.root.dmat[8:12, 18:22]
#data = h5file.root.dmat[0:10,:]
h5file.close()


with tb.open_file('your_file.h5', mode='r') as h5file:
    table = h5file.root.entry_name
    # 'entry_name' is name of your HDF5 file
    
    # read into a NumPy array
    data = table.read()

tic1 = time.perf_counter()
f = h5py.File(dahdf, mode='r')
darr = f['dmat'][0:307,:]
f.close()
toc1 = time.perf_counter()
print(" took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)

# New hdf in pytables
import numpy as np
import tables
from contextlib import closing

FILTERS = tables.Filters(complib='zlib', complevel=5)

FILTERS = tables.Filters(complib='blosc2', shuffle=True, complevel=5)



data = np.zeros(10**7)

with closing(tables.open_file('compressed', mode='w', filters=FILTERS)) as hdf:
    hdf.create_carray('/', 'array', obj=data)

with closing(tables.open_file('uncompressed', mode='w')) as hdf:
    hdf.create_array('/', 'array', obj=data)

filters = Filters(complevel=1, complib='blosc', fletcher32=True)
arr = fileh.create_earray(fileh.root, 'earray', atom, (0,2),
                         "A growable array", filters=filters)

# Append several rows in only one call
arr.append(np.array([[1., 2.],
                     [2., 3.],
                     [3., 4.]], dtype=np.float32))

# Print information on that enlargeable array
print("Result Array:")
print(repr(arr))
fileh.close()



fileh = tb.open_file('test5.h5', mode='w')
atom = Float32Atom()
filters = Filters(complevel=1, complib='blosc', fletcher32=True)
arr = fileh.create_earray(fileh.root, 'earray', atom, (0,2),
                         "A growable array", filters=filters, shape = (307, 6006711))

# WORKING
# Create 32 float array in hdf file
afile = r"C:/Users/pj276/Projects/UNICOR_arch/unicor/kaza/testpyt.h5"
fileh = tables.open_file(afile, mode='w')
filters = tables.Filters(complib='blosc2', shuffle=True, complevel=5)
atom = tables.Float32Atom()
a1 = fileh.create_carray(fileh.root, 'dset', atom, shape = (307, 6006711))
fileh.close()

# WRite data
fileh = tables.open_file(afile, mode='w')
fileh[0:10,:] = darr[0:10,:]

table = fileh.root.entry_name





# Append several rows in only one call
arr.append(np.array([[1., 2.],
                     [2., 3.],
                     [3., 4.]], dtype=np.float32))

# Print information on that enlargeable array
print("Result Array:")
print(repr(arr))
fileh.close()

# HDF NOTES
import numpy as np
import tables as tb

afile = r"C:/Users/pj276/Projects/UNICOR_arch/unicor/kaza/testpyt.h5"
fileh = tables.open_file(afile, mode='w')
filters = tables.Filters(complib='blosc2', shuffle=True, complevel=5)
atom = tables.Float32Atom()
a1 = fileh.create_carray(fileh.root, 'dset', atom, shape = (307, 6006711))
fileh.close()


# PYTABLES
fileName = r"C:/Users/pj276/Projects/UNICOR_arch/unicor/kaza/testpyt2.h5"
shape = (1500, 6006711)
atom = tb.Float32Atom()
filters = tb.Filters(complib='blosc2', shuffle=True, complevel=5, fletcher32=False)

h5f = tb.open_file(fileName, 'w')
ca = h5f.create_carray(h5f.root, 'dset', atom, shape,
                       filters=filters)

# Fill a hyperslab in ``ca``.
ca[0:307, :] = darr
h5f.close()

# Re-open a read another hyperslab
inds = np.unique(np.sort(np.random.randint(0,300,307)))

h5f = tb.open_file(fileName)
print(h5f)
#print(h5f.root.dset[0:2, 18:22])
#atest = h5f.root.dset[[0,1,2,3,4,5,6,7,8,9], :]
#atest = h5f.root.dset[0:8, :]
for i in inds:
    atest = h5f.root.dset[i, :]
atest = h5f.root.dset[0:307, :]
h5f.close()
f
# H5PY
afile = r"C:/Users/pj276/Projects/UNICOR_arch/unicor/kaza/testpy5.h5"
f = h5py.File(afile, mode='a')
f.create_dataset("dmat", (1500, 6006711),
                 compression="lzf", shuffle=True,
                 chunks=True, dtype='float32')
f.close()

with h5py.File(afile, 'a', ) as f:
    f['dmat'][0:307,:] = darr

f = h5py.File(afile, mode='r')
for i in inds:
    ppArr = f['dmat'][i, :]
#ppArr = f['dmat'][inds,:]
f.close()


# GROUPING NOTES
def batchesCalc8(reO, maxL):
    alist = []
    cc = np.array([]).astype('int')
    for index, row in reO.iterrows():
        cc = np.hstack((cc,row))
        if len(np.unique(cc)) >= maxL:
            alist.append(np.unique(cc))
            cc = np.array([]).astype('int')
    return alist

grp1 = grp1.apply(lambda x: x)
grp2 = grp2.apply(lambda x: x)

tic1 = time.perf_counter()
atest = batchesCalc8(grp2, 308)

grp1 = grp1.reset_index()
grptest = grp1.groupby('s1').len()

# ORDERING NOES
aa = np.array([[18,19],[18,20],[18,21],[18,22],[18,23],[100,23]])
bb = np.array([[19,20],[19,21],[19,22],[19,23],[19,24],[19,26],[19,27],[100,23],[100,24]])

aau = np.unique(aa)
bbu = np.unique(bb)

# Nodes from 1 also in 2
ix = np.intersect1d(aau,bbu)
# Nodes from 2 not in 1
diff = bbu[~np.isin(bbu,aau)]
# Stack
ixdiff = np.hstack((ix,diff))
# Order
order = ixdiff.argsort()
# Use order to sort the rows of the
# stacked distance array


tic1 = time.perf_counter()
f = h5py.File(dahdf, mode='r')
darr = f['dmat'][0:236,:]
f.close()
toc1 = time.perf_counter()
print("Calculating batches took " +  f"{toc1 - tic1:0.4f} seconds", flush=True)
