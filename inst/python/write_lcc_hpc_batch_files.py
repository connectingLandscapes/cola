import io
import numpy as np

# Change slurm resource requests in the triple quoted text below.
# Make sure the file paths in the srun statement correspond to the files
# on your system. 

# Make sure the filename, f, is relevant for your use case.

# Number of corridors. This can be obtained in the log of the first script
# for running corridors on hpc systems. It can also be obtained by checking
# the number of lines in the reorder.csv file, one of the outputs of the
# 1st script.
maxi = 1498359

# Number of batches you want to run. E.g. if you want to cut procesing time
# in half, run two batches. 
nBatches = 20

# Batch size. Calculate this by dividing the number of corridors by the number
# of batches you want to run.
bsize = int(np.ceil(maxi/nBatches)) #74918

for i in range(1, nBatches + 1):
    # Get end index
    ei = i*bsize
    # Get start index
    si = ei-bsize
    # For last batch, reset end index to number of corridors 
    if ei >= maxi:
        print("last batch")
        ei = maxi
        print((si,ei))
        break
    print((si,ei))
    # Write batch text to file
    # Left justify everything
    aa = f"""#!/bin/bash
#SBATCH --job-name=dhole{i}
#SBATCH --chdir=/scratch/pj276/colatest/output
#SBATCH --output=/scratch/pj276/colatest/logs/job_%A_%a.log
#SBATCH --time=06:00:00
#SBATCH --partition=core
#SBATCH --cpus-per-task=10
#SBATCH --mem 42000

module load anaconda3

srun ~/.conda/envs/colamonsoon/bin/python /home/pj276/projects/cola-main/inst/python/lcc_zarr_cache_test.py /scratch/pj276/usfsip_connectivity/dhole_scens/Accessibility_layer_ssp1_2050_M.tif /scratch/pj276/usfsip_connectivity/dhole_scens/lcc_ssp1_2050_M_r{i}.tif /scratch/pj276/usfsip_connectivity/dhole_scens/dazarr.zarr 0 0 10 None /scratch/pj276/usfsip_connectivity/dhole_scens/reorder.csv /scratch/pj276/usfsip_connectivity/dhole_scens/nodeids.csv {si} {ei}
"""
    f = io.open(f'C:/Users/pj276/Downloads/dhole_scens/lcc2_hpc_job_submit_{i}.sh', 'w', newline='\n')
    f.write(aa)
    f.close()
    

