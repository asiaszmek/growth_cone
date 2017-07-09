import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('No filename')
    f = h5py.File(sys.argv[1])
    grid = f['model']['grid']
    grid_no = f['model']['output']['__main__']['elements']
    species = f['model']['output']['__main__']['species']
    res = f['trial0']['output']['__main__']['population']
    shape = (res.shape[0],res.shape[-1]+1)
    #region name is in -3 of every row grid_no
    regions = f['model']['regions']
    # for g in grid:
    #     print g
    vols = {}
    voxels = {}
    for reg in regions:
        vols[reg] = 0
        voxels[reg] = []
    for reg in regions:
        for i,g in enumerate(grid):
            if regions[g["region"]] in reg:
                vols[reg] += g["volume"]
                voxels[reg].append(i)
            

    header = "time"
    for i,specie in enumerate(species):
        header += " " + specie 
    voli = 0
    tots = []
    for reg in regions:
        if reg == 'default':
            continue
        voli += vols[reg]
        for fi in f:
            if 'trial' in fi:
                file_dend = sys.argv[1]+'_'+fi+'_conc_'+reg
                res = f[fi]['output']['__main__']['population']
                shape = (res.shape[0],res.shape[-1]+1)
                out = np.zeros(shape)
                out[:,0] = f[fi]['output']['__main__']['times']
                for i,specie in enumerate(species):
                    out[:,i+1] = res[:,voxels[reg],i].sum(axis=1)/vols[reg]/6.022*10
      
                np.savetxt(file_dend,out,comments="",header=header)

