from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import h5py

title = {'axon1':'axon','axon2':'gc bottom','axon3':'gc top','ERmembrane12':'ER','ERmembrane21':'ER','ERmembrane11':'ER','ERmembrane22':'ER',}
n = 5
if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('No filename')
    f = h5py.File(sys.argv[1])
    grid = f['model']['grid']
    grid_no = f['model']['output']['__main__']['elements']
    species_idx = []
    for sp in f['model']['output']['__main__']['species']:
        species_idx.append(sp)
    res = f['trial0']['output']['__main__']['population']
    shape = (res.shape[0],res.shape[-1]+1)
    #region name is in -3 of every row grid_no
    
    # for g in grid:
    #     print g
    volumes = {}
    voxels = {}
    regions = []
    vox_in_row = {}
    volumes_in_row = {}
    for reg in f['model']['regions']:
        if reg != 'default' and 'ER' not in reg:
            regions.append(reg)
            
    fname = sys.argv[1].split('h5')[0]
    for reg in regions:
        voxels[reg] = {'submembrane1':[],'submembrane2':[]}
        volumes[reg] = {'submembrane1':0,'submembrane2':0}
        vox_in_row[reg] = {}
        volumes_in_row[reg] = {}
    for reg in regions:
        row = 0
        
        vox_in_row[reg][row] = {}
        volumes_in_row[reg][row] = {}
        
        for i,g in enumerate(grid):
            if g[15] == reg:
                if g[17] == 'submembrane':
                    if i:
                        if grid[i-1][17] == 'submembrane':
                            voxels[reg]['submembrane1'].append(i)
                            volumes[reg]['submembrane1']+=g["volume"]
                            vox = 0
                            vox_in_row[reg][row] = {}
                            volumes_in_row[reg][row] = {}
                        
                        else:
                            voxels[reg]['submembrane2'].append(i)
                            volumes[reg]['submembrane2']+=g["volume"]
                            
                            row+=1
                    else:
                        vox_in_row[reg][row] = {}
                        volumes_in_row[reg][row] = {}

                        voxels[reg]['submembrane1'].append(i)
                        volumes[reg]['submembrane1']+=g["volume"]
                        vox = 0
                else:
                    vox += 1
                    
                    vox_in_row[reg][row][vox] = i
                    volumes_in_row[reg][row][vox] = g["volume"]


    for reg in regions:
        if len(vox_in_row[reg]) == 1:
            for vox in vox_in_row[reg][0]:
                voxels[reg][vox] = [vox_in_row[reg][0][vox]]
                volumes[reg][vox] = volumes_in_row[reg][0][vox]
        else:
            mini = 2000000
            min_row = 0
            row_len = {}
            for row in vox_in_row[reg]:
                if len(vox_in_row[reg][row])<mini:
                    mini = len(vox_in_row[reg][row])
                    min_row = row
                row_len[row] = len(  vox_in_row[reg][row])
            the_same = True
            for  row in row_len:
                if row_len[row] != mini:
                    the_same = False
                    
            if not the_same:
                for vox in vox_in_row[reg][min_row]:
                    voxels[reg][vox] = []
                    volumes[reg][vox] = 0
                for  row in row_len:
                    if row_len[row] == mini:
                        diff = 0
                    else:
                        
                        diff = int((row_len[row]-mini)/2)
                    
                    first = [i+1 for i in range(diff)]
                    last = [row_len[row]-i for i in range(diff)]
                
                    for vox in vox_in_row[reg][row]:
                        
                        if vox in first:
                            which = 1
                        elif vox in last:
                            which = mini
                        else:
     
                            which = vox-diff

                            
                        voxels[reg][which].append(vox_in_row[reg][row][vox])
                        volumes[reg][which] += volumes_in_row[reg][row][vox] 
                            
            else:
                for  i,row in enumerate(row_len.keys()):
                   
                    if not i:
                        for vox in vox_in_row[reg][row]:
                            voxels[reg][vox] = [vox_in_row[reg][row][vox]]
                            volumes[reg][vox] = volumes_in_row[reg][row][vox]
                    else:
                        
                        for vox in vox_in_row[reg][row]:
                            voxels[reg][vox].append(vox_in_row[reg][row][vox])
                            volumes[reg][vox] += volumes_in_row[reg][row][vox]
                            
                    
        
    species = {'CK':['CKpCamCa4','CKpCamCa2','CKp','CKCamCa2','CKCamCa4'],'Epac':['Epac1cAMP'],'PKA':['PKAc','I1PKAc','PKAcPDE4','PKAc_PDE4_cAMP'],'PP2B':['PP2BCamCa4','Ip35PP2BCamCa4','Ip35PP1PP2BCamCa4'],'cAMP':['cAMP'], 'Ca':['Ca']}
    
    results = {}
    maxi = {'CK':12000,'Epac':800,'PKA':150,'PP2B':500, 'cAMP':1300, 'Ca':2500}
    mini = {'CK':0,'Epac':0,'PKA':0,'PP2B':0,'cAMP':0,'Ca':0}
    for i in species:
        results[i] = {}
        # maxi[i] = 0
        # mini[i] = 20000000
    
    #results['PP2BCaMKII'] = {}
    print results
    for k,reg in enumerate(regions):
        for fi in f:
            if 'trial' in fi:
                #             file_dend = sys.argv[1]+'_'+fi+'_conc_'+reg
                res = f[fi]['output']['__main__']['population']
                time = f[fi]['output']['__main__']['times']
                shape = (res.shape[0],len(voxels[reg]))
                print species
                for i,sp in enumerate(species):
                    out = np.zeros(shape)
                    if reg not in results[sp]:
                        results[sp][reg]= {}
                    for j,specie in enumerate(species[sp]):
                        idx = species_idx.index(specie)
                        for vox in voxels[reg]:
                            if vox == 'submembrane1':
                                out[:,0] += res[:,voxels[reg]['submembrane1'],idx].sum(axis=1)/volumes[reg]['submembrane1']/6.022*10
                            elif vox == 'submembrane2':
                                out[:,-1] += res[:,voxels[reg]['submembrane2'],idx].sum(axis=1)/volumes[reg]['submembrane2']/6.022*10
                            else:
                                out[:,vox] += res[:,voxels[reg][vox],idx].sum(axis=1)/volumes[reg][vox]/6.022*10
                   
                    results[sp][reg][fi] = out
            
            
              #   outPP2BCaMKII_ratio  = np.zeros((res.shape[0]/n,len(voxels[reg])))
                
    #             for j in range(res.shape[0]/n):
                    
    #                 outPP2BCaMKII_ratio[j,:] =   results['CK'][reg][j*n:(j+1)*n,:].mean(axis=0)/results['PP2B'][reg][j*n:(j+1)*n,:].mean(axis=0)
    #             if fi not in results['PP2BCaMKII'][reg][fi]:
    #                 results['PP2BCaMKII'][reg][fi] = outPP2BCaMKII_ratio

    #             if maxi['PP2BCaMKII'] < outPP2BCaMKII_ratio.max():
    #                 maxi['PP2BCaMKII'] = outPP2BCaMKII_ratio.max()
    #             if mini['PP2BCaMKII']> outPP2BCaMKII_ratio.min():
    #                 mini['PP2BCaMKII'] = outPP2BCaMKII_ratio.min()
    # print mini,maxi
    average_results = {}
    average_results['PP2BCaMKII'] = {}
    for sp in results:
        average_results[sp] = {}
        for reg in results[sp]:
            try:
                
                average_results[sp][reg] = sum([results[sp][reg][fi] for fi in results[sp][reg]])/len(results[sp][reg])
            except ValueError:
                average_results[sp][reg] = results[sp][reg]['trial0'].copy()
                i = 1
                for fi in results[sp][reg]:
                    if fi != 'trial0':
                        if len( average_results[sp][reg]) == len(results[sp][reg][fi]):
                            average_results[sp][reg]+= results[sp][reg][fi]
                            i += 1
    
                average_results[sp][reg] /= i
                                            
            # if maxi[sp] < average_results[sp][reg].max():
            #     maxi[sp] = average_results[sp][reg].max()
            # if mini[sp]> average_results[sp][reg].min():
            #     mini[sp] = average_results[sp][reg].min()
    mini['PP2BCaMKII'] = 0
    maxi['PP2BCaMKII'] = 500
    for reg in regions:
        outPP2BCaMKII_ratio  = np.zeros((int(average_results['CK'][reg].shape[0]/n),len(voxels[reg])))
                
        for j in range(int(average_results['CK'][reg].shape[0]/n)):
                    
            outPP2BCaMKII_ratio[j,:] =   average_results['CK'][reg][j*n:(j+1)*n,:].mean(axis=0)/average_results['PP2B'][reg][j*n:(j+1)*n,:].mean(axis=0)
            
        average_results['PP2BCaMKII'][reg] = outPP2BCaMKII_ratio

        # if maxi['PP2BCaMKII'] < outPP2BCaMKII_ratio.max():
        #     maxi['PP2BCaMKII'] = outPP2BCaMKII_ratio.max()
        # if mini['PP2BCaMKII']> outPP2BCaMKII_ratio.min():
        #     mini['PP2BCaMKII'] = outPP2BCaMKII_ratio.min()

    for sp in average_results:
        fig = plt.figure()
        ax = []
        out_name = fname+'_'+sp
        #cax = fig.add_axes([0.1, 0.1, 0.03, 0.8])
        print sp
        for k,reg in enumerate(regions):
            if reg == 'default':
                continue
            
            ax.append(fig.add_subplot(1,len(regions),k+1))
            ax[k].set_title(title[reg])
            im = ax[k].imshow(average_results[sp][reg],aspect='auto',cmap='CMRmap',vmin=mini[sp],vmax=maxi[sp])#,extent=[0,0,80.,2])
            if k:
                ax[k].axes.get_xaxis().set_ticklabels([])
                ax[k].axes.get_yaxis().set_ticklabels([])
            else:
                ax[k].axes.get_xaxis().set_ticklabels([])
                ax[k].axes.get_yaxis().set_ticklabels([-10,0,10,20,30,40,50,60,70,80])
                ax[k].set_ylabel('time (s)')
        cax = fig.add_axes([.1, 1., 0.8, .1])

        fig.colorbar(im,cax=cax, orientation='horizontal')
        plt.savefig(out_name+'.png',format='png', bbox_inches='tight',pad_inches=0.1,dpi=300)
        plt.savefig(out_name+'.eps',format='eps', bbox_inches='tight',pad_inches=0.1,dpi=300)
                

#    plt.show()
