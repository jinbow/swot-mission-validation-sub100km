# plot all swot data in a large plot 

# Load modules
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
xrod= xr.open_dataset
# Load data
def plot_mooring_locs(ax):
    ## load mooring positions
    fn='../data/mooring.data/mooring_positions.csv'
    positions=pd.read_csv(fn)
    #select mooring name has letter P
    p=positions[positions['name'].str.contains(f'P')]
    ax.scatter(p['lon'],p['lat'],color='red',marker='.',s=0.2,zorder=10)
    s=positions[positions['name'].str.contains(f'S')]
    ax.scatter(s['lon'],s['lat'],color='blue',marker='.',s=0.2,zorder=10)
    return 

def plot_karin(lon,lat,ax,keep_mean=False):
    if keep_mean:
        ax.pcolormesh(lon,lat,ssha_karin,cmap='Spectral_r',vmin=-0.1,vmax=0.1)
    else:
        ax.pcolormesh(lon,lat,ssha_karin-np.nanmean(ssha_karin,axis=0,keepdims=True),
                      cmap='Spectral_r',vmin=-0.05,vmax=0.05)

for keep_mean in [True,False]:
    for fnn in [13,26]:
        fn=f'../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_{fnn:03d}_sub_lat-30-40.nc'
        kk=xrod(fn)
        print(kk)
        nt=kk.cycle.size

        ntt=nt//10+1
        fig,ax=plt.subplots(ntt,10,figsize=(20,20),sharex=True,sharey=True,gridspec_kw={'hspace':0.05,'wspace':0.05})
        axs=ax.flatten()
        mk=(kk.latitude[:,20]>35)&(kk.latitude[:,20]<36.5)
        for i in range(nt):
            ssha_karin=(kk.ssha_karin+kk.height_cor_xover+kk.internal_tide_hret)[i,...].values[mk,:]
            ssha_karin_qual=kk.ssha_karin_qual[i,...].values[mk,:]
            ssha_karin=np.where(ssha_karin_qual==0,ssha_karin,np.nan)
            lon,lat=kk['longitude'][mk,:]-360,kk['latitude'][mk,:]
            if keep_mean:
                axs[i].pcolormesh(lon,lat,ssha_karin,cmap='Spectral_r',vmin=-0.1,vmax=0.1)
            else:
                ssha_karin=ssha_karin-np.nanmean(ssha_karin,axis=0,keepdims=True)
                axs[i].pcolormesh(lon,lat,ssha_karin,cmap='Spectral_r',vmin=-0.05,vmax=0.05)            
            if fnn==26:
                axs[i].plot(lon[:,16],lat[:,16],'-',color='gray',lw=0.5,alpha=0.9)
            time=kk['time'][i,200].item()*np.timedelta64(1,'s') + np.datetime64('2023-01-01')
            
            axs[i].text(0.02,0.9,f'{i}: {time}',fontsize=8,transform=axs[i].transAxes)
            plot_mooring_locs(axs[i])
            axs[i].set_ylim(35.1,36.7)
        plt.tight_layout()
        if keep_mean:
            plt.savefig(f'../figures/swot_images_pass_{fnn:03d}_w_mean.png',dpi=300)
        else:
            plt.savefig(f'../figures/swot_images_pass_{fnn:03d}_wo_mean.png',dpi=300)
        del fig, ax