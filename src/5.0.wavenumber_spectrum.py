# calcualte the wavenumber spectrum of steric and karin ssh 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import utils
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy import signal


def spectrum_karin_swath_center(method='welch'):
    
    """ 
    data_fn: file name to the glider or mooring steric height data table
    """
    from scipy import signal
    import xarray as xr
    
    xrod= xr.open_dataset
    # Load data
    fn=f'../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc' # the descending pass
    kk=xrod(fn)
    lon=kk['longitude'][:,15:18]-360 # the center of the descending pass along the mooring array is at num_pixels=16
    lat=kk['latitude'][:,15:18]
    msk=(lat[:,1]>35-0.5)&(lat[:,1]<36.5)
    ssha=100*(kk['ssha_karin']+kk['height_cor_xover'])[:,msk,16].values
    
    nt,ny=ssha.shape
    
    m=np.isnan(ssha).sum(axis=1)==0
    ssha=ssha[m,...]
    print(ssha.shape,ny)
    ssha=ssha-np.nanmean(ssha,axis=0,keepdims=True)
    ssha=ssha-np.nanmean(ssha,axis=1,keepdims=True)
    a,b=signal.welch(ssha.flatten(),fs=1/2,nperseg=ny,nfft=ny,noverlap=0,window='hann',detrend='linear')
    if method=='welch':
        return a[1:],b[1:]  
    elif method=='lombscargle':
        b=0
        for i in range(ssha.shape[0]):
            ss=signal.detrend(ssha[i,:])
            var=ss.var()
            b0=signal.lombscargle(np.arange(ny)*2,ss*np.hanning(ny),2*np.pi*a[1:],normalize=False) 
            b0*=var/(b0.sum()*(a[1]-a[0]))
            b+=b0 #*2*np.pi*len(ssha[i,:])
        b/=ssha.shape[0]
        return a[1:],b


def wavenumber_spectrum(da,dx,varn,extend_to=120,method='interpolate'):
    """ 
    x is the alongtrack distance
    y is the ssha data
    """
    from scipy import fftpack,signal,stats,interpolate
    ddx=utils.distance_between_points(da['lon'].values[0:-1],da['lon'].values[1:],da['lat'].values[0:-1],da['lat'].values[1:])/1e3
    if type(varn)==str:
        s1=da[varn].values
    elif type(varn)==type([]): #if varn is a list, then take the difference of the two variables
        s1=da[varn[0]].values-da[varn[1]].values
    dx2=(extend_to-100)/2 # dx2 km to the left and right of the data
    x=np.r_[0,dx2,ddx].cumsum() # extend the data to the left by dx2 km
    x=np.r_[x,extend_to] # extend the data to the right by dx2 km
    
    
    y=np.r_[0,s1,0] # pad the data with zeros
       
    if method=='interpolate':
        xx=np.arange(0,extend_to+1,dx) # dx-km resolution
        nn=xx.size # number of points
        yy=interpolate.interp1d(x,y)(xx)
        fs,lsxx=signal.welch(yy,fs=1/(xx[1]-xx[0]),nperseg=nn,nfft=nn,noverlap=0,
                        window='hann',detrend='linear')
        return fs, lsxx
    elif method=='lombscargle':
        ffs=1/extend_to; nyquist=1/20
        fs=np.linspace(1/100.,nyquist,10)
        hann=signal.windows.hann(extend_to+1)
        winx=interpolate.interp1d(np.linspace(0,extend_to,extend_to+1),hann)(x)
        lsxx=signal.lombscargle(x,y*winx,2*np.pi*fs,normalize=False)
        lsxx*=y.var()/(lsxx.sum()*(fs[1]-fs[0])) # normalize the spectrum to match the variance of the data
        return fs, lsxx
    return 

def calc_spectrum(data,dx=2,extend_to=110,method='lombscargle'):
    time=data['time_karin'].unique()
    print(f'number of snapshots {len(time)}')
    #mids=['S1', 'P1', 'P2', 'S2', 'P3', 'P4', 'S3', 'P5', 'P6', 'S4', 'P7']
    for i in tqdm(range(len(time)),desc='calculating wavenumber spectrum'):
        da=data[data['time_karin']==time[i]].sort_values(by='lat',ascending=False)
        fs,s=wavenumber_spectrum(da,dx,'steric_linear',extend_to,method)
        fs,k=wavenumber_spectrum(da,dx,'ssha_karin',extend_to,method)
        fs,d=wavenumber_spectrum(da,dx,['ssha_karin','steric_linear'],extend_to,method)
        if i==0:
            ss=s.copy();kk=k.copy();diff=d.copy()
        else:
            ss=np.vstack([ss,s])
            kk=np.vstack([kk,k])
            diff=np.vstack([diff,d])
            
    return fs,ss,kk,diff

def plot_spectrum(data,fs0=[1/110,1/10]):
    fs,ss,kk,diff=calc_spectrum(data,dx=2,method='lombscargle',extend_to=120)
    fm=(fs>=fs0[0]) & (fs<=fs0[1])
    print(1./fs[fm])
    fig,ax=plt.subplots(1,1,figsize=(6,4))
    
    ax.loglog(fs[fm],ss.mean(axis=0)[fm],'r-',lw=2,label='steric')
    ax.loglog(fs[fm],kk.mean(axis=0)[fm],'b-',lw=2,label='karin')
    ax.loglog(fs[fm],diff.mean(axis=0)[fm],'g-',lw=2,label='(steric-karin)')
    ax.loglog(fs[fm],diff.mean(axis=0)[fm]/2,'g--',lw=2,label='(steric-karin)/2')
    
    a,b=spectrum_karin_swath_center(method='lombscargle')
    plt.loglog(a,b,'k-',lw=2,label='karin (2km pixel)')
    #a,b=spectrum_karin_swath_center(method='welch')
    #plt.loglog(a,b,'k--',lw=2,label='karin swath center')    
    plt.xlim(1/200,1/10)
    plt.xlabel('Wavenumber (cycle/km, cpkm)')
    plt.ylabel('PSD (cm$^2$/cpkm)')
    plt.legend()
    
    ax.loglog(a,7.5+0.00125*a**(-2),'--',color='gray',lw=2,zorder=0,label='7.5+0.00125k$^{-2}$')
    ax.text(0.083,10,'7.5',color='gray',fontsize=12)

    #ax.loglog(a,(7.5+0.00125*a**(-2))/3,'--',color='gray',lw=1,zorder=0,label='(7.5+0.00125k$^{-2})/4$')
    #ax.text(0.083,10,'7.5',color='gray',fontsize=12)
    
    ax.loglog(a,2+0.00125*a**(-2),':',color='gray',lw=2,zorder=0,label='2.0+0.00125k$^{-2}$')
    ax.text(0.083,2.5,'2.0',color='gray',fontsize=12)
    plt.legend()

    ci_low_factor = 0.825
    ci_high_factor = 1.27
    ax.vlines(1/150,100*ci_low_factor,100*ci_high_factor,linestyle='-',color='gray',lw=2)
    ax.text(1/120, 100, "95%", color='gray', 
        ha='center', va='center', backgroundcolor='white',fontsize=15)
    ax.set_ylim(0.9,500)

    plt.tight_layout()
    print(f'save figure to ../figures/wavenumber_spectrum_{-args.bottom_depth}.pdf')
    plt.savefig(f'../figures/wavenumber_spectrum_{-args.bottom_depth}.pdf')
    plt.savefig(f'../figures/wavenumber_spectrum_{-args.bottom_depth}.png',dpi=150)
    return 

def parse_arguments():
    import argparse
    """Parse command-line arguments for all constants."""
    parser = argparse.ArgumentParser(description="colocation parameters.")
    parser.add_argument("--bottom_depth", type=int, default=-500, help="Bottom depth (negative value in meters).")
    parser.add_argument("--top_depth", type=int, default=0, help="Top depth (negative value or 0 in meters).")
    parser.add_argument("--deltat_threshold", type=int, default=100, help="Maximum allowable time difference between steric (mooring profile) and SWOT flyover time (in minutes).")
    parser.add_argument("--valid_points", type=int, default=9, help="Minimum number of valid points required to remove trends.")

    return parser.parse_args()

if __name__=='__main__':
    """ 
    compute the wavenumber spectrum of the steric height and karin-ssha and their difference.
    input arguments: integration_depth [bottom depth, top depth]
    
    Example usage:
    python 5.0.wavenumber_spectrum.py -500 0
    """
    import sys
    #make args global
    global args
    args=parse_arguments()
    fn_in=f'../data/4.0.colocated_data_karin_moorings_gliders_depth_{-args.bottom_depth:3d}.valid_points.{args.valid_points:02d}.deltat_{args.deltat_threshold:03d}m.clean.csv'

    data=pd.read_csv(fn_in)
    
    plot_spectrum(data,fs0=[1/110,1/5])
    
    plt.show()
    