# # Calculate density from moorings and gliders and put them in one file
import xarray as xr
from glob import glob 
import pylab as plt 
import numpy as np
import pandas as pd 


xrod=xr.open_dataset 

def c_rho(s, t, p, ispressure=True, lat=35.5, lon=-125, pr=0):
    """
    Calculate potential density from practical
    salinity, in-situ temperature, and pressure.

    Parameters
    ===========
    s: array-like
       Practical salinity
    t: array-like
       In-situ temperature (degrees Celsius)
    p: array-like
       Pressure (dbar) or depth (meters)
    ispressure: bool, optional (default=True)
       If False, interpret `p` as depth (meters) and convert to pressure
    lat: float, optional (default=36.7)
       Latitude of the measurement (degrees)
    lon: float, optional (default=-122.0335)
       Longitude of the measurement (degrees)
    pr: float, optional (default=0)
       Reference level pressure for potential density calculation (dbar)

    Returns
    =======
    rho: array-like
       Potential density (kg/m^3)
    """
    import gsw
    import numpy as np

    p = np.abs(p)
    if not ispressure:
        p = gsw.p_from_z(-p, lat)
    sa = gsw.SA_from_SP(s, p, lon, lat)
    rho = gsw.pot_rho_t_exact(sa, t, p, pr)

    return rho

def process_moorings(fnout):
    plotit=False
    if plotit:
        fig,ax=plt.subplots(7,1,figsize=(10,10),sharey=True,sharex=True)
    dout={}
    for ii in range(1,8):
        files=sorted(glob(f"{path_input}/level-2/*P{ii:1d}_CTD*nc")) 
        for i,fn in enumerate(files):
            d=xr.open_dataset(fn)
            #display(d.keys())
            msk=(d['TEMP_QC']==0)&(d['PSAL_QC']==0)
            temp=np.where(msk,d['TEMP'],np.nan).flatten()
            salt=np.where(msk,d['PSAL'],np.nan).flatten()
            pres=np.where(msk,d['PRES'],np.nan).flatten()
            #depth=np.where(msk,d['DEPTH'],np.nan)
            density=c_rho(salt,temp,pres)
            #print(np.nanmin(density),np.nanmax(density))
            time=d['TIME'].data.flatten()
            if plotit:
                cs=ax[ii-1].scatter(time[::30],-pres[::30],c=density[::30]-1000,s=0.5,cmap=plt.cm.jet,vmin=22,vmax=27)
                ax[ii-1].set_ylabel(f"P{ii:d}")
            msk=np.isfinite(pres)&np.isfinite(density)&np.isfinite(time)
            pres=pres[msk]
            density=density[msk]
            time=time[msk]
            #load the buoy location information 
            locs=pd.read_csv(f'{path_input}/PMEL_GPS_latlon_hrly/SWO{ii:1d}_hrly.txt',parse_dates=['time'])
            locs.set_index('time', inplace=True)
            lon=np.interp(time.astype('int64'), locs.index.astype('int64'), locs['longitude'])
            lat=np.interp(time.astype('int64'), locs.index.astype('int64'), locs['latitude'])
            plt.scatter(lon[::50],lat[::50],s=2)
            #print(lat)
            if 'FIXED' in fn:
                dd=np.r_[lat.reshape(1,-1),lon.reshape(1,-1),pres.reshape(1,-1),density.reshape(1,-1)]
                dout[f"P{ii:1d}_F"]=xr.DataArray(data=dd,dims=("variable",f"time_P{ii:1d}_F"),
                                                coords={"variable":["lat",'lon',"pressure","density"],f"time_P{ii:1d}_F":time},
                                                attrs={'comment':f"in-situ density from fixed CTD on mooring P{ii:1d}"})
            if 'PROF' in fn:
                dd=np.r_[lat.reshape(1,-1),lon.reshape(1,-1),pres.reshape(1,-1),density.reshape(1,-1)]
                dout[f"P{ii:1d}_P"]=xr.DataArray(data=dd,dims=("variable",f"time_P{ii:1d}_P"),
                                                coords={"variable":["lat",'lon',"pressure","density"],f"time_P{ii:1d}_P":time},
                                                attrs={'comment':f"in-situ density from Prawler CTD on mooring P{ii:1d}"})
        print(f"P{ii:1d}",dd.shape)

    #S moorings
    plotit=False
    fig,ax=plt.subplots(4,1,figsize=(10,8),sharey=True,sharex=True)

    for ii in range(1,5):
        files=sorted(glob(f"{path_input}/level-2/*S{ii:1d}*nc")) 
        for i,fn in enumerate(files):
            d=xr.open_dataset(fn)   
            if 'PROFILER' in fn:
                pn=d['PROF_NUM'].values.flatten()
                pn=np.r_[0,np.diff(pn)]
                pn[pn!=0]=1 
                pn=np.cumsum(pn) #profiler number
                
                #display(d)
                #msk=(d['QC_FLAG']<=1)
                time=d['TIME'].values.flatten()
                lat=d['LATITUDE_SURFACE_BUOY'].interp(TIME_MIDPROFILE=time).values
                lon=d['LONGITUDE_SURFACE_BUOY'].interp(TIME_MIDPROFILE=time).values
                pres=d['PRES'].values.flatten()
                salt=d['PSAL'].values.flatten()
                temp=d['TEMP'].values.flatten()
                
                density=c_rho(salt,temp,pres)
                
                #print(np.nanmin(density),np.nanmax(density))
                if plotit:
                    cs=ax[ii-1].scatter(time[::50],-pres[::50],c=density[::50]-1000,s=0.5,cmap=plt.cm.jet,vmin=22,vmax=27)
                    ax[ii-1].set_ylabel(f"S{ii:d}")
                
                
                msk=np.isfinite(pres)&np.isfinite(density)&np.isfinite(time)&(d['QC_FLAG'].values.flatten()==0)
                
                dd=np.r_[lat[msk].reshape(1,-1),lon[msk].reshape(1,-1),pres[msk].reshape(1,-1),density[msk].reshape(1,-1)]
                dout[f"S{ii:1d}_P"]=xr.DataArray(data=dd,dims=("variable",f"time_S{ii:1d}_P"),
                                                coords={"variable":['lat','lon',"pressure","density"],f"time_S{ii:1d}_P":time[msk]},
                                                attrs={'comment':f"in-situ density from WireWalker CTD on mooring S{ii:1d}"})
                
            else:
                msk=(d['TEMP_QC']==0)&(d['PSAL_QC']==0)
                temp=np.where(msk,d['TEMP'],np.nan)
                salt=np.where(msk,d['PSAL'],np.nan)
                pres=np.where(msk,d['PRES'],np.nan)
                #depth=np.where(msk,d['DEPTH'],np.nan)
                density=c_rho(salt,temp,pres)
                #print(np.nanmin(density),np.nanmax(density))
                time=d['TIME'].values
                nd,nt=temp.shape
                if plotit:
                    for j in range(nd):
                        ax[ii-1].scatter(time,-pres[j,:],c=density[j,:]-1000,s=0.5,cmap=plt.cm.jet,vmin=22,vmax=40)
                    ax[ii-1].set_ylabel(f"S{ii:d}")
                lat=d['LATITUDE_SURFACE_BUOY'].data.flatten()
                lon=d['LONGITUDE_SURFACE_BUOY'].data.flatten()
                
                for j in range(nd):
                    
                    msk=np.isfinite(pres[j,:])&np.isfinite(density[j,:])&np.isfinite(time[:])&np.isfinite(lon)&np.isfinite(lat)
                    
                    dd=np.r_[lat[msk].reshape(1,-1),lon[msk].reshape(1,-1),pres[j,:][msk].reshape(1,-1),density[j,:][msk].reshape(1,-1)]
                    
                    dout[f"S{ii:1d}_F_{j+1:1d}"]=xr.DataArray(data=dd,dims=("variable",f"time_S{ii:1d}_F_{j+1:1d}"),
                                                    coords={"variable":['lat','lon',"pressure","density"],f"time_S{ii:1d}_F_{j+1:1d}":time[msk]},
                                                    attrs={'comment':f"in-situ density from the {j+1:1d}th fixed CTD on mooring S{ii:1d}"})
        print(f"S{ii:1d}",dd.shape)    

    dout=xr.Dataset(dout)
    dout.attrs={'history':"Created by the script 1.density_all_moorings_gliders_level-2.py",
                'creator':'J. Wang',
                'time':str(np.datetime64('now'))}
    print(dout)
    dout.to_netcdf(fnout)
    return dout

def remove_outliers(fn_in):
    df=xr.open_dataset(fn_in)
    z_range = np.arange(-550,1,5)
    mids=['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7']

    for mid in mids:
        dff=np.abs(df[f'{mid}_P'].sel(variable='pressure').values)
        dens=df[f'{mid}_P'].sel(variable='density').values
        for z0 in z_range:
            zmin,zmax=z0,z0+5
            msk= (-dff>zmin) & (-dff<zmax)
            
            if msk.sum()>100:
                den=dens[msk]
                dena=den-np.nanmean(den)
                den_sigma=np.nanstd(dena)
                mk=np.abs(dena)>5*den_sigma
                den[mk]=np.nan
                df[f'{mid}_P'].data[3,msk]=den
                print(f"remove {mk.sum()} outliers, from {msk.sum()} data points at depth {z0} for mooring {mid}")
                #print(np.isnan(dens[msk]).sum(), np.isnan(den).sum(),np.isnan(df[f'{mid}_P'].data[3,msk]).sum())
        #df[f'{mid}_P']=dff.T
    df.to_netcdf(fn_in.replace('.nc','_removed_outliers.nc'))
    print(f'save data to {fn_in.replace(".nc","_removed_outliers.nc")}')
    del df 
    return
    
if __name__=='__main__':
    
    path_input="../data/mooring.data/"
    
    do_moorings=True
    do_remove_outliers=True
    fn='../data/mooring.data/density_all_moorings_level-2.nc'
    if do_moorings:
        dout=process_moorings(fn)
        print(dout)
    if remove_outliers:
        remove_outliers(fn)