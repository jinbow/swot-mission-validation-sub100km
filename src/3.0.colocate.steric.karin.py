# co-locate steric and karin 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from tqdm import tqdm 
from scipy.interpolate import griddata
from scipy import signal 

xrod=xr.open_dataset
ff=open('../data/3.0.colocate.steric.karin.log.txt',mode='w')

def colocate(karin_fn,steric):
    ff.write(f'run colocate({karin_fn},steric)')
    # clean steric, dropna for steric, and surface_time
    steric=steric.dropna(subset=['steric','surface_time'])
    
    ff.write("load swot data from {karin_fn}\n")
    kk=xrod(karin_fn)
    n0=len(kk)
    kk=kk.dropna(subset=['steric','surface_time'])
    ff.write(f'total number of rows {n0}, total usable rows {len(kk)}\n')
    #limit the latitude range to cover only the calval region
    mk=(np.nanmean(kk.latitude,axis=1)>34.5) & (np.nanmean(kk.latitude,axis=1)<36.5)
    dis=np.nanmean(kk.latitude,axis=1)-35.5
    ilat=np.where(np.abs(dis)==np.nanmin(np.abs(dis)))[0][0]    
    time=kk.time[:,ilat].values 
    time=time*np.timedelta64(1,'s')+np.datetime64('2023-01-01')
    
    msk=(time>=np.datetime64('2023-04-01')) & (time<np.datetime64('2023-07-12')) & (np.isnat(time)==False)
    time_karin=time[msk]
    ff.write(f"SWOT passing time {time_karin}\n")
    swh=kk.swh_model[msk,...].values[:,mk,:] #select the time and latitude range
    ssha_karin=100*(kk.ssha_karin+kk.height_cor_xover+kk.internal_tide_hret)[msk,...].values[:,mk,:] 
    
    ssha_karin_qual=kk.ssha_karin_qual[msk,...].values[:,mk,:]
    ssha_karin=np.where(ssha_karin_qual==0,ssha_karin,np.nan) #mask the bad data
    
    ssha_karin=ssha_karin-np.nanmean(ssha_karin,axis=0,keepdims=True) #remove the time mean
    ssha_karin=ssha_karin-np.nanmean(ssha_karin,axis=1,keepdims=True) #remove the alongtrack mean
    
    lat_karin=kk.latitude.values[mk,:]
    lon_karin=kk.longitude.values[mk,:]-360
    pass_num=int(karin_fn.split('pass_')[1][:3])
    
    for i in tqdm(range(len(time_karin)), desc='colocating steric height from all moorings with SWOT data'):
        ssha=ssha_karin[i,...].flatten()
        mk=(np.isfinite(ssha)) & (np.isfinite(lat_karin.flatten())) & (np.isfinite(lon_karin.flatten()))
        
        ssha=ssha[mk]
        lon1d=lon_karin.flatten()[mk]
        lat1d=lat_karin.flatten()[mk]
        swh1d=swh[i,...].flatten()[mk]
        
        ktime=np.datetime64(time_karin[i],'ns')
        
        mooring_ids=['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7']#,'ru32','ru38']
        for mid in mooring_ids:
    
            ss=steric[steric['Mooring_ID']==mid]
            #sort everything according to surface_time
            ss=ss.sort_values(by='surface_time')
            
            stime=np.datetime64('2023-01-01')+pd.to_timedelta(ss['surface_time'],'s')
            
            ss_d=xr.DataArray(ss['steric'],coords={'time':stime},dims=['time']).interp(time=ktime,method='nearest').values
            ss_d_linear=xr.DataArray(ss['steric'],coords={'time':stime},dims=['time']).interp(time=ktime,method='linear').values

            ss_lon=xr.DataArray(ss['lon'],coords={'time':stime},dims=['time']).interp(time=ktime,method='nearest').values
            ss_lat=xr.DataArray(ss['lat'],coords={'time':stime},dims=['time']).interp(time=ktime,method='nearest').values
            
            if (lon1d.size==0) or (lat1d.size==0) or np.isnan(ss_lon.item()) or np.isnan(ss_lat.item()):
                ssha0=np.nan
                swh0=np.nan
            else:
                ssha0=griddata((lon1d,lat1d),ssha,(ss_lon.item(),ss_lat.item()),method='linear')
                swh0=griddata((lon1d,lat1d),swh1d,(ss_lon.item(),ss_lat.item()),method='linear')
                
            delta_t = (stime-ktime)/np.timedelta64(1,'m')
            
            time_left=np.max(delta_t[delta_t<0])
            time_right=np.min(delta_t[delta_t>=0])

            dktime=(ktime-np.datetime64('2023-01-01'))/np.timedelta64(1,'s')
            dout0=np.r_[dktime,ss_lon,ss_lat,time_left,time_right,ssha0,ss_d,ss_d_linear,swh0,pass_num]
            dout0=dout0.reshape(1,10)
            dout0=pd.DataFrame(dout0,columns=['time_karin','lon','lat','time_delta_left','time_delta_right','ssha_karin','steric','steric_linear','swh','pass_num'])
            dout0=dout0.assign(**{'mooring_id':mid})
            try:
                dout=pd.concat([dout,dout0])
            except:
                dout=dout0.copy()
    return dout

def parse_arguments():
    import argparse
    """Parse command-line arguments for all constants."""
    parser = argparse.ArgumentParser(description="Steric height calculation depth parameters.")
    parser.add_argument("--bottom_depth", type=int, default=-500, help="Bottom depth (negative value in meters).")
    parser.add_argument("--top_depth", type=int, default=0, help="Top depth (negative value in meters).")
    
    return parser.parse_args()

if __name__=='__main__':
    """
    This script collocates the steric height from all moorings with the SWOT data.
    Example usage:
        python 3.0.colocate.steric.karin.py --bottom_depth=-500 --top_depth=0
    """
    args=parse_arguments()
    integration_depth=[args.bottom_depth,args.top_depth]
  
    print('integration depth: ',integration_depth)
    fn_mooring=f'../data/mooring.data/all_mooring_steric_heights_depth_{-args.bottom_depth:3d}.csv'
    fn_glider=f'../data/rutgers/ru32_ru38_steric_heights_depth_{-args.bottom_depth:3d}.csv'
    
    steric=pd.concat([pd.read_csv(fn_mooring),pd.read_csv(fn_glider)])
  
    dd13=colocate('../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_013_sub_lat-30-40.nc',steric)
    print('pass 013 is done',dd13)
    dd26=colocate('../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc',steric)
    print('pass 026 is done',dd26)

    ddd=pd.concat([dd13,dd26])
    print('concatenated results from two passes',ddd)
    ddd=ddd.dropna(subset=['steric','steric_linear','ssha_karin'])
    print('total number of collocated points: ',len(ddd))
    fnout=f'../data/3.0.colocated_data_karin_moorings_gliders_depth_{-args.bottom_depth:3d}.csv'
    print(f'save the collocated data to {fnout}')
    ddd.to_csv(fnout,index=False, float_format='%.5f',na_rep='NaN')
