# profile statistics 

# import commonly used packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
xrod=xr.open_dataset

def stats(data,mid,bottom_depth):
    nlines=len(data.copy())
    data=data.dropna(subset=['steric'])
    nlines_valid=len(data)
    
    mdt=np.mean(data['time_max']-data['time_min'])/60 
    npt=np.mean(data['num_points'])
    mmin=np.mean(data['depth_min'])
    mmax=np.mean(data['depth_max'])
    dd=pd.DataFrame([[mid,nlines,nlines_valid/nlines,mdt,npt,mmin,mmax,bottom_depth]],
                    columns=['ID','total_profiles','valid_profiles','mean_delta_T','n_obs',
                             'mean_bottom_z','mean_surface_z','bottom_depth'])
    return dd 

for bottom_depth in [-500, -400, -300]:
    # import the data
    steric=pd.read_csv(f'../data/mooring.data/all_mooring_steric_heights_depth_{-bottom_depth:3d}.csv')
    steric_glider=pd.read_csv(f'../data/rutgers/ru32_ru38_steric_heights_depth_{-bottom_depth:3d}.csv')
    # select steric according Mooring_ID='S1'
    ss=pd.concat([steric,steric_glider])

    steric_time=np.datetime64('2023-01-01')+pd.to_timedelta(steric['surface_time'],'s')

    #fig,ax=plt.subplots(1,1,figsize=(8,4))
    #steric=steric[(steric['time_max']-steric['time_min'])<3600]
    mooring_ids=['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7','ru32','ru38']

    for a in mooring_ids:
        da=stats(ss[ss['Mooring_ID']==a],a,bottom_depth)
        if a=='S1' and bottom_depth==-500:
            dout=da.copy()
        else:
            dout=pd.concat([dout,da])
print(dout.sort_values(by=['ID','bottom_depth']))
dout.to_csv('../data/2.1.profile.statistics.output.csv',float_format='%7.2f')