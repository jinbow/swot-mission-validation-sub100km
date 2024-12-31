# profile statistics 

# import commonly used packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
xrod=xr.open_dataset


# import the data
steric=pd.read_csv('../data/mooring.data/all_mooring_steric_heights.csv')
steric_glider=pd.read_csv('../data/rutgers/ru32_ru38_steric_heights.csv')
# select steric according Mooring_ID='S1'
ss=pd.concat([steric,steric_glider])

steric_time=np.datetime64('2023-01-01')+pd.to_timedelta(steric['surface_time'],'s')

fig,ax=plt.subplots(1,1,figsize=(8,4))
#steric=steric[(steric['time_max']-steric['time_min'])<3600]
for a in ['S','P','ru']:
    msk=ss['Mooring_ID'].str.contains(a)
    print(f'{a} mooring has {msk.sum()} profiles.')
    deltat=(ss[msk]['time_max']-ss[msk]['time_min'])/60
    print(f'the mean profile time is {deltat[deltat<200].mean()} minutes')
    a,b=np.histogram(deltat,
                     bins=100,density=False,range=(0,100))
    
    ax.bar(b[:-1],a,width=b[1]-b[0],label=a)

ax.set_xlabel('Profiler time duration (min)')
ax.set_title('PDF of the Prawler profile duration')
ax.set_ylabel('Number of profiles')
ax.set_xlim(0,100)

plt.savefig('../figures/profile_duration.png',dpi=300)
plt.show()