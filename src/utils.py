import numpy as np
import xarray as xr
import pandas as pd
from glob import glob
xrod=xr.open_dataset

def plot_science_requirement(ax,k,color='gray',lw=2,grid_size=7.5):
    """
    Plot the SWOT science requirement in a log-log plot.

    Parameters:
    -----------
    ax: matplotlib.axes.Axes object
        Axes object to plot the scientific requirement.
    k: array_like
        Array of wavenumber frequency [cycles/km] to be plotted.
    grid_size: float, optional
        Cross-swath distance in averaging
    Returns:
    --------
    None
    """
    ax.loglog(k,2/grid_size*7.5+0.00125*k**(-2),color=color,lw=lw)
    #ax.vlines(1/15,1e1,1e3,color='gray')
    #ax.text(1/15,1e3,'15 km',fontsize=16)
    return


def distance_between_points(lon0, lons, lat0, lats):
    import numpy as np
    # Convert latitude and longitude to 
    # spherical coordinates in radians.
    degrees_to_radians = np.pi/180.0
        
    # phi = 90 - latitude
    phi1 = lat0*degrees_to_radians
    phi2 = lats*degrees_to_radians
    dphi = phi1-phi2
    
    # theta = longitude
    theta1 = lon0*degrees_to_radians
    theta2 = lons*degrees_to_radians
    dtheta=theta1-theta2
    # Compute spherical distance from spherical coordinates.
        
    # For two locations in spherical coordinates 
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) = 
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length
    
    #co = (np.sin(phi1)*np.sin(phi2) + 
    #       np.cos(phi1)*np.cos(phi2)*np.cos(theta1-theta2))
    #arc = np.arccos( co )
    
    #The haversine formula
    co = np.sqrt(np.sin(dphi/2)**2 + np.cos(phi1)*np.cos(phi2)*np.sin(dtheta/2.0)**2)
    arc = 2* np.arcsin(co)
    dist = arc*6371.0e3
    
    return dist


def plot_glider_locs(ax):
    # Load modules
    import numpy as np
    import xarray as xr
  
    import matplotlib.pyplot as plt
    import utils 

    xrod= xr.open_dataset
    # Load data
    fn=f'../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc' # the descending pass
    kk=xrod(fn)
    lon=kk['longitude'][:,15:18]-360 # the center of the descending pass along the mooring array is at num_pixels=16
    lat=kk['latitude'][:,15:18]
    msk=(lat[:,1]>35)&(lat[:,1]<36.5)
    lon,lat=lon[msk,:],lat[msk,:]
    # load glider data 
    fn='../data/rutgers/ru32_ru38_steric_heights_depth_500.csv'
    gliders=pd.read_csv(fn)
    #gliders['time']=pd.to_datetime(gliders['time'])
    #gliders=gliders.set_index('time')
    glider_lon,glider_lat=gliders['lon'],gliders['lat']
    #iterate every row of glider data and subselect the glider data with distance smaller than 0.05 degrees to the SWOT track given by lon and lat
    iid=[]
    for i in range(len(glider_lon)):
        dist=np.sqrt((lon[:,1]-glider_lon[i])**2+(lat[:,1]-glider_lat[i])**2)
        if np.min(dist)<0.025:
            iid.append(i)
    glider_lon,glider_lat=glider_lon[iid],glider_lat[iid]
    ax.plot(lon[:,1],lat[:,1],'-',lw=2,color='gray',label='SWOT')
    ax.plot(lon[:,0],lat[:,0],'--',lw=2,color='gray')
    ax.plot(lon[:,-1],lat[:,-1],'--',lw=2,color='gray')
    ax.plot(glider_lon,glider_lat,'.',markersize=4,color='m',label='Glider')
    return 


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

def remove_linear_trend_with_nan(datain):
      """ 
      remove the linear trend from the data with nan values
      Parameters
      ----------
      datain : numpy.ndarray
         input data with nan values
   
      Returns
      -------
      dataout : numpy.ndarray
         data with linear trend removed
      """
      import numpy as np
      
      msk=np.isfinite(datain)
      nn=datain.size
      if msk.sum()<5:
         return np.ones_like(datain)*np.nan
      else:  
         x=np.arange(nn)
         a = np.polyval(np.polyfit(x[msk],datain[msk],1),x)      
         return datain-a 
     

pass_num='026'
dt64=np.datetime64

def get_deep_steric(mn):
    fn = glob(f'/mnt/flow/swot/calval/Data/insitu_mooring/JPL_QC/FIXED_CTD/*-{mn}_CTD*.nc')

    d=xrod(fn[0])

    sh = d['STERIC_HEIGHT_BY_LAYER'].sel(TIME=slice('2023-04-05','2023-07-10'))
    if mn=='S1':
        aa=np.zeros_like(sh[:,0]) #np.nanmean(sh[:,:1],axis=1)
    else:
        aa=np.nanmean(sh[:,:3],axis=1)
    sh = xr.DataArray(aa,dims=['TIME'],coords={'TIME':sh.TIME})
                   
    return sh
#def clean_up_time(dain,time_dim_name):

def get_karin(lon,lat,ssh_varname='',fix_loc=True,karin_version='ver02'):
    """
    interpolate the karin data to the lon,lat,time

    Parameters
    ----------
    lon : xarray.DataArray
        longitude of the mooring with time as the dimension
    lat : xarray.DataArray
        latitude of the mooring with time as the dimension
    fix_loc : bool, optional
        if True, use the nearest karin data to the fixed mooring location. The default is True.

    Returns
    -------
    karin : numpy.ndarray
        karin data interpolated to the mooring location
    timec : numpy.ndarray
        time of the karin data
    """
    from scipy.interpolate import LinearNDInterpolator
    #print('lon,lat std',np.nanstd(lon),np.nanstd(lat))
    # load the karin data
    d=xrod(f'{swot_subset_path}/data_cube_Expert_{pass_num}_california_bias_removed_{karin_version}.nc')
    if pass_num=='026':
        i0,i1=0,35
    else:
        i0,i1=35,70
    lons,lats=d['longitude_avg_ssh'][1,:,i0:i1].data.flatten(),d['latitude_avg_ssh'][1,:, i0:i1 ].data.flatten()
    karin=d[ssh_varname].sel(timec=slice('2023-04-05','2023-07-10'))[:,:,i0:i1]
    timec=d['timec'].sel(timec=slice('2023-04-05','2023-07-10'))
    
    karin=karin.data.reshape(-1,lons.size)
    #karin -= np.nanmean(karin,axis=0,keepdims=True)

    #subset the karin data to the region of interest by lat 35N-36.4N enclosed by the mooring locations
    msk=(lats>35)&(lats<36.4)
    lons,lats=lons[msk],lats[msk]
    msk0=np.isnan(lons)|np.isnan(lats)
    lons,lats=lons[~msk0],lats[~msk0]
    karin=karin[:,msk][:,~msk0]
    
    # use the number of valid points as the quality control
    msk = (np.isfinite(karin).sum(axis=1)<2000)
    #var = np.nanstd(karin,axis=1) 
    #msk = var>threshold
    print(30*'*','Total number of valid snapshots',(~msk).sum())
    karin[msk,:]=np.nan

    dis=np.sqrt((lons-np.nanmean(lon))**2+(lats-np.nanmean(lat))**2)

    if fix_loc:
        karin = karin.reshape(-1,lons.size)[:,dis.argmin()]
    else:
        aa=[]
        for i in range(karin.shape[0]):
            msk2 = np.isfinite(karin[i,:])
            if msk2.sum()>2000:
                lo=lon.interp(time=timec[i],method='nearest').data
                la=lat.interp(time=timec[i],method='nearest').data
                aa.append(LinearNDInterpolator(np.c_[lons[msk2], lats[msk2]],karin[i,msk2])((lo,la)) )
            else:
                aa.append(np.nan)
        karin=np.array(aa)
    #
    #karin[var>0.05]=np.nan
    #print(karin, d['timec'].data.size)
    
    return karin,timec.data

def get_steric(fn,fix_loc=True,with_deep=False,ax=None,
               ssh_varname='ssh_karin_2_bias_removed',
               karin_version='ver02',):
    """ 
    get the steric height from the mooring data
    Parameters
    ----------
    fn : str
        filename of the mooring data
    fix_loc : bool, optional
        if True, use the nearest karin data to the fixed mooring location. The default is True.
    with_deep : bool, optional
        if True, add the deep steric height. The default is False.
    ax : matplotlib.axes, optional
        if not None, plot the steric height and karin data. The default is None.

    Returns
    -------
    sh : numpy.ndarray
        steric height
    karin : numpy.ndarray
        karin data interpolated to the mooring location
    lon : numpy.ndarray
        longitude of the mooring
    lat : numpy.ndarray
        latitude of the mooring
    t : datetime64 array
        time of the mooring data
    sh_out : xarray.DataArray
        steric height as xarray.DataArray
    """
    import xarray as xr
    xrod = xr.open_dataset
    d=xrod(fn)
    
    reference_date = np.datetime64('1950-01-01T00:00:00')  # The time in Luke's file is the days since 01/01/1950
    time = d['TIME_STERIC_HEIGHT'].data

    tt=reference_date+np.timedelta64(1,'s')*(d['TIME_STERIC_HEIGHT'].data*86400).astype('int')
    time = pd.to_datetime(tt)
    time=np.where(time!=0,time,np.datetime64('NaT'))
    msk = np.isnat(time)|np.isnan(d['STERIC_HEIGHT_ANOMALY'].data)

    time=time[~msk]
    ii = np.argsort(time)
 
    sh = d['STERIC_HEIGHT_ANOMALY'].data[~msk][ii]
    lon,lat=d['LONGITUDE_GPSSB_STERIC_HEIGHT'].data[~msk][ii]+360,d['LATITUDE_GPSSB_STERIC_HEIGHT'].data[~msk][ii]
    amin,amax=np.abs(d['STERIC_HEIGHT_MIN_PROF_DEPTH'].data[~msk][ii]),np.abs(d['STERIC_HEIGHT_MAX_PROF_DEPTH'].data[~msk][ii])

    #d=d.assign_coords(TIME_STERIC_HEIGHT=time)
    #d=d.dropna(dim='TIME_STERIC_HEIGHT',how='all')
    #d=d.sortby('TIME_STERIC_HEIGHT')
    #print(30*'-')
    #d=d.sel(TIME_STERIC_HEIGHT=slice('2023-04-05','2023-07-10'))
    #print(np.isnat(d['TIME_STERIC_HEIGHT']).sum())
    #amin,amax=np.abs(d['STERIC_HEIGHT_MIN_PROF_DEPTH'].data),np.abs(d['STERIC_HEIGHT_MAX_PROF_DEPTH'].data)
    thickness=amax-amin
    amin[np.isnan(amin)]=20 
    thickness[np.isnan(thickness)]=400
    msk= (amin>20)|(thickness <400)
    msk[np.isnan(msk)]=True
    time=time[~msk];sh=sh[~msk];lon=lon[~msk];lat=lat[~msk]
    
    #sh = np.where(msk,np.nan,d['STERIC_HEIGHT_ANOMALY'].data)
    #sh=sh[~mskt]
    #time=time[~mskt]

    
    mn = fn.split('/')[-1].split('-')[1][:2]
  
    if with_deep: # add the deep steric height
        if 'S' in mn: #only SIO moorings have deep steric height before PMEL recovery
            sh_deep=get_deep_steric(mn).interp(TIME=time)
            sh+=sh_deep

    sh_out=xr.DataArray(sh,dims=['time'],coords={'time':time}) #steric height
    lat=xr.DataArray(lat,dims=['time'],coords={'time':time}) #steric height latitude
    lon=xr.DataArray(lon,dims=['time'],coords={'time':time}) #steric height longitude


    if ax != None: #plot the steric height and karin data if ax is not None
        if np.isfinite(sh).sum()>110: #only plot if there are more than 110 snapshots

            karin,t_karin=get_karin(lon,lat,fix_loc=fix_loc, #whether fixed location or not 
                              ssh_varname=ssh_varname,
                              karin_version=karin_version,
                              threshold=karin_threshold,)
            msk = (t_karin>dt64('2023-04-05'))&(t_karin<dt64('2023-07-10')) #only use the data between April 5 and July 10
            t_karin=t_karin[msk]
            karin=karin[msk]
            karin-=np.nanmean(karin) #remove the time mean of karin data
            
            msk=(time>dt64('2023-04-05'))&(time<dt64('2023-07-10')) #perform the same time selection for the mooring data
            time=time[msk]
            sh=sh[msk]
            sh-=np.nanmean(sh)
            sh=xr.DataArray(sh,dims=['time'],coords={'time':time})
            
            ax.plot(time,sh,lw=0.3,color='gray',
                    label=mn,zorder=0)

            #interpolate steric to KaRIn time
            sh=sh.interp(time=t_karin)
            dif = np.nanmean(sh-karin)

            karin+=dif
            ax.scatter(t_karin,sh,s=20,c='r',marker='+')
            ax.scatter(t_karin,karin,s=20,c='b',marker='+',
                       facecolors='none',edgecolors='b',
                       zorder=1,linewidth=1)
            dif = sh-karin
            print(mn,)
            ax.legend()
    return  sh.data,karin,lon.interp(time=t_karin),lat.interp(time=t_karin),t_karin,sh_out



def remove_mean_density(pp,density,lat0=35.5,is_depth=False):
    """
    Remove the time-mean density to the mooring depth, skipping NaNs.

    Parameters:
    pp : pressure or depth (depth is negative, pressure is positive)
    density: in-situ density 
    is_depth : bool
        True if pp is depth, False if pp is pressure.

    Returns:
    numpy.ndarray
        density anomaly
    """
    import numpy as np
    from scipy.interpolate import interp1d
    import gsw
    import xarray as xr
    dmean=xr.open_dataset('../data/density_timemean.nc')

    # Extract data from dataset
    if is_depth:
        pp = -np.abs(pp)
    else:
        pp = gsw.z_from_p(np.abs(pp), lat0)

    # Mask NaNs in time-mean density data
    msk = np.isfinite(dmean['density_timemean_S'].values[:])
    # Interpolate the time-mean density profile to the mooring depths
    ff = interp1d(-dmean['z'][msk], dmean['density_timemean_S'][msk], fill_value='extrapolate')
    # Compute the difference between the actual density and the interpolated density
    den_anomaly = density - ff(pp)

    return den_anomaly
