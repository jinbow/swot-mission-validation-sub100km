# clean up and select colocated data 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import xarray as xr


def remove_trend(data,var_name,valid_points=5,
                 temporal=True,spatial=True):
    """ 
    Remove linear trend from data with NaNs 
    data: pandas dataframe 
    var_name: variable name to remove trend, 'steric_linear' or 'steric'. 
    valid_points: minimum number of valid points to remove trend, default is 5, less than 5 points will be dropped
    temporal: remove temporal trend, default is True
    spatial: remove spatial trend, default is True
    
    """
    def remove_linear(x,y,nn):
        if np.isfinite(y).sum()<nn:
            return y-np.nanmean(y)
        p=np.polyfit(x,y,1)
        yy=np.polyval(p,x)
        return y-yy
    #drop all rows with valid points less than valid_points
    data0=data.copy()
    unique_time=data0['time_karin'].unique()
    for i in range(len(unique_time)):
        ss=data0[data0['time_karin']==unique_time[i]]
        if len(ss)<valid_points:
            data0=data0.drop(data0[data0['time_karin']==unique_time[i]].index)
    for i in [1,2]: #loop through twice to remove the trend
        if temporal:    
            for i,mid in enumerate(mids):
                #print the rmsd 
                da=data0[data0['mooring_id']==mid].sort_values(by='time_karin')
                ss=remove_linear(da['time_karin'],da[var_name],10)
                da.loc[:,var_name]=ss
                kk=remove_linear(da['time_karin'],da['ssha_karin'],10)
                da.loc[:,'ssha_karin']=kk
                data0[data0['mooring_id']==mid]=da
        if spatial:
            unique_time=data0['time_karin'].unique()
            for i in range(len(unique_time)):
                ss=data0[data0['time_karin']==unique_time[i]]
                st=remove_linear(ss['lat'],ss[var_name],valid_points)
                ss.loc[:,var_name]=st
                st=remove_linear(ss['lat'],ss['ssha_karin'],valid_points)
                ss.loc[:,'ssha_karin']=st
                data0.loc[data0['time_karin']==unique_time[i]]=ss
    return data0


def select_along_karin_swath_center(df,dis=0.025):
    """ 
    data_fn: file name to the glider or mooring steric height data table
    """
    xrod= xr.open_dataset
    # Load data
    fn=f'../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc' # the descending pass
    kk=xrod(fn)
    lon=kk['longitude'][:,15:18]-360 # the center of the descending pass along the mooring array is at num_pixels=16
    lat=kk['latitude'][:,15:18]
    msk=(lat[:,1]>35)&(lat[:,1]<36.5)
    lon_karin,lat_karin=lon[msk,:],lat[msk,:]
    
    p_lon,p_lat=df['lon'].values,df['lat'].values
    #iterate every row of glider data and subselect the glider data with distance smaller than 0.05 degrees to the SWOT track given by lon and lat
    iid=[]
    for i in range(len(p_lon)):
        dist=np.sqrt((lon_karin[:,1]-p_lon[i])**2+(lat_karin[:,1]-p_lat[i])**2)
        if np.min(dist)<dis:
            iid.append(i)
    
    new_p=df.iloc[iid]
    #print(f'number of selected points: {len(iid)}')
    #return the selected glider data and the SWOT track
    return new_p,lon_karin,lat_karin
# load the colocated data

def select_colocated_data(dd,
                        max_deltat_left, # maximum time difference between steric and karin in minutes
                        max_deltat_right, # maximum time difference between steric and karin in minutes
                        distance_to_center, # maximum distance to the center of the swath in degrees
                        verbose=True):
    """ 
    Select colocated data from the input data table according to the criteria. 
    
    Parameters:
        dd: pandas dataframe for colocated data
        max_deltat_left: maximum time difference between steric and karin in minutes
        max_deltat_right: maximum time difference between steric and karin in minutes
        distance_to_center: maximum distance to the center of the swath in degrees
    Returns:
        ddd: selected colocated data
        lon_karin: longitude of the SWOT track
        lat_karin: latitude of the SWOT track    
    """
    data0=dd.copy()
    if verbose:
        ff.write(f'number of data points with gliders: {len(data0)} \n')
    data0= data0.drop(data0[data0['mooring_id'].str.contains('ru')].index)
    if verbose:
        ff.write(f'number of data points without gliders: {len(data0)} \n')

    #drop the rows where steric_linear and steric difference is larger than 0.5 cm
    data0= data0.drop(data0[np.abs(data0['steric_linear']-data0['steric'])>0.5].index)
    mk0=(data0['time_delta_right']<=max_deltat_right) | (np.abs(data0['time_delta_left'])<=max_deltat_left)
    ddd,lon_karin,lat_karin=select_along_karin_swath_center(data0[mk0],distance_to_center)
    if verbose:
        ff.write(f'number of data points with good steric: {len(ddd)} \n')

    return ddd,lon_karin,lat_karin

def plot_locations(dd,
                   max_deltat_left=60, # maximum time difference between steric and karin in minutes
                   max_deltat_right=60, # maximum time difference between steric and karin in minutes
                   distance_to_center=0.02, # maximum distance to the center of the swath in degrees
                   ):
    """ 
    Plot the locations of the colocated data.
    """
    fig,ax=plt.subplots(1,1,figsize=(2.5,5))
    #select data if mooring_id has ru 
    ddd=dd,lon_karin,lat_karin=select_colocated_data(dd,max_deltat_left,max_deltat_right,distance_to_center)

    gg=ddd[ddd['mooring_id'].str.contains('ru')]
    ax.plot(gg['lon'],gg['lat'],'g+',markersize=3,label='glider')
    pp=ddd[ddd['mooring_id'].str.contains('P')]
    ax.plot(pp['lon'],pp['lat'],'b+',markersize=3,label='P mooring')
    ss=ddd[ddd['mooring_id'].str.contains('S')]
    ax.plot(ss['lon'],ss['lat'],'r+',markersize=5,label='S mooring')
    
    ax.scatter(lon_karin,lat_karin,c='gray',s=3)
    ax.grid(True)
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    plt.ylim(35.2,36.3)
    plt.tight_layout()
    ff.write(f'save figure to ../figures/4.0.colocated_data_locations_{append_str}.png\n')

    plt.savefig(f'../figures/4.0.colocated_data_locations_{append_str}.png')

    return 

def plot_latitudes(dd,
                   max_deltat_left=60, # maximum time difference between steric and karin in minutes
                   max_deltat_right=60, # maximum time difference between steric and karin in minutes
                   distance_to_center=0.025, # maximum distance to the center of the swath in degrees
                   ):
    """ 
    Plot the latitudes of the colocated data.
    """
    fig,ax=plt.subplots(1,1,figsize=(9,3))
    #select data if mooring_id has ru 
    ddd,_,_=select_colocated_data(dd,max_deltat_left,max_deltat_right,distance_to_center)

    gg=ddd[ddd['mooring_id'].str.contains('ru')]
    ax.plot(gg['time_karin'],gg['lat'],'go',markersize=5,label='glider')
    pp=ddd[ddd['mooring_id'].str.contains('P')]
    ax.plot(pp['time_karin'],pp['lat'],'bx',markersize=5,label='P mooring')
    ss=ddd[ddd['mooring_id'].str.contains('S')]
    ax.plot(ss['time_karin'],ss['lat'],'r+',markersize=5,label='S mooring')
    
    ax.grid(True)
    ax.set_xlabel('Time')
    ax.set_ylabel('Latitude')
    plt.ylim(35.2,36.3)
    for i,mid in enumerate(['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7']):
        ax.text(np.datetime64('2023-03-28'),36.2-i*0.094,mid,color='black',fontsize=12,fontweight='bold')
    #plt.legend(loc='bottom left')
    plt.tight_layout()
    ff.write(f'save figure to ../figures/4.0.colocated_data_latitudes_time_{append_str}.png\n')

    plt.savefig(f'../figures/4.0.colocated_data_latitudes_time_{append_str}.png')

    return 

def plot_profiles_space(dd,var_name='steric_linear'):
    """
    Plot the spatial profiles of the colocated data. 
    """
    
    #ddd,_,_=select_colocated_data(dd,max_deltat_left,max_deltat_right,distance_to_center)
    #ddd=remove_trend(ddd,var_name,temporal=True,spatial=True)
    ddd=dd.copy()
    time=ddd['time_karin'].unique()
    nn=len(time)
    if nn//5*5==nn:
        n_row=nn//5
    else:    
        n_row=nn//5+1
    fig,ax=plt.subplots(n_row,5,figsize=(10/1.3,15/1.5),sharex=True,sharey=True,gridspec_kw={'hspace': 0.1, 'wspace': 0.1})
    axs=ax.flatten()
    for i in range(nn):
        sss=ddd[ddd['time_karin']==time[i]]
        ss=sss.sort_values(by='lat')
        axs[i].plot(ss['lat'],ss[var_name],'b-+')
        axs[i].plot(ss['lat'],ss['ssha_karin'],'r-x')
        axs[i].set_title('')
        axs[i].grid(True)
        axs[i].set_ylim(-3,3)
    ax[n_row//2,0].set_ylabel('ssha (cm)')
    ax[-1,2].set_xlabel('Latitude')
    ax[0,2].set_title('Mooring SHA (blue) and SWOT SSHA (red) spatial profiles')
    plt.tight_layout()
    figfn=f'../figures/4.0.colocated_data_steric_spatial_profiles_{append_str}.png'
    ff.write(f'save figure to {figfn}\n')
    plt.savefig(figfn,dpi=300)
    return

def remove_linear(x,y):
    ix=np.argsort(x)
    p=np.polyfit(x[ix],y[ix],1)
    yy=np.polyval(p,x[ix])
    return y[ix]-yy


def plot_time_series(data0,var_name):
    fig,ax=plt.subplots(11,1,figsize=(10,10),sharex=True,sharey=True,
                        gridspec_kw={'hspace': 0.1, 'wspace': 0.1})
    data=data0#remove_trend(data0,var_name,temporal=True,spatial=True)
    mids=['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7']
    ccc=[]
    for i,mid in enumerate(mids):
        #print the rmsd 
        da=data[data['mooring_id']==mid].sort_values(by='time_karin')
        ss=da[var_name]
        kk=da['ssha_karin']
        #calculate the coorelation coefficient
        cc=np.corrcoef(ss,kk)[0,1]
        print(f'correlation coefficients {mid} {cc}')
        ccc.append(cc)
        time=da['time_karin']*np.timedelta64(1,'s')+np.datetime64('2023-01-01')
        ax[i].plot(time,ss,'b-+',markersize=2)
        ax[i].plot(time,kk,'r-x',markersize=2)
        ax[i].text(np.datetime64('2023-03-28'),1.6,mid,color='black',fontsize=12,fontweight='bold')
        ax[i].set_ylabel('cm')
        ax[i].set_ylim(-3.1,3.1)
    print('correlation coefficient mean',np.array(ccc).mean(),np.array(ccc).std())
    plt.tight_layout()
    ax[0].set_title('Mooring SHA (blue) and SWOT SSHA (red) time series')
    plt.savefig(f'../figures/4.0.colocated_data_time_series_{append_str}.png',dpi=300)
    ff.write(f'save figure to ../figures/4.0.colocated_data_time_series_{append_str}.png\n')
    return 

def error_stats(data,var_name,plotit=False,verbose=True):
    """
    Calculate the error statistics for the colocated data.
    """
    
    mids=['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7']
    mean_diff=[];std=[];rmsd=[];num_points=[]
    for i, mid in enumerate(mids):
        #print the rmsd 
        da=data[data['mooring_id']==mid]
        ss=da[var_name]
        kk=da['ssha_karin']
        dif=np.abs(ss-kk) # absolute difference
        mean_diff.append(np.mean(dif)) # mean absolute difference
        std.append(np.std(dif)) # standard deviation of the difference 
        rmsd.append(np.var((ss-kk))**0.5) # root mean square difference
        num_points.append(len(ss))
    rmsdd=pd.DataFrame(np.r_[np.array(mean_diff).reshape(1,-1), 
                            np.array(std).reshape(1,-1),
                            np.array(rmsd).reshape(1,-1),
                            np.array(num_points).reshape(1,-1)],
                            index=['mad','std','rmsd','num_points'],columns=mids)
    rmsdd['total']=np.r_[np.mean(np.abs(data[var_name]-data['ssha_karin'])),
                        np.std(np.abs(data[var_name]-data['ssha_karin'])),
                        np.var((data[var_name]-data['ssha_karin']))**0.5,
                        np.sum(rmsdd.loc['num_points',:])]
    
    #rmsd['total']=np.r_[np.nanstd((data[var_name]-data['ssha_karin'])),np.sum(rmsd.loc['num_points',:])]
    #rmsd['max']=np.nanmax((data[var_name]-data['ssha_karin']))
    if verbose:
        ff.write(f'{5*' '}statistics \n{20*'-'} \n {rmsdd}\n')
    
    if plotit:
        fig=plt.figure(figsize=(5,3))
        plt.plot(rmsdd.columns,rmsdd.loc['mad',:],'o-',color='blue',label='Mean Absolute Difference')
        plt.errorbar(rmsdd.columns,rmsdd.loc['mad',:],yerr=rmsdd.loc['std',:] ,color='blue')
        plt.plot(rmsdd.columns,rmsdd.loc['rmsd',:],'o-',color='red',label='Root Mean Square Difference')
        plt.ylabel('cm')
        #plot legend to the top
        plt.grid(True)
        plt.legend()
        plt.ylim(0,1.2)
        ax1=plt.gca().twinx()
        ax1.plot(rmsdd.columns[:-1],rmsdd.loc['num_points',:][:-1],'d',color='y',markersize=5)
        ax1.set_ylabel('Number of points',color='y')
        ax1.set_ylim(0,120)
        ax1.set_yticks(np.arange(0,121,20))
        ax1.set_yticklabels(ax1.get_yticks(),color='y')
        plt.tight_layout()
        if verbose:
            ff.write(f'save figure to ../figures/4.0.colocated_data_rmsd_{append_str}.png\n')
        plt.savefig(f'../figures/4.0.colocated_data_rmsd_{append_str}.png',dpi=300)
        ### AD vs SWH
        plt.figure(figsize=(5,3))
        da=data.copy()
        ss=da[var_name]
        kk=da['ssha_karin']
        swh=da['swh']
        dif=np.abs(ss-kk)**2
        nn=len(swh)
        plt.scatter(swh,dif,s=5)
        #plot bin average
        swh_grid=np.arange(0.5,5.6,0.5)
        ad_grid=np.zeros(len(swh_grid))
        ad_std_grid=np.zeros(len(swh_grid))
        for i in range(len(swh_grid)):
            mk=(swh>=swh_grid[i]-0.5) & (swh<swh_grid[i]+0.5)
            if (np.sum(mk)>20):
                ad_grid[i]=np.nanmean(dif[mk])
                ad_std_grid[i]=np.nanstd(dif[mk])
            else:
                ad_grid[i]=np.nan
                ad_std_grid[i]=np.nan
        plt.errorbar(swh_grid,ad_grid,yerr=ad_std_grid,fmt='o',color='r')
        
        plt.xlabel('SWH (m)')
        plt.ylabel('MSD (cm$^2$)')
        plt.title('MSD vs SWH')
        plt.xlim(0,5.5)
        plt.ylim(0,3)
        plt.xticks(np.arange(0,5.6,0.5))
        plt.yticks(np.arange(0,3,1))
        plt.grid(True)
        plt.tight_layout()
        figfn=f'../figures/4.0.colocated_data_swh_rmsd_{append_str}.png'
        plt.savefig(figfn,dpi=300)
        if verbose:
            ff.write(f'save figure to {figfn}\n')
                    
        plt.figure(figsize=(6,6))
        da=data.copy()
        ss=da[var_name]
        kk=da['ssha_karin']
        plt.scatter(ss,kk,s=5,c='r',marker='+')
        plt.xlabel('SHA (cm)')
        plt.ylabel('SSHA_KaRIn (cm)')
        # 2d histogram 
        # hist, x_edges, y_edges = np.histogram2d(ss,kk, bins=[40,40])
        # plt.pcolormesh(x_edges, y_edges, hist.T, cmap='coolwarm')                           
        vv=3
        plt.plot([-vv,vv],[-vv,vv],'k--',zorder=0)
        plt.xticks(np.arange(-vv,vv+0.1,1))
        plt.yticks(np.arange(-vv,vv+0.1,1))
        plt.xlim(-vv,vv)
        plt.ylim(-vv,vv)
        plt.axis('equal')
        plt.grid(True)
        plt.tight_layout()
        figfn=f'../figures/4.0.colocated_data_scatter_karin_vs_steric_{append_str}.png'
        plt.savefig(figfn,dpi=300)
        if verbose:
            ff.write(f'save figure to {figfn}\n')
        
        # plot the CDF of rmsd 
        plt.figure(figsize=(5,3))

        dif=dif[np.isfinite(dif)]
        a,b=np.histogram(dif,30)
        plt.plot(b[:-1],a,'r-d',label='PDF')
        plt.ylabel('number of pairs',color='r')
        plt.ylim(0,a.max())
        plt.yticks(np.arange(0,a.max(),20),color='r')
        plt.xlabel('Absolute Difference (cm)')
        a=np.cumsum(a)/np.sum(a)
        ax1=plt.gca().twinx()
        ax1.plot(b[:-1],a*100,'b-o',label='CDF')
        ax1.set_title('PDF and CDF of mean absolute difference')
        ax1.set_xticks(np.arange(0,2.0,0.2))
        ax1.set_yticks(np.arange(0,110,10))
        ax1.tick_params(axis='y', colors='b')
        ax1.set_ylim(0,100)
        ax1.set_ylabel('CDF Percentage (%)',color='b')
        plt.xlim(0,2)
        plt.grid(True)
        plt.tight_layout()
        figfn=f'../figures/4.0.colocated_data_rmsd_cdf_{append_str}.png'
        plt.savefig(figfn,dpi=300)
        if verbose:
            ff.write(f'save figure to {figfn}\n')

    return rmsdd

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
    Main function to process colocated data and compute statistics.
    
    This function processes colocated data, calculates statistics, and handles 
    filtering based on user-specified parameters.

    Example usage:
        python 4.0.process.colocated.data.py --bottom_depth=-500 --top_depth=0 --deltat_threshold=100 --valid_points=9
    """
    #make ff global variable
    import sys
    global ff,mids,args,append_str
    mids=['S1','P1','P2','S2','P3','P4','S3','P5','P6','S4','P7']
    
    args=parse_arguments()
    deltat_threshold=args.deltat_threshold # maximum time difference between steric and karin in minutes
    valid_points=args.valid_points # minimum number of valid points to remove trend, default is 5; less than 5 points will be dropped
    # the final parameters used in the paper are: -500 0 100 9
    append_str=f'depth_{-args.bottom_depth:3d}.valid_points.{valid_points:02d}.deltat_{deltat_threshold:03d}m'
    log_fn=f'../data/4.0.log.colocated.data.stats_{append_str}.txt'
    
    plotit=True
    
    if True:
        print(f'write statistics to {log_fn}')
        ff = open(log_fn, 'w')
        ff.write(f'python {" ".join(sys.argv)}\n')
        
        max_deltat_left=deltat_threshold # maximum time difference between steric and karin in minutes
        max_deltat_right=deltat_threshold # maximum time difference between steric and karin in minutes
        distance_to_center=0.1 # maximum distance to the center of the swath in degrees
        var_name='steric_linear'    
        fn_in=f'../data/3.0.colocated_data_karin_moorings_gliders_depth_{-args.bottom_depth:3d}.csv'
        
        ff.write(f'load data from {fn_in}\nuse varname {var_name}\n')
        
        data0=pd.read_csv(fn_in)
        
        data0,_,_=select_colocated_data(data0,max_deltat_left,max_deltat_right,distance_to_center)

        data0=remove_trend(data0,var_name,valid_points=valid_points,
                        temporal=True,spatial=True)
        fn_out=f'../data/4.0.colocated_data_karin_moorings_gliders_{append_str}.clean.csv'
        data0.to_csv(fn_out,index=False) 
        #number of rows in data0
        rmsd=error_stats(data0,var_name,plotit)
        ff.write(f'number of rows in data0: {len(data0)}\n')
        ff.write(f'number of unique time: {len(data0['time_karin'].unique())}\n')
        tt=data0['time_karin'].unique()
        ff.write(f'number of days {(tt.max()-tt.min())/86400})\n')
        ff.write(f'save cleaned data to {fn_out}\n')
        
        if plotit:
            data1=plot_profiles_space(data0,var_name)
            data2=plot_time_series(data0,var_name)
    
    ff.close()
    
    
