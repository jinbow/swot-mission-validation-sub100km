# Calibrate gliders against S moorings
# Only use steric height as the bias will have an effect on steric height by a constant 

# Load modules
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import utils 


def select_along_karin_swath_center(data_fn,dis=0.025):
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
    lon,lat=lon[msk,:],lat[msk,:]

    df=pd.read_csv(data_fn)
    
    p_lon,p_lat=df['lon'],df['lat']
    #iterate every row of glider data and subselect the glider data with distance smaller than 0.05 degrees to the SWOT track given by lon and lat
    iid=[]
    for i in range(len(p_lon)):
        dist=np.sqrt((lon[:,1]-p_lon[i])**2+(lat[:,1]-p_lat[i])**2)
        if np.min(dist)<dis:
            iid.append(i)
    
    new_p=df.iloc[iid]
    #return the selected glider data and the SWOT track
    return new_p,lon,lat 

def plot_and_save_gliders_along_center():
    xrod= xr.open_dataset
    # Load data
    fn=f'../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc' # the descending pass
    kk=xrod(fn)
    lon=kk['longitude'][:,15:18]-360 # the center of the descending pass along the mooring array is at num_pixels=16
    lat=kk['latitude'][:,15:18]
    msk=(lat[:,1]>35)&(lat[:,1]<36.5)
    lon,lat=lon[msk,:],lat[msk,:]

    # load glider data 
    fn='../data/rutgers/ru32_ru38_steric_heights.csv'
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
    plt.plot(lon[:,1],lat[:,1],'-',lw=2,color='gray',label='SWOT')
    plt.plot(lon[:,0],lat[:,0],'--',lw=2,color='gray')
    plt.plot(lon[:,-1],lat[:,-1],'--',lw=2,color='gray')


    plt.plot(glider_lon,glider_lat,'.',markersize=4,color='m',label='Glider')
    utils.plot_mooring_locs(plt.gca())

    plt.show()

    new_gliders=gliders.iloc[iid]
    new_gliders.to_csv('../data/rutgers/ru32_ru38_steric_heights_along_center.csv',index=False)
    return new_gliders 

def cross_calibrate(integration_depth=[-420,0]):
    gliders=pd.read_csv(f'../data/rutgers/ru32_ru38_steric_heights_{-integration_depth[0]:3d}.csv')
    moorings=pd.read_csv(f'../data/mooring.data/all_mooring_steric_heights_{-integration_depth[0]:3d}.csv')

    ii=0
    for index, glider in gliders.iterrows():
        #for ss in ['P1','P2','P3','P4','P5','P6','P7']:
        #for ss in ['P1']:
        for ss in ['S2','S3','S4']:
            mooring=moorings[moorings['Mooring_ID']==ss]
            if len(mooring)>0:
                #time_mooring=(mooring['time_min']+mooring['time_max'])/2.0
                time_mooring=mooring['surface_time'].values
                #time_glider=(glider['time_min']+glider['time_max'])/2.0

                time_glider=glider['surface_time']
                dt = np.abs(time_mooring-time_glider)/60 #convert to minutes
                if dt.min()<60:
                    mk=np.argmin(dt)
                    #calculate the distance between mooring and glider at the mk index
                    dis=np.sqrt((mooring['lon'].iloc[mk]-glider['lon'])**2+(mooring['lat'].iloc[mk]-glider['lat'])**2)
                    if dis<0.01:
                        a=[glider['lat'],glider['lon'],glider['time_min'],glider['time_max'],
                            glider['steric'],glider['Mooring_ID'],mooring['lat'].iloc[mk],mooring['lon'].iloc[mk],
                            mooring['time_min'].iloc[mk],mooring['time_max'].iloc[mk],mooring['steric'].iloc[mk],
                            mooring['Mooring_ID'].iloc[mk] ]
                        if 'dout' not in locals():
                            dout=[a]
                        else:
                            dout.append(a)
        ii+=1 
    if 'dout' in locals():                          
        dout=pd.DataFrame(dout,columns=['glider_lat','glider_lon',
                                     'glider_time_min','glider_time_max',
                                     'glider_steric','glider_ID',
                                     'mooring_lat','mooring_lon',
                                     'mooring_time_min','mooring_time_max',
                                     'mooring_steric','mooring_ID'])
    else:
        raise ValueError('No matching mooring and glider data found.')
    return dout 
if __name__ == '__main__':
    #plot_and_save_gliders_along_center()
    dd=cross_calibrate()
    bias=(dd['mooring_steric']-dd['glider_steric']).mean()
    mooring_steric=dd['mooring_steric'].values-bias
    plt.scatter(mooring_steric,-dd['glider_steric'],s=10)
    print(np.std(mooring_steric+dd['glider_steric']))
    plt.plot([-0.1,0.1],[-0.1,0.1],'k--')
    plt.title(f'bias={np.std(mooring_steric+dd['glider_steric'])*100:.2f}')
    plt.show()