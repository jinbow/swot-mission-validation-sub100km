import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import gsw
import pandas as pd
#import dask
#import dask.dataframe as dd
#from dask import delayed
from tqdm import tqdm
from glob import glob
#from dask.diagnostics import ProgressBar
import sys
import argparse
import pathlib  # For handling paths
from pathlib import Path

xrod = xr.open_dataset

# remove spikes
def despike_profile(x, y, sigma=5):
    """
    Remove 6-sigma spikes from a profile.

    Parameters:
        x (array-like): depth points
        y (array-like): Values at time points
        sigma (int, optional): Sigma value for spike removal. Defaults to 5.

    Returns:
        tuple: x, y, and a mask for missing values
    """
    import numpy as np
    dx = np.diff(x) 
    dx[dx == 0] = np.inf # avoid division by zero

    diff1 = np.r_[np.diff(y) / dx, 0]
    diff2 = np.r_[0, np.diff(y) / dx]

    mm = np.nanstd(np.abs(diff1 * diff2)) 

    mk = (diff1 * diff2) < -mm * sigma #spikes has large gradients on both sides of the point with opposite signs
    z= np.where(mk, np.nan, y)
    mk=np.isfinite(z)
    return x[mk], y[mk], z[mk]

def remove_spikes(df, sigma=5,varname='steric'):
    """
    Remove spikes from a time series in a time series in pandas DataFrame.

    Parameters:
        df (pandas.DataFrame): Input DataFrame with columns 'teric', 'surface_time'
        sigma (int, optional): Sigma value for spike removal. Defaults to 5.

    Returns:
        tuple: steric, time, and a mask for missing values
    """
    # Calculate differences in surface_time and steric
    dt = np.diff(df['surface_time'])
    diff1 = np.r_[np.diff(df[varname]) / dt, 0]
    diff2 = np.r_[0, np.diff(df[varname]) / dt]

    # Calculate standard deviation of gradients
    mm = np.nanstd(np.abs(diff1 * diff2))

    # Remove spikes with gradients that are more than 6 times the standard deviation
    mk = (diff1 * diff2) < -mm * sigma

    # Print number of spikes being removed
    print(f"Removing {mk.sum()} spikes from {len(df['steric'])} total points")

    # Return cleaned steric, surface_time, and mask for missing values
    return df[varname], df['surface_time'], np.where(mk, np.nan, df[varname])


def calculate_steric_height_from_profile(depth, density, time,
                    integration_depth=[-500,0],
                    max_spacing=10,
                    time_span_minimum=10 * 60,
                    time_span_maximum=120 * 60,
                    rho0=1027.5,
                    lat=35.5,
                    surface_time_depth=-100,
                    bottom_boundary_tolerance=30,
                    top_boundary_tolerance=10,
                    surface_depth=-8,
                    is_depth=True):

    """
    depth, density are from the Rutgers glider data. No need to clean before passing to this function.
    This routine is profile based, different from the mooring routine which deals with long 1D arrays organized by time.
    """ 
    import utils
    if is_depth:    #make sure depth is negative
        depth=-np.abs(depth)
    else:
        depth=np.abs(depth) #make sure pressure is positive
    
    # clean the data
    msk = np.isfinite(depth) & np.isfinite(density) & (depth >= integration_depth[0]) & (depth <= integration_depth[1])
    
    if msk.sum() < (integration_depth[1]-integration_depth[0])/5:
        print(f"Less than {(integration_depth[1]-integration_depth[0])/5} points in the profile. Return NaN.")
        return 0,np.nan,np.nan,np.nan,np.nan
    
    depth, density,time = depth[msk], density[msk], time[msk]
    if type(time[0]) is np.datetime64:
        #convert to seconds since 2023-01-01
        time = (time - np.datetime64('2023-01-01T00:00:00')) / np.timedelta64(1, 's')
    mk = np.argsort(depth) #sort the data according to depth [-500 -- 0]
    depth, density, time = depth[mk], density[mk], time[mk]
    
    #if the upper {top_boundary_tolerance} meters below the integration depth at the top is not sampled, return NaN. Depth is garrenteed negative at this point.
    if depth[-1]<integration_depth[1]-top_boundary_tolerance: 
        print(f"The shallowest sampling depth {depth[-1]} is below the threshold  {integration_depth[1]-top_boundary_tolerance}. Return NaN.")
        return 0,np.nan,np.nan,np.nan,np.nan # return NaN if the upper 15 meters are not sampled
    # if the deepest sampling is above the integration depth, return NaN
    if depth[0]>integration_depth[0]+bottom_boundary_tolerance: #remember depth is negative
        print(f"The deepest sampling depth {depth[0]} is above the threshold  {integration_depth[0]+bottom_boundary_tolerance}. Return NaN.")
        return 0,np.nan,np.nan,np.nan,np.nan
    
    if np.nanmax(np.diff(depth))>max_spacing:
        print(f"Gaps in the data is larger than the {max_spacing}-meter threshold. Return NaN.")
        return 0,np.nan,np.nan,np.nan,np.nan

    #remove the mean density profile to get density anomaly. This matters significantly for the fixed CTDs that move vertically. 
    density=utils.remove_mean_density(depth,density,is_depth=is_depth,)
    
    surface_time = time[depth>surface_time_depth].mean() #the time when the shallowest depth is sampled
    
    # replace the upper 6 meters with the median density of the upper 8 meters  
    if integration_depth[1]>surface_depth:
        msk = depth > surface_depth 
        if msk.sum()>0:
            density[msk] = np.median(density[msk])
    
    depth,_,density_nospikes=despike_profile(depth,density)
    depth_min,depth_max=depth[0],depth[-1]
    num_points=len(depth)
    if depth[-1]<integration_depth[1]: #if the shallowest depth is above the integration   
        depth=np.r_[depth,integration_depth[1]] #add the upper bound of the integration depth
        density_nospikes=np.r_[density_nospikes,density_nospikes[-1]] #add the density at the shallowest depth to the density array
    if depth[0]>integration_depth[0]: #if the deepest depth is below the integration depth
        depth=np.r_[integration_depth[0],depth] #add the lower bound of the integration depth
        density_nospikes=np.r_[density_nospikes[0],density_nospikes] #add the density at the deepest depth to the density array
    # Calculate the steric height
    steric = -np.trapz(density_nospikes, depth) / rho0 #note the negative sign in front of the integral, corresponding to $steric=-\int_{h_0}^{h_1} \rho/rho0 dz$ 
    
    return int(num_points),int(surface_time),steric*100,depth_min,depth_max

def process_glider(glider_name,
                    integration_depth=[-500,0],
                    max_spacing=10,
                    time_span_minimum=10 * 60,
                    time_span_maximum=120 * 60,
                    rho0=1027.5,
                    surface_time_depth=-100,
                    bottom_boundary_tolerance=30,
                    top_boundary_tolerance=10,
                    surface_depth=-8):

    from glob import glob 
    import utils
    print(f"processing {glider_name}")
    path_input=f"../data/rutgers/{glider_name}_files/*nc"
    files=sorted(glob(path_input))

    dout=np.zeros((1,10))
    ivalid=0
    for i,fn in enumerate(files):
        data=xrod(fn)
        
        depth=data['depth'].values
        density=data['density_adjusted'].values
        lat=data['profile_lat'].values
        lon=data['profile_lon'].values
        time=(data['time'].values - np.datetime64('2023-01-01T00:00:00'))/np.timedelta64(1,'s')
        
        try:
            time_min=np.nanmin(time[np.abs(depth)<-integration_depth[0]])
            time_max=np.nanmax(time[np.abs(depth)<-integration_depth[0]])
        except:
            time_min=np.nan
            time_max=np.nan
            
        nn,surface_time,steric,depth_min,depth_max=calculate_steric_height_from_profile(depth,density,time,lat=lat,
                    integration_depth=integration_depth,
                    max_spacing=max_spacing,
                    time_span_minimum=time_span_minimum,
                    time_span_maximum=time_span_maximum,
                    rho0=rho0,
                    surface_time_depth=surface_time_depth,
                    bottom_boundary_tolerance=bottom_boundary_tolerance,
                    top_boundary_tolerance=top_boundary_tolerance,
                    surface_depth=surface_depth)
        
        if np.isfinite(steric):
            ivalid+=1
        dout=np.r_[dout,np.array([lat,lon,time_min,time_max,surface_time,
                                  depth_min,depth_max,
                                  nn,steric,glider_name]).reshape(1,10)]
        del data
        print(fn, f'steric={steric:.2f}')
    dout=pd.DataFrame(dout[1:,:],columns=['lat','lon',
                                          'time_min','time_max',                                     
                                          'surface_time',
                                          'depth_min',
                                          'depth_max',
                                          'num_points',
                                          'steric',
                                          'Mooring_ID'])
    #drop rows that have NaNs in steric
    #dout=dout.dropna(subset=['steric'])
    print(f'finish processing {glider_name}, total lines {len(dout)} with {ivalid} valid.')
    #dout.to_csv('../data/rutgers/ru32_ru38_steric_heights.csv')
    return dout

def calculate_mooring_steric(time, data, mooring_name, 
                     integration_depth=[-500,0], # depth range for steric height calculation
                     max_spacing=10, # minimum gap threshold for identifying a good profile
                     time_span_minimum=10 * 60, # minimum time span for steric height calculation
                     time_span_maximum=120 * 60, # two hours is a long time for a profile to be taken without anything wrong 
                     rho0=1027.5,
                     surface_time_depth=-100, # the depth range between surface_time_depth and surface used to estimate mean time
                     bottom_boundary_tolerance=30, # maximum depth gap on the edge of profiles at the bottom
                     top_boundary_tolerance=10, # maximum depth gap on the edge of profiles at the top
                     surface_depth=-8, # the depth of the surface layer, in which density will be replaced by the median of the layer. This is to remove spikes and outliers often found in the surface layer.
                     ):
    """
    Calculate Steric Height from a 1D Mooring Profile Data

    This function calculates the steric height from 1D mooring profile data.
    All measurements are flattened to a 1D array in chronological order.

    Parameters:
        time (float): Time in seconds since January 1, 2023
        data (DataFrame): Mooring data table created by density_all_moorings_gliders_level_2.py
        mooring_name (str or int): Mooring identifier (S[1-4] or P[1-7])
        integration_depth (list, optional): Depth range for steric height calculation. Defaults to [-500, 0]
        max_spacing (int, optional): Minimum gap threshold for identifying a good profile. Defaults to 10 meters. Profiles with gaps larger than this value will be ignored.
        time_span_minimum (int, optional): Minimum time span for steric height calculation. Defaults to 25*60
        rho0 (float, optional): Density at sea level. Defaults to 1027.5
        surface_time_depth (int, optional): Depth range between surface_time_depth and surface used to estimate mean time. Defaults to -100
        bottom_boundary_tolerance (int, optional): Maximum depth gap on the edge of profiles at the bottom. Defaults to 30 meters. The missing data will be filled with last valid data, i.e, if integration depth is [-500,0] and the deepest point is -471 meters, the missing data between -471 and -500 meters will be filled with the density at -471 meters. This will introduce a small error in the steric height calculation. 
    Returns:
        DataFrame: Table with steric height and other information.
        ['lat', 'lon', 'time_min', 'time_max', 'surface_time', 'num_points', 'steric', 'Mooring_ID']
    """


    def process_profile(time1,lat,lon,depth1,den1):
        lat1, lon1 = np.nanmean(lat), np.nanmean(lon)
        
        num_points,surface_time,steric,depth_min,depth_max=calculate_steric_height_from_profile(depth1, den1, time1,
                    integration_depth=integration_depth,
                    max_spacing=max_spacing,
                    time_span_minimum=time_span_minimum,
                    time_span_maximum=time_span_maximum,
                    rho0=rho0,
                    lat=lat1,
                    surface_time_depth=surface_time_depth,
                    bottom_boundary_tolerance=bottom_boundary_tolerance,
                    top_boundary_tolerance=top_boundary_tolerance,
                    surface_depth=surface_depth,
                    is_depth=True)
        try:
            time_min = np.nanmin(time1[np.abs(depth1)<-integration_depth[0]])
            time_max = np.nanmax(time1[np.abs(depth1)<-integration_depth[0]])
        except:
            time_min = np.nan
            time_max = np.nan
            
        return [lat1, lon1, time_min, time_max, surface_time, depth_min, depth_max, num_points, steric]
    
    delta = (time - np.datetime64('2023-01-01T00:00:00')) / np.timedelta64(1, 's')
    pids = np.r_[0, np.diff(delta.flatten())]
    pids[pids < 30] = 0
    pids[pids > 30] = 1
    pids = pids.cumsum() #split profiles through the time gaps

    lat = data.sel(variable='lat').values # [0, :], data[1, :], data[2, :], data[3, :]
    lon = data.sel(variable='lon').values # [1, :]
    pressure = data.sel(variable='pressure').values # [2, :]
    density = data.sel(variable='density').values # [3, :]
    
    depth = gsw.z_from_p(pressure, np.nanmean(lat))
    depth = -np.abs(depth)  # make sure depth is negative
    num_profiles = int(pids.max()) + 1 # number of profiles, The "+1" is because pids starts from 0
    dout=np.zeros((1,9))
    for i in tqdm(range(num_profiles)):
        mk = pids == i
        if mk.sum()>10:
            a=process_profile(delta[mk],lat[mk],lon[mk],depth[mk],density[mk])
            dout=np.r_[dout,np.array(a).reshape(1,9)]
    dout=dout[1:,:]
    print(dout.shape)
    df = pd.DataFrame(dout, columns=['lat', 'lon', 'time_min', 
                                     'time_max', 'surface_time', 
                                     'depth_min', 'depth_max',
                                     'num_points', 'steric'])
    df['Mooring_ID'] = mooring_name[:2]
    a,tt,b=remove_spikes(df) 
    df['steric']=b
    #df=df.dropna(subset=['steric'])
    return df

def parse_arguments():
    """Parse command-line arguments for all constants."""
    parser = argparse.ArgumentParser(description="Steric height calculation parameters.")
    
    # Depth range arguments
    parser.add_argument("--bottom_depth", type=int, default=-500, help="Bottom depth (negative value in meters).")
    parser.add_argument("--top_depth", type=int, default=0, help="Top depth (negative value in meters).")
    
    # Steric height calculation parameters
    parser.add_argument("--max_spacing", type=int, default=10, help="Minimum gap threshold for identifying a good profile (km).")
    parser.add_argument("--time_span_minimum", type=int, default=10*60, help="Minimum time span for steric height calculation (seconds).")
    parser.add_argument("--time_span_maximum", type=int, default=120*60, help="Maximum time span for steric height calculation (seconds).")
    parser.add_argument("--rho0", type=float, default=1027.5, help="Density at sea level (kg/m^3).")
    parser.add_argument("--surface_time_depth", type=int, default=-100, help="Depth range for surface time calculation (m).")
    parser.add_argument("--bottom_boundary_tolerance", type=int, default=80, help="Maximum depth gap at the bottom (m).")
    parser.add_argument("--top_boundary_tolerance", type=int, default=10, help="Maximum depth gap at the top (m).")
    parser.add_argument("--surface_depth", type=int, default=-8, help="Depth of the surface well-mixed layer (m).")
    
    return parser.parse_args()

def validate_depth_range(integration_depth):
    """Ensure the depth range is within a reasonable range."""
    bottom_depth, top_depth = integration_depth
    if bottom_depth >= 0 or top_depth > 0:
        print("Error: Depth values must be negative, representing depths below the surface.")
        sys.exit(1)
    if abs(bottom_depth - top_depth) < 10:  # Depth difference should be at least 10m
        print("Error: The depth difference is too small.")
        sys.exit(1)

def print_parameters(args):
    """Print the configured parameters for the calculation."""
    print(f"Configured with depth range: {args.bottom_depth}m to {args.top_depth}m")
    print(f"Minimum gap threshold for identifying a good profile: {args.max_spacing} km")
    print(f"Minimum time span: {args.time_span_minimum / 60} minutes")
    print(f"Maximum time span: {args.time_span_maximum / 60} minutes")
    print(f"Density at sea level: {args.rho0} kg/m^3")
    print(f"Surface time depth: {args.surface_time_depth} m")
    print(f"Bottom boundary tolerance: {args.bottom_boundary_tolerance} m")
    print(f"Top boundary tolerance: {args.top_boundary_tolerance} m")
    print(f"Surface well-mixed layer depth: {args.surface_depth} m")

def save_to_csv(dataframe, filename):
    """
    Save the results to a CSV file.
    """
    
    dataframe.to_csv(filename, index=False, float_format='%.4f')
    print(f"Results saved to {filename}")
    

def process_mooring_data(mooring_ids, ds, integration_depth, 
                         max_spacing, time_span_minimum, 
                         time_span_maximum, rho0, surface_time_depth, 
                         bottom_boundary_tolerance, surface_depth):
    """
    Process mooring data to calculate steric height.
    """
    results = []  # Initialize results as a list to hold individual DataFrame results
    
    for var_name in mooring_ids:
        print(f"{'='*20}\nProcessing {var_name}\n{'='*18}")
        time = ds[f'time_{var_name}'].values
        data = ds[var_name]
        
        # Calculate steric for each mooring
        aa = calculate_mooring_steric(time, data, var_name, 
                                      integration_depth=integration_depth, 
                                      max_spacing=max_spacing, 
                                      time_span_minimum=time_span_minimum, 
                                      time_span_maximum=time_span_maximum, 
                                      rho0=rho0, 
                                      surface_time_depth=surface_time_depth, 
                                      bottom_boundary_tolerance=bottom_boundary_tolerance, 
                                      surface_depth=surface_depth)
        
        # Append results to list
        results.append(aa)
    
    # Concatenate all results into a single DataFrame
    return pd.concat(results, ignore_index=True)

def main():
    """Main function to execute the steric height calculation."""
    # Step 1: Parse and validate command-line arguments
    args = parse_arguments()
    integration_depth = [args.bottom_depth, args.top_depth]
    validate_depth_range(integration_depth)

    # Step 2: Print configured parameters
    print_parameters(args)
    
    # Flags to control whether to calculate for gliders or moorings
    calculate_glider = True 
    calculate_mooring = True

    # Step 3: Process gliders if the flag is set
    if calculate_glider:
        # Process both gliders
        glider_ids = ['ru32', 'ru38']
        glider_dataframes = []
        for glider_id in glider_ids:
            df=process_glider(glider_id,
                    integration_depth=integration_depth,
                    max_spacing=args.max_spacing,
                    time_span_minimum=args.time_span_minimum,
                    time_span_maximum=args.time_span_maximum,
                    rho0=args.rho0,
                    surface_time_depth=args.surface_time_depth,
                    bottom_boundary_tolerance=args.bottom_boundary_tolerance,
                    surface_depth=args.surface_depth)
            glider_dataframes.append(df)
        
        # Concatenate all glider data and save it to CSV
        glider_df = pd.concat(glider_dataframes)
        #glider_df = glider_df.dropna(subset=['steric'])  # Drop rows where steric height is NaN
        glider_filename = f'../data/rutgers/ru32_ru38_steric_heights_depth_{-integration_depth[0]}.csv'
        save_to_csv(glider_df, glider_filename)

    # Step 4: Process moorings if the flag is set
    if calculate_mooring:
        # List of mooring IDs to process
        mooring_ids = [f'S{i}_P' for i in range(1, 5)] + [f'P{i}_P' for i in range(1, 8)] 
        
        # Load the dataset
        ds_path = Path('../data/mooring.data/density_all_moorings_level-2_removed_outliers.nc')
        if not ds_path.exists():
            print(f"Error: Dataset {ds_path} not found.")
            sys.exit(1)

        ds = xr.open_dataset(ds_path)
        
        # Process mooring data
        mooring_df = process_mooring_data(mooring_ids, ds, integration_depth, args.max_spacing, 
                                          args.time_span_minimum, args.time_span_maximum, args.rho0, 
                                          args.surface_time_depth, args.bottom_boundary_tolerance, 
                                          args.surface_depth)
        
        # Save the results to CSV
        mooring_filename = f'../data/mooring.data/all_mooring_steric_heights_depth_{-args.bottom_depth:3d}.csv'
        save_to_csv(mooring_df, mooring_filename)

if __name__ == '__main__':
    main()