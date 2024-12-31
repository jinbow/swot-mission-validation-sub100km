import xarray as xr 
from cartopy import crs as ccrs
import cartopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Load SSH data from SWOT and NeurOST
#fn="../data/NeurOST_SSH-SST_20230403_20240722.nc"
#dd=xr.open_dataset(fn).sel(longitude=slice(360-128,360-120)).sel(latitude=slice(30,40))
dd=xr.open_dataset('../data/NeurOST_SSH-SST_20230403_20240722_sub_lat-30-40.nc')
print('load NeurOST data from ../data/NeurOST_SSH-SST_20230403_20240722_sub_lat-30-40.nc')
sla_neurost=dd.sla[...,0]
sla_neurost-=np.nanmean(sla_neurost)
lon_neurost=dd.longitude
lat_neurost=dd.latitude

timec='2023-04-03'
fn='../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_013_sub_lat-30-40.nc'
dd=xr.open_dataset(fn).sel(cycle=479)
print('load swot data from ',fn)

#.sel(timec=timec)
karin=(dd['ssha_karin']+dd['height_cor_xover']+dd['internal_tide_hret'])
karin-=np.nanmean(karin)
qual=dd['ssha_karin_qual']

karin=np.where(qual==0,karin,np.nan)
lon_karin=dd.longitude
lat_karin=dd.latitude
swot=[[lon_karin, lat_karin, karin]]

fn='../data/SWOT_L2_LR_SSH_2.0_combined_calval_orbit_pass_026_sub_lat-30-40.nc'
dd=xr.open_dataset(fn).sel(cycle=479)
print('load swot data from ',fn)

karin=(dd['ssha_karin']+dd['height_cor_xover']+dd['internal_tide_hret'])
karin-=np.nanmean(karin)
qual=dd['ssha_karin_qual']
karin=np.where(qual==0,karin,np.nan)
lon_karin=dd.longitude
lat_karin=dd.latitude

swot.append([lon_karin, lat_karin, karin])

from cartopy import crs as ccrs
import cartopy
import utils
#set a basemap using map projection

fig=plt.figure(figsize=(7,7))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()
#set continent color to grey
ax.add_feature(cartopy.feature.LAND, zorder=2
               , edgecolor='black', facecolor='#f2f2f2')
ax.set_extent([-128,-120,33,39])
# remove the labels on top and right side
ax.set_xticks(np.arange(-128,-119,1), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(33,40,1), crs=ccrs.PlateCarree())
#lon_formatter = cartopy.mpl.ticker.LongitudeFormatter()
#lat_formatter = cartopy.mpl.ticker.LatitudeFormatter()
#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

#add gridlines
#ax.gridlines(draw_labels=True)
vmin,vmax=-0.15,0.15
cs=ax.contourf(lon_neurost,lat_neurost,sla_neurost, levels=np.linspace(vmin,vmax,21), extend='both', cmap='Spectral_r',vmin=vmin,vmax=vmax, 
              transform=ccrs.PlateCarree(),z_order=1)
cs=ax.pcolormesh(swot[0][0],swot[0][1],swot[0][2],cmap='Spectral_r',vmin=vmin,vmax=vmax,
                 transform=ccrs.PlateCarree(),zorder=2)
cs=ax.pcolormesh(swot[1][0],swot[1][1],swot[1][2],cmap='Spectral_r',vmin=vmin,vmax=vmax,
                 transform=ccrs.PlateCarree(),zorder=2)

plt.colorbar(cs,orientation='horizontal',label='Sea Surface Height Anomaly (meter)',shrink=0.99,extend='both',pad=0.1,aspect=50)

for i in [0,1]:
    ax.plot(swot[i][0][:,[2,30,36,-4]],swot[i][1][:,[2,30,36,-4]],color='gray',marker='.',markersize=1,
            transform=ccrs.PlateCarree(),zorder=10,alpha=0.5)

## load mooring positions
fn='../data/mooring.data/mooring_positions.csv'
print(f'load mooring positions from {fn}')
positions=pd.read_csv(fn)
#select mooring name has letter P
p=positions[positions['name'].str.contains(f'P')]
plt.scatter(p['lon'],p['lat'],color='red',marker='.',s=0.2,transform=ccrs.PlateCarree(),zorder=10)
s=positions[positions['name'].str.contains(f'S')]
plt.scatter(s['lon'],s['lat'],color='blue',marker='.',s=0.2,transform=ccrs.PlateCarree(),zorder=10)

ax.set_extent([-127,-120,33,39])

plt.savefig('../figures/figure1_locations.png',dpi=600,bbox_inches='tight')
print('save figure to ../figures/figure1_locations.png')

for i in range(1,8):
    p=positions[positions['name'].str.contains(f'P{i:1d}')]
    plt.scatter(p['lon'],p['lat'],color='red',marker='.',s=0.2,transform=ccrs.PlateCarree(),zorder=10)
    mean_lon=np.mean(p['lon'])
    mean_lat=np.mean(p['lat'])
    #plt.text(mean_lon+0.05,mean_lat-0.02,f'P{i:1d}', fontsize=15, color='r',transform=ccrs.PlateCarree(),zorder=10)

for i in range(1,5):
    s=positions[positions['name'].str.contains(f'S{i:1d}')]
    plt.scatter(s['lon'],s['lat'],color='blue',marker='.',s=0.2,transform=ccrs.PlateCarree(),zorder=10)
    mean_lon=np.mean(s['lon'])
    mean_lat=np.mean(s['lat'])
    #plt.text(mean_lon-0.1,mean_lat-0.03,f'S{i:1d}', fontsize=15, color='b',transform=ccrs.PlateCarree(),zorder=10)

utils.plot_glider_locs(plt.gca())

ax.set_extent([-125.5,-124.5,35,36.3])

plt.savefig('../figures/figure1_locations_zoom.png',dpi=300,bbox_inches='tight')
print('save figure to ../figures/figure1_locations_zoom.png')
#plt.show()