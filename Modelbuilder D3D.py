# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 14:42:08 2023

@author: jvoor
"""

# import packages
import os
import matplotlib.pyplot as plt
plt.close('all')
import dfm_tools as dfmt
import hydrolib.core.dflowfm as hcdfm
import xarray as xr
import pandas as pd
import numpy as np
import contextily as ctx

# user input
model_name = 'Saba'
dir_output = os.path.abspath(f'./{model_name}_model')
# path_style = 'windows' # windows / unix
overwrite = False # used for downloading of forcing data. Always set to True when changing the domain
is_geographic = True # spherical (True) or cartesian (False) coordinates
crs = 'EPSG:4326' #EPSG:3857 # coordinate reference system

# domain and resolution
if model_name=='Saba':
    lon_min, lon_max, lat_min, lat_max = -63.27, -63.21, 17.60, 17.66 #for EPSG:4326
    #lon_min, lon_max, lat_min, lat_max = -7043117, -7035516, 1990500, 1998332 #for EPSG:3857
dxy = 0.001#fine for EPSG:4326 0.001

#dates as understood by pandas.period_range(). ERA5 has freq='M' (month) and CMEMS has freq='D' (day)
date_min = '2023-03-01'
date_max = '2023-03-30'
ref_date = '2023-01-01'

# make directories and list all files
os.makedirs(dir_output, exist_ok=True)
dir_output_data = os.path.join(dir_output, 'data')
os.makedirs(dir_output_data, exist_ok=True)

#####################################################################################################################################

# generate spherical regular grid
mk_object = dfmt.make_basegrid(lon_min, lon_max, lat_min, lat_max, dx=dxy, dy=dxy, is_geographic=is_geographic)

# generate plifile from grid extent and coastlines
bnd_gdf = dfmt.generate_bndpli_cutland(mk=mk_object, res='h', buffer=0.01)
# bnd_gdf['name'] = f'{model_name}_bnd_{bnd_gdf.index}'  # Appending index to ensure uniqueness
# bnd_gdf['name'] = f'{model_name}_bnd1'
bnd_gdf_interp = dfmt.interpolate_bndpli(bnd_gdf, res=0.003)


poly_file = os.path.join(dir_output, f'{model_name}.pli')
pli_polyfile = dfmt.geodataframe_to_PolyFile(bnd_gdf_interp)
pli_polyfile.save(poly_file)

# plot basegrid and polyline
fig, ax = plt.subplots()
mk_object.mesh2d_get().plot_edges(ax,zorder=1)
bnd_gdf_interp.plot(ax=ax, edgecolor='r')
ctx.add_basemap(ax=ax, crs=crs, attribution=False) #=topo background layer
# dfmt.plot_coastlines(ax=landboundaryshp, crs=crs) #plotting landboundary

landboundaryshp = r'C:\Users\jvoor\Delft University of Technology\MSc thesis - General\06 - Modelling\Delft3D FM\03 Pre-post processing D3D\D3Dfm_preprocessing\Saba_model\Saba_ldb.shp'
dfmt.meshkernel_delete_withshp(mk=mk_object, coastlines_shp=landboundaryshp)

# ax.plot(x_coordinates, y_coordinates, color='black', linestyle = '-', linewidth=1)  # Plotting coastlines from ldb
# dfmt.plot_coastlines(ax=h, crs=crs)  ### plotting coastline (very course) from dfmt tools
"loading own landboundary file"
ldb = np.loadtxt(r'C:\Users\jvoor\Delft University of Technology\MSc thesis - General\06 - Modelling\Delft3D FM\03 Pre-post processing D3D\D3Dfm_preprocessing\output_coordinates_4326.txt', delimiter = ',')
x_coordinates = ldb[:, 0]
y_coordinates = ldb[:, 1]
plt.plot(x_coordinates, y_coordinates, color='black', linestyle = '-', linewidth=1)  # Plotting coastlines from ldb


"still change rest of the code to use own landboundary file instead of coastlines file"
# connect to a coarse version of the GEBCO_2022 dataset on OPeNDAP
# alternatively download your own cutout from https://download.gebco.net (use a buffer of e.g. 1 degree)
file_gebco = f'https://opendap.deltares.nl/thredds/dodsC/opendap/deltares/Delft3D/netcdf_example_files/GEBCO_2022/GEBCO_2022_coarsefac08.nc'
data_bathy_sel = xr.open_dataset(file_gebco).sel(lon=slice(lon_min-1,lon_max+1),lat=slice(lat_min-1,lat_max+1))

# refine grid
min_edge_size = 100 # in meters
dfmt.refine_basegrid(mk=mk_object, data_bathy_sel=data_bathy_sel, min_edge_size=min_edge_size)

# plot
# fig, ax = plt.subplots()
# mk_object.mesh2d_get().plot_edges(ax,zorder=1)
# ctx.add_basemap(ax=ax, crs=crs, attribution=False)
# #dfmt.plot_coastlines(ax=ax, crs=crs)

# # remove land with GSHHS coastlines
# dfmt.meshkernel_delete_withcoastlines(mk=mk_object, res='h')

# # plot
# fig, ax = plt.subplots()
# mk_object.mesh2d_get().plot_edges(ax,zorder=1)
# ctx.add_basemap(ax=ax, crs=crs, attribution=False)
# dfmt.plot_coastlines(ax=ax, crs=crs)

# convert to xugrid
xu_grid_uds = dfmt.meshkernel_to_UgridDataset(mk=mk_object, crs=crs)

# interpolate bathymetry onto the grid
data_bathy_interp = data_bathy_sel.interp(lon=xu_grid_uds.obj.mesh2d_node_x, lat=xu_grid_uds.obj.mesh2d_node_y).reset_coords(['lat','lon'])
xu_grid_uds['mesh2d_node_z'] = data_bathy_interp.elevation.clip(max=10)

# plot bathymetry and grid
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.mesh2d_node_z.ugrid.plot(ax=ax,center=False)
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='white',alpha=0.2)
ctx.add_basemap(ax=ax, crs=crs, attribution=False)
dfmt.plot_coastlines(ax=ax, crs=crs)

# write xugrid grid to netcdf
netfile = os.path.join(dir_output, f'{model_name}_net.nc')
xu_grid_uds.ugrid.to_netcdf(netfile)

# ###############################################################################################################
# #FROM HERE EXTERNAL B.C.

# generate new format external forcings file (.ext): initial and open boundary condition
ext_file_new = os.path.join(dir_output, f'{model_name}_new.ext')
ext_new = hcdfm.ExtModel()

# interpolate tidal components to boundary conditions file (.bc)
tidemodel = 'tpxo80_opendap' # tidemodel: FES2014, FES2012, EOT20, GTSMv4.1, GTSMv4.1_opendap
ForcingModel_object = dfmt.interpolate_tide_to_bc(tidemodel=tidemodel, file_pli=poly_file, component_list=None)
file_bc_out = os.path.join(dir_output,f'tide_{model_name}_{tidemodel}.bc')
ForcingModel_object.save(filepath=file_bc_out)
boundary_object = hcdfm.Boundary(quantity='waterlevelbnd',
                                  locationfile=poly_file,
                                  forcingfile=ForcingModel_object)
ext_new.boundary.append(boundary_object)

# CMEMS - download spatial fields of salinity, temperature, currents and sea surface height
# you can also add WAQ variables like 'no3' and 'phyc'
# check dfmt.get_conversion_dict() for an overview of parameter/quantity names
dir_output_data_cmems = os.path.join(dir_output_data, 'cmems')
os.makedirs(dir_output_data_cmems, exist_ok=True)
for varkey in ['so','thetao','uo','vo','zos']: #so=salinity, thetao=temperature, uo=advectionvelocity, vo=advectionvelocity, zos = waterlvlboundary
    dfmt.download_CMEMS(varkey=varkey,
                        longitude_min=lon_min, longitude_max=lon_max, latitude_min=lat_min, latitude_max=lat_max,
                        date_min=date_min, date_max=date_max,
                        dir_output=dir_output_data_cmems, file_prefix='cmems_', overwrite=overwrite)

# CMEMS - boundary conditions file (.bc) and add to ext_bnd
# you can also add WAQ variables like 'tracerbndNO3' and 'tracerbndPON1'
# check dfmt.get_conversion_dict() for an overview of parameter/quantity names
list_quantities = ['waterlevelbnd','salinitybnd','temperaturebnd','uxuyadvectionvelocitybnd']
dir_pattern = os.path.join(dir_output_data_cmems,'cmems_{ncvarname}_*.nc')
ext_new = dfmt.cmems_nc_to_bc(ext_bnd=ext_new,
                              refdate_str=f'minutes since {ref_date} 00:00:00 +00:00',
                              dir_output=dir_output,
                              list_quantities=list_quantities,
                              tstart=date_min,
                              tstop=date_max, 
                              file_pli=poly_file,
                              dir_pattern=dir_pattern)
dfmt.get
#save new ext file
ext_new.save(filepath=ext_file_new) # ,path_style=path_style)

# plot downloaded CMEMS data
file_cmems = os.path.join(dir_output_data,'cmems','*.nc')
ds_cmems = xr.open_mfdataset(file_cmems)
ds_cmems

# plot
fig, (ax1,ax2) = plt.subplots(1, 2, figsize=(10,5))
ds_cmems.so.isel(time=0, depth=0).plot(ax=ax1)
dfmt.plot_coastlines(ax=ax1, crs=crs)
ds_cmems.thetao.isel(time=0, depth=0).plot(ax=ax2)
dfmt.plot_coastlines(ax=ax2, crs=crs)
fig.tight_layout()

# plot interpolated CMEMS data (boundary conditions in .bc)
file_bc_sal = os.path.join(dir_output, f'salinitybnd_{model_name}_CMEMS.bc')
bc_obj_sal = hcdfm.ForcingModel(file_bc_sal)
forcing_xr_sal = dfmt.forcinglike_to_Dataset(bc_obj_sal.forcing[0], convertnan=True)

file_bc_uxuy = os.path.join(dir_output,f'uxuyadvectionvelocitybnd_{model_name}_CMEMS.bc')
bc_obj_uxuy = hcdfm.ForcingModel(file_bc_uxuy)
forcing_xr_uxuy = dfmt.forcinglike_to_Dataset(bc_obj_uxuy.forcing[0], convertnan=True)

# plot
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(10,8))
forcing_xr_sal['salinitybnd'].T.plot(ax=ax1)
forcing_xr_uxuy['ux'].T.plot(ax=ax2)
forcing_xr_uxuy['uy'].T.plot(ax=ax3)
ax1.set_ylim(xu_grid_uds.mesh2d_node_z.min(), None)

# generate old format external forcings file (.ext): spatial data
ext_file_old = os.path.join(dir_output, f'{model_name}_old.ext')
ext_old = hcdfm.ExtOldModel()

# CMEMS - initial conditions
ext_old = dfmt.cmems_nc_to_ini(ext_old=ext_old,
                                dir_output=dir_output,
                                list_quantities=list_quantities,
                                tstart=date_min,
                                dir_pattern=dir_pattern)

# ERA5 - download spatial fields of air pressure, wind speeds and Charnock coefficient
dir_output_data_era5 = os.path.join(dir_output_data, 'ERA5')
os.makedirs(dir_output_data_era5, exist_ok=True)
    
varlist_list = [['msl','u10n','v10n','chnk']]

for varlist in varlist_list:
    for varkey in varlist:
        dfmt.download_ERA5(varkey, 
                            longitude_min=lon_min, longitude_max=lon_max, latitude_min=lat_min, latitude_max=lat_max,
                            date_min=date_min, date_max=date_max,
                            dir_output=dir_output_data_era5, overwrite=overwrite)

# ERA5 meteo - convert to netCDF for usage in Delft3D FM
ext_old = dfmt.preprocess_merge_meteofiles_era5(ext_old=ext_old,
                                                varkey_list=varlist_list,
                                                dir_data=dir_output_data_era5,
                                                dir_output=dir_output,
                                                time_slice=slice(date_min, date_max))

ext_old.save(filepath=ext_file_old) # , path_style=path_style)

# plot converted ERA5 data
file_era5 = os.path.join(dir_output,'data','ERA5','*.nc')
ds_era5 = xr.open_mfdataset(file_era5)
ds_era5

# plot
fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,5))
ds_era5.u10n.isel(time=0).plot(ax=ax1)
dfmt.plot_coastlines(ax=ax1, crs=crs)
ds_era5.v10n.isel(time=0).plot(ax=ax2)
dfmt.plot_coastlines(ax=ax2, crs=crs)
fig.tight_layout()


###########################################################################################
'Now add observation points'

"adjust script to add own list of observation points. Either directly or from a textfile"

# generate obspoints on all grid faces
xpts = xu_grid_uds.grid.face_x
ypts = xu_grid_uds.grid.face_y
npts = [f'x{x:.2f}_y{y:.2f}'.replace('.','p') for x,y in zip(xpts,ypts)]
obs_pd = pd.DataFrame(dict(x=xpts,y=ypts,name=npts))

# subselect n arbitary obspoints and plot
n = 26
ipts = np.random.randint(0, len(obs_pd), n)
obs_pd = obs_pd.iloc[ipts]
print(obs_pd)
fig, ax = plt.subplots(figsize=(8,4))
xu_grid_uds.grid.plot(ax=ax,linewidth=0.5,color='k',alpha=0.2)
# ax.plot(obs_pd['x'],obs_pd['y'],'rx')
# dfmt.plot_coastlines(ax=ax, crs=crs)

# save obsfile
file_obs = os.path.join(dir_output, f'{model_name}_obs.xyn')
obs_pd.to_csv(file_obs, sep=' ', header=False, index=False, float_format='%.6f')

###########################################################################################

'mdu-file generation'

# initialize mdu file and update settings
mdu_file = os.path.join(dir_output, f'{model_name}.mdu')
mdu = hcdfm.FMModel()

# add the grid (_net.nc, network file)
mdu.geometry.netfile = netfile

# support for initial sal/tem fields via iniwithnudge, this requires 3D model
# mdu.geometry.kmx = 5
# mdu.physics.iniwithnudge = 2

# add the external forcing files (.ext)
mdu.external_forcing.extforcefile = ext_file_old
mdu.external_forcing.extforcefilenew = ext_new

# update time settings
mdu.time.refdate = pd.Timestamp(ref_date).strftime('%Y%m%d')
mdu.time.tunit = 'S'
mdu.time.dtmax = 30
mdu.time.startdatetime = pd.Timestamp(date_min).strftime('%Y%m%d%H%M%S')
mdu.time.stopdatetime = pd.Timestamp(date_max).strftime('%Y%m%d%H%M%S')
mdu.time.autotimestep = 3

# update output settings
#mdu.output.obsfile = file_obs
mdu.output.hisinterval = [60]
mdu.output.mapinterval = [1800]#[86400]
mdu.output.rstinterval = [0]
mdu.output.statsinterval = [3600]

# save .mdu file
mdu.save(mdu_file) # ,path_style=path_style)

# visualize the model tree
mdu.show_tree()

################
nproc = 1 # number of processes
dimrset_folder = r"C:\Program Files (x86)\Deltares\Delft3D Flexible Mesh Suite HMWQ (2021.03)\plugins\DeltaShell.Dimr\kernels" #alternatively r"p:\d-hydro\dimrset\weekly\2.25.17.78708"
dfmt.create_model_exec_files(file_mdu=mdu_file, nproc=nproc, dimrset_folder=dimrset_folder)

###################3

# list all files
print(f"your model is now ready in {dir_output}:")
os.listdir(dir_output)



