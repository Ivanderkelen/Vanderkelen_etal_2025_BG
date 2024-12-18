"""
Python module with all utils and functions for analysing CESM output on piz daint
"""

# --------------------------------------------------------------------
# Import modules 
# -------------------------------------------------------------------
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import regionmask
import cartopy.crs as ccrs
import warnings
import numpy as np
from ctsm_py.utils import *

# --------------------------------------------------------------------
# Settings - define variables necessary for functions
# --------------------------------------------------------------------


# set directories
outdir = '/scratch/snx3000/ivanderk/'
outdir = '/project/s1207/ivanderk/scratch/'


# Define directory where processing is done -- subject to change
procdir =  outdir + 'processing_4p1000/' 

# go to processing directory 
os.chdir(procdir)

# ignore all runtime warnings
warnings.filterwarnings('ignore')


# case settings
case_ctl = 'I2000Clm51Sp.hcru_hcru_mt13.CTL'

block  = 'lnd' 
stream = 'h0' 

start_year = 2004

# ---------------------------------------------------------------------
# 1. Functions to open datasets
# ---------------------------------------------------------------------


# get list of files with their full paths 
def listdir_fullpath(d, stream=stream):
    return [os.path.join(d, f) for f in os.listdir(d) if stream in f]

# get a list of files to open ds from
def get_filelist(case, stream=stream, block=block):

    # Define directory where timeseries data is stored
    tseriesdir = outdir + 'archive/' + case + '/' + block + '/hist/'
    
    # define filename
    filelist = listdir_fullpath(tseriesdir, stream=stream)
    
    return filelist

# import dataset for specific case variables

def import_case(
    case,
    myVars=None,
    myVegtypes=None,
    timeSlice=None,
    myVars_missing_ok=[],
    only_active_patches=False,
    rename_lsmlatlon=False,
    chunks=None,
    alter_lons=True,
    stream=stream
):

    # get file list
    filelist = get_filelist(case, stream=stream, block=block)
    
    
    # import ds (using function from ctsm_py)
    ds = import_ds(
        filelist,
        myVars=myVars,
        myVegtypes=None,
        timeSlice=timeSlice,
        myVars_missing_ok=[],
        only_active_patches=False,
        rename_lsmlatlon=False,
        chunks=None,
    )

    if alter_lons:
        # alter lons to span -180 to 180
        ds = lon_pm2idl(ds, fail_silently=False)

    return ds


# regrid column level output from vector to grid
# Adjusted version of the function from the ctsm toolbox to also work with columns output (for cases run with all pfts on single columns)

def grid_one_variable_col(this_ds, thisVar, fillValue=None, **kwargs):
    
    # Get this Dataset's values for selection(s), if provided
    this_ds = xr_flexsel(this_ds, \
        **kwargs)
    
    # Get DataArrays needed for gridding
    thisvar_da = get_thisVar_da(thisVar, this_ds)
    vt_da = None
    
    if "patch" in thisvar_da.dims:
        spatial_unit = "patch"
        xy_1d_prefix = "patches"
        if "patches1d_itype_veg" in this_ds:
            vt_da = get_thisVar_da("patches1d_itype_veg", this_ds)
    elif "column" in thisvar_da.dims:
        spatial_unit = "column"
        xy_1d_prefix = "patches"
        if "patches1d_itype_veg" in this_ds:
            vt_da = get_thisVar_da("patches1d_itype_veg", this_ds)
    elif "gridcell" in thisvar_da.dims:
        spatial_unit = "gridcell"
        xy_1d_prefix = "grid"
    else:
        raise RuntimeError(f"What variables to use for _ixy and _jxy of variable with dims {thisvar_da.dims}?")
    ixy_da = get_thisVar_da(xy_1d_prefix + "1d_ixy", this_ds)
    jxy_da = get_thisVar_da(xy_1d_prefix + "1d_jxy", this_ds)
    
    
    if not fillValue and "_FillValue" in thisvar_da.attrs:
        fillValue = thisvar_da.attrs["_FillValue"]
    
    # Renumber vt_da to work as indices on new ivt dimension, if needed.
    ### Ensures that the unique set of vt_da values begins with 1 and
    ### contains no missing steps.
    if "ivt" in this_ds and vt_da is not None:
        vt_da.values = np.array([np.where(this_ds.ivt.values == x)[0][0] for x in vt_da.values])
    
    # Get new dimension list
    new_dims = list(thisvar_da.dims)
    ### Remove "[spatial_unit]".
    if spatial_unit in new_dims:
        new_dims.remove(spatial_unit)
    #  Add "ivt_str" (vegetation type, as string). This needs to go at the end, to avoid a possible situation where you wind up with multiple Ellipsis members of fill_indices.
    if "ivt" in this_ds and (spatial_unit=="patch" or spatial_unit=="column"):
        new_dims.append("ivt_str")
    ### Add lat and lon to end of list
    new_dims = new_dims + ["lat", "lon"]

    # Set up empty array
    n_list = []
    for dim in new_dims:
        if dim == "ivt_str":
            n = this_ds.sizes["ivt"]
        elif dim in thisvar_da.coords:
            n = thisvar_da.sizes[dim]
        else:
            n = this_ds.sizes[dim]
        n_list = n_list + [n]

    thisvar_gridded = np.empty(n_list)
    if fillValue:
        thisvar_gridded[:] = fillValue
    else:
        thisvar_gridded[:] = np.nan

    # Fill with this variable
    fill_indices = []
    for dim in new_dims:
        if dim == "lat":
            fill_indices.append(jxy_da.values.astype(int) - 1)
        elif dim == "lon":
            fill_indices.append(ixy_da.values.astype(int) - 1)
        elif dim == "ivt_str":
            fill_indices.append(vt_da)
        elif not fill_indices:
        # I.e., if fill_indices is empty. Could also do "elif len(fill_indices)==0".
            fill_indices.append(Ellipsis)
    try:
        thisvar_gridded[tuple(fill_indices[:len(fill_indices)])] = thisvar_da.values
    except:
        thisvar_gridded[tuple(fill_indices[:len(fill_indices)])] = thisvar_da.values.transpose()
    if not np.any(np.bitwise_not(np.isnan(thisvar_gridded))):
        if np.all(np.isnan(thisvar_da.values)):
            print('Warning: This DataArray (and thus map) is all NaN')
        else:
            raise RuntimeError("thisvar_gridded was not filled!")
    
    # Assign coordinates, attributes and name
    thisvar_gridded = xr.DataArray(thisvar_gridded, \
        dims=tuple(new_dims),
        attrs=thisvar_da.attrs)
    for dim in new_dims:
        if dim == "ivt_str":
            values = this_ds.vegtype_str.values
        elif dim in thisvar_da.coords:
            values = thisvar_da[dim]
        else:
            values = this_ds[dim].values
        thisvar_gridded = thisvar_gridded.assign_coords({dim: values})
    thisvar_gridded.name = thisVar
    
    # Add FillValue attribute
    if fillValue:
        thisvar_gridded.attrs["_FillValue"] = fillValue

    return thisvar_gridded


# load auxiliary control case -- not for direct analysis
def load_case_ctl(variables): 
    
    # case settings
    case_ctl = 'IHistClm51Sp.hcru_hcru_mt13.4p1000.spunup'
    case   = case_ctl
    block  = 'lnd' 
    stream = 'h0' 

    # discard 2004 until 2008 as spin up years
    start_year, end_year = 2001, 2001 
    time_slice = slice(str(start_year)+"-01-01", str(end_year)+"-12-01")
    
    # import the control case
    ds_ctl = import_case(case_ctl, myVars=variables, timeSlice=time_slice)
    
    return ds_ctl




# get masks for sandy and clay soils based on given thresholds in pct
def get_texture_masks(pct_sand_treshold_max,pct_sand_treshold_min,  pct_clay_threshold, ds_input, levsoi_thickness):


    # sand and clay datasets 
    da_pct_sand = ds_input['PCT_SAND'].assign_coords({"nlevsoi": (levsoi_thickness[:10])})
    da_pct_clay = ds_input['PCT_CLAY'].assign_coords({"nlevsoi": (levsoi_thickness[:10])})

    da_pct_sand_2d = da_pct_sand.weighted(da_pct_sand.nlevsoi).mean('nlevsoi')
    da_pct_clay_2d = da_pct_clay.weighted(da_pct_clay.nlevsoi).mean('nlevsoi')

    # weight masks over soil profile
    da_sand_mask = (da_pct_sand_2d > pct_sand_treshold_min) * (da_pct_sand_2d < pct_sand_treshold_max)
    da_clay_mask = da_pct_clay_2d > pct_clay_threshold
    
    return da_sand_mask, da_clay_mask

# load gridded pft variables -- from postprocessed dirgrid_multiple_pfts_multiple_vars
def load_da_gridded_pft(variable, pft_int, case, flag_daily=False): 
    
    if flag_daily: 
        text_daily = '_daily'
    else: 
        text_daily = ''
        
    da = xr.open_dataset(procdir+'/postprocessing/pft_gridding/'+variable+'_'+str(pft_int)+'.'+case+text_daily+'.nc', engine='netcdf4')[variable+'_'+str(pft_int)]
    return da

# ---------------------------------------------------------------------
# 2. Functions to perform calculations
# ---------------------------------------------------------------------



# get manually  soil tickness and depths as defined in CLM
def get_soildepths(): 
    
    # tickness at soil layer - delta zi - from CLM user's guide Table 2.3 
    levsoi_thickness = [0.02,0.04,0.6,0.8,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.54,0.64,0.74,0.84,0.94,1.04,1.14]

    # depth at soil layer zh,i  - from CLM user's guide Table 2.3 
    levsoi_depth     = [0.02,0.06,0.12,0.2,0.32,0.48,0.68,0.92,1.2,1.52,1.88,2.28,2.72,3.26,3.9,4.64,5.48,6.42,7.46,8.6]

    # tickness at soil layer - delta zi - from CLM user's guide Table 2.3 
    levgrnd_thickness = [0.02,0.04,0.6,0.8,0.12,0.16,0.2,0.24,0.28,0.32,0.36,0.4,0.44,0.54,0.64,0.74,0.84,0.94,1.04,1.14,2.39,4.676,7.635,11.14,15.115]

    # depth at soil layer zh,i  - from CLM user's guide Table 2.3 
    levgrnd_depth = [0.02,0.06,0.12,0.2,0.32,0.48,0.68,0.92,1.2,1.52,1.88,2.28,2.72,3.26,3.9,4.64,5.48,6.42,7.46,8.6,10.99,15.666,23.301,34.441,49.556]

    return levsoi_thickness, levsoi_depth, levgrnd_thickness , levgrnd_depth 


# load inputdata and fix longitudes and latitudes
def load_inputdata(scenario, input_variables, case_ctl): 

    inputdir = procdir+'surfdata_4p1000/'
    surfat_filename = 'surfdata_360x720cru_16pfts_Irrig_CMIP6_simyr2000_c170824_CTL.nc'
    surfdat_filename_scen = 'surfdata_360x720cru_16pfts_Irrig_CMIP6_simyr2000_c170824_'+scenario+'.nc'

    # load input data and altered input data
    ds_input = import_ds(inputdir+surfat_filename, myVars=input_variables,rename_lsmlatlon=True)
    ds_input_scen = import_ds(inputdir+surfdat_filename_scen, myVars=['ORGANIC'],rename_lsmlatlon=True)

    # load for lat and lons
    ds_lonlat = import_case(case_ctl, myVars=['EFLX_LH_TOT'], timeSlice=slice(str(1995)+"-01-01", str(1995)+"-02-01"), alter_lons=False)

    # adapt lons and lats
    ds_input['lon'] = ds_lonlat['lon']
    ds_input['lat'] = ds_lonlat['lat']

    ds_input_scen['lon'] = ds_lonlat['lon']
    ds_input_scen['lat'] = ds_lonlat['lat']

    ds_input = lon_pm2idl(ds_input)
    ds_input_scen = lon_pm2idl(ds_input_scen)

    # load individual variables
    landmask = ds_input['PFTDATA_MASK']

    return ds_input, ds_input_scen, landmask


# load precalculated wilting point datasets (in postprocessing) 
def load_wilting_h2osoi(scenario): 

    ds_wilting_h2osoi = xr.open_dataset(procdir+'/postprocessing/wilting_volumetric_soil_water_CTL.nc')
    ds_wilting_h2osoi_scen = xr.open_dataset(procdir+'/postprocessing/wilting_volumetric_soil_water_'+scenario+'.nc')

    da_wilting_h2osoi = ds_wilting_h2osoi['H2OSOI_WILT']
    da_wilting_h2osoi_scen = ds_wilting_h2osoi_scen['H2OSOI_WILT']
    
    return da_wilting_h2osoi, da_wilting_h2osoi_scen

# calculate number of gridcells in region based on mask (eg for sandy and clay soils)
def calc_ngridcells_ofmask_inregion(da_mask, mask, region_ids): 
    da_region_list = []
    regions = []
    region_abbrevs = []

    for i,region in enumerate(region_ids):

        da_region = da_mask.where(mask==region).sum(dim=('lat','lon'))

        # only continue if pft is found in region
        if not da_region.isnull().sum().values ==da_region.size: 

            region_name = regionmask.defined_regions.ar6.land[region].name
            region_abbrev = regionmask.defined_regions.ar6.land[region].abbrev

            regions.append(region)
            region_abbrevs.append(region_abbrev)
            da_region_list.append(da_region)

    da_regions = xr.concat(da_region_list, 'region')
    da_regions['region'] = regions
    
    return da_regions


# calculare regional averages per PFT for different textures -- only for work with pfts
def calc_reldelta_regionalmean(da, da_scen, region_ids, mask_regions, texture_mask=None):

    da_region_list = []
    regions = []
    for i,region in enumerate(region_ids):


        if len(texture_mask)>1: 

            da_scen= da_scen.where(texture_mask)
            da = da.where(texture_mask)

        da_region = da.where(mask_regions==region).mean(dim=('lat','lon'))

        # only continue if pft is found in region
        if not da_region.isnull().sum().values ==da_region.size: 
            da_scen_region = da_scen.where(mask_regions==region).mean(dim=('lat','lon'))

            region_name = regionmask.defined_regions.ar6.land[region].name
            region_abbrev = regionmask.defined_regions.ar6.land[region].abbrev


            da_delta_region = (da_scen_region - da_region)/da_region *100
            regions.append(region)
            da_region_list.append(da_delta_region)

    ds_delta_regions = xr.concat(da_region_list, 'region')
    ds_delta_regions['region'] = regions

    return ds_delta_regions

# calculate delta and relative delta for pft output gridded variables
def calc_delta_and_rel_pft(case_ctl,case_scen,variable,pft_list):

    list_da_pft = []
    list_da_pft_scen = []

    for pft_int in pft_list: 

        list_da_pft.append(load_da_gridded_pft(variable, pft_int, case_ctl).rename(variable))
        list_da_pft_scen.append(load_da_gridded_pft(variable, pft_int, case_scen).rename(variable))

    ds_ctl = xr.concat(list_da_pft, dim='ivt_str').rename({'ivt_str':'pft'})
    ds_ctl['pft'] = pft_list

    ds_scen = xr.concat(list_da_pft_scen, dim='ivt_str').rename({'ivt_str':'pft'})
    ds_scen['pft'] = pft_list

    ds_delta = ds_scen - ds_ctl
    ds_delta_rel = (ds_delta/ ds_ctl)*100

    return ds_delta, ds_delta_rel


# calculate regional mean relative differences of  two data arrays -- only for work with pfts
def calc_regionalmeans(da, region_ids, mask):
    
    da_region_list = []
    regions = []
    for i,region in enumerate(region_ids):
        
        print('processing region: '+str(i+1) +' of '+str(len(region_ids)), end='\r')

        da_region = da.where(mask==region).mean(dim=('lat','lon'))

        # only continue if pft is found in region

        region_name = regionmask.defined_regions.ar6.land[region].name
        region_abbrev = regionmask.defined_regions.ar6.land[region].abbrev

        regions.append(region)
        da_region_list.append(da_region)

    ds_regions = xr.concat(da_region_list, 'region')
    ds_regions['region'] = regions
    
    return ds_regions

# calculate delta and relative delta for pft output gridded variables -- for a single pft variable
def calc_delta_and_rel_pft_3d(case_ctl, case_scen, variable, pft_int, flag_singleyear = False):
    
    da_pft_scen = load_da_gridded_pft(variable, pft_int, case_scen).rename(variable)
    da_pft_ctl = load_da_gridded_pft(variable, pft_int, case_ctl).rename(variable)
    
    
    if 'ivt_str' not in da_pft_scen.coords: 
        da_pft_scen = da_pft_scen.assign_coords({"ivt_str": "weighted sum of c3 and c4 grasses"})
    
    if 'ivt_str' not in da_pft_ctl.coords: 
        da_pft_ctl = da_pft_ctl.assign_coords({"ivt_str": "weighted sum of c3 and c4 grasses"})    
        
        
    # select only single year to save memory
    if not isinstance(flag_singleyear, bool): 
        da_pft_scen = da_pft_scen.loc[str(flag_singleyear[0])+"-01-01":str(flag_singleyear[1])+"-12-01"]
        da_pft_ctl  = da_pft_ctl.loc[str(flag_singleyear[0])+"-01-01":str(flag_singleyear[1])+"-12-01"]

    da_delta = da_pft_scen - da_pft_ctl
    da_delta_rel = (da_delta/ abs(da_pft_ctl))*100
    
    del da_pft_scen, da_pft_ctl
    return da_delta, da_delta_rel


# grid different variables and different pfts into one dataset
def grid_multiple_pfts_multiple_vars(ds,  variables, pft_list, case, flag_daily=False): 
    
    for variable in variables: 
        
        print('processing variable '+str(variable))
        
        for pft_int in pft_list:# range(0,len(pftlist)): 

            print('processing pft '+str(pft_int))

            da_gridded   = grid_one_variable_col(ds, variable, fillValue=None, vegtype=pft_int)

            da_gridded   = lon_pm2idl(da_gridded, fail_silently=False)

            ds_gridded = xr.Dataset()

            ds_gridded[variable+'_'+str(pft_int)] = da_gridded

            if flag_daily: 
                text_daily = '_daily'
            else: 
                text_daily = ''

            procdir = '/capstor/scratch/cscs/ivanderk/processing_4p1000/' #'/project/s1207/ivanderk/scratch/' 

            ds_gridded.to_netcdf(procdir+'/postprocessing/pft_gridding/'+variable+'_'+str(pft_int)+'.'+case+text_daily+'.nc')
            print(procdir+'/postprocessing/pft_gridding/'+variable+'_'+str(pft_int)+'.'+case+text_daily+'.nc')

            del da_gridded, ds_gridded            

            
# calculate cumulative water stress for an average year
def calc_annual_cumulative_stress(da, da_wilting_h2osoi, mask=False): 

    
    da_wilting_delta = (da_wilting_h2osoi - da)
    if not mask==False: 
        da_wilting_delta = da_wilting_delta.where(mask)
    
    da_stress = da_wilting_delta.where(da_wilting_delta>0)
    da_stress_cumulative = da_stress.groupby('time.year').sum()
    
    return da_stress_cumulative


# per pft, calculate and save regional mean delta per soil level and region
def save_regionalmean_delta_and_delta_rel(variable, pft_list, case_ctl, case_scen, region_ids, mask_regions, scenario): 
    
        
    for pft_int in pft_list: 
        print('processing pft: '+str(pft_int))

        # calculate delta and relative delta of saturation ratio. 
        da_delta, da_delta_rel = calc_delta_and_rel_pft_3d(case_ctl, case_scen, variable, pft_int)
        
        if variable =='SOILPSI': 

            da_delta = da_delta.isel(levgrnd=slice(0,10)).rename({'levgrnd':'levsoi'}).squeeze()
            da_delta_rel = da_delta_rel.isel(levgrnd=slice(0,10)).rename({'levgrnd':'levsoi'}).squeeze()

        # calculate regional means. 
        ds_regions = calc_regionalmeans(da_delta, region_ids, mask_regions)

        del da_delta

        ds_regions.to_netcdf(procdir+'/postprocessing/pft_gridding/regional_means.'+variable+'.'+str(pft_int)+'.'+scenario+'.delta_change.nc')

        # calculate regional means. 
        ds_regions_rel = calc_regionalmeans(da_delta_rel, region_ids, mask_regions)

        del da_delta_rel
        ds_regions_rel.to_netcdf(procdir+'/postprocessing/pft_gridding/regional_means.'+variable+'.'+str(pft_int)+'.'+scenario+'.delta_change_rel.nc')

        print()

# load postprocessed annual drydays
def load_da_drydays_ymean(variable, case, stream, start_year, end_year, save_text): 
    filename_drydays = procdir+'postprocessing/dryday_ymean/'+variable+'.'+case+'.'+stream+'.'+str(start_year)+"-"+str(end_year)+'_'+save_text+"_ymean.nc"
        
    return xr.open_dataset(filename_drydays)[variable]


# combine the different da's in one dataset
def load_ds_drydays_ymean(variables, case, stream, start_year, end_year, save_text='dryday'): 
    for i, variable in enumerate(variables): 
        da  = load_da_drydays_ymean(variable, case, stream, start_year, end_year, save_text=save_text)
        if i ==0: 
            ds = da.to_dataset(name=variable)
        else: 
            ds[variable] = da
    return ds

# calculate or load if not existing the time mean of the delta and control dataset
def load_timmean_dataset(case_ctl, case_scen, variables, stream, start_year, end_year, flag_recalulate=False):
    
    filename_delta_timmean  = procdir+'postprocessing/timmean/'+case_ctl+'_delta.'+stream+'.'+str(start_year)+"-"+str(end_year)+"_timmean.nc"
    filename_ctl_timmean = procdir+'postprocessing/timmean/'+case_ctl+'.'+stream+'.'+str(start_year)+"-"+str(end_year)+"_timmean.nc"

    # do calculation
    if (not os.path.isfile(filename_delta_timmean) and not os.path.isfile(filename_ctl_timmean)) or flag_recalulate:
       
        time_slice = slice(str(start_year)+"-01-01", str(end_year)+"-12-31")
        ds_ctl = import_case(case_ctl, myVars=variables, timeSlice=time_slice)
        ds_scen = import_case(case_scen, myVars=variables, timeSlice=time_slice)

        da_delta = ds_scen - ds_ctl
        ds_delta_timmean = (ds_delta).mean('time').compute() 
        ds_ctl_timmean = ds_ctl.mean('time').compute()

        
        ds_delta_timmean.to_netcdf(filename_delta_timmean)
        ds_ctl_timmean.to_netcdf(filename_ctl_timmean)

    # load data
    else: 
        
        ds_delta_timmean = xr.open_dataset(filename_delta_timmean)
        ds_ctl_timmean = xr.open_dataset(filename_ctl_timmean)
    
    return ds_delta_timmean, ds_ctl_timmean

# check if in hydrol variable list and if so, convert units
def conv_hydrol_var(ds_in, variables):
    ds_out = ds_in
    for variable in variables: 
        da_in = ds_in[variable]

        if da_in.units == 'mm/s':
            ds_out[variable] = conv_mm_s_to_m_year(da_in)

        elif da_in.units == 'W/m^2': 
            ds_out[variable] = conv_W_per_m2_to_m_year(da_in)
            
        elif da_in.units =='kg/m2' or da_in.units == 'mm': 
            da_out_values = da_in.values  * 0.001 # m
            ds_out[variable] = xr.DataArray(da_out_values, coords=da_in.coords,
         dims=da_in.dims)
            ds_out[variable].attrs['units'] = 'm'
            

        else: 
            ds_out[variable] = da_in

    return ds_out

# conversion function
def conv_mm_s_to_mm_day(da_in):

    da_out_values = da_in * 86400  
    # update attributes and change units
    da_out = xr.DataArray(da_out_values, coords=da_in.coords,
         dims=da_in.dims)
    da_out.attrs= da_in.attrs
    da_out.attrs['units'] = 'mm/day' 

    return da_out

# conversion function
def conv_W_per_m2_to_m_year(da_in):
    
    # latent heat of vaporization
    lamda = 2.501 * 10**6 #J/kg
    da_out_values = da_in.values *(1/(lamda * 1000) )  * 86400 * 365 # m/year
    da_out = xr.DataArray(da_out_values, coords=da_in.coords,
         dims=da_in.dims)
    # update attributes and change units
    da_out.attrs = da_in.attrs
    da_out.attrs['units'] = 'm/year' 

    return da_out

# conversion function
def conv_mm_s_to_m_year(da_in):

    da_out_values = da_in.values * 86400  * 365 * 0.001

    da_out = xr.DataArray(da_out_values, coords=da_in.coords,
         dims=da_in.dims)
    # update attributes and change units
    da_out.attrs= da_in.attrs
    da_out.attrs['units'] = 'm/year' 

    return da_out

# ---------------------------------------------------------------------
# 3. Functions to plot
# ---------------------------------------------------------------------

def set_plot_param():
    """Set my own customized plotting parameters"""
    
    import matplotlib as mpl
    mpl.rc('axes',edgecolor='grey')
    mpl.rc('axes',labelcolor='dimgrey')
    mpl.rc('xtick',color='dimgrey')
    mpl.rc('xtick',labelsize=12)
    mpl.rc('ytick',color='dimgrey')
    mpl.rc('ytick',labelsize=12)
    mpl.rc('axes',titlesize=14)
    mpl.rc('axes',labelsize=12)
    mpl.rc('legend',fontsize='large')
    mpl.rc('text',color='dimgrey')


def plot_ts(da):
    """ plot timeseries (of spatial mean) """
    da_ts = da.mean(dim=('lon','lat'))
    da_ts.plot()
    plt.title(da.long_name + '('+da.units+')')
    plt.xlim(da_ts.time[0].values,da_ts.time[-1].values)


def plot_ymean(da):
    """calculate and plot annual mean timeseries """
    # check if input da is already ymean, otherwise do calculation 
    if len(da) < 500: 
        da_ymean = da.mean(dim=('lon','lat'))
    else: 
        da_ymean = da.mean(dim=('lon','lat')).groupby('time.year').mean('time')
    
    xlims = (da_ymean.year[0].values,da_ymean.year[-1].values)
    da_tseries = da_ymean.plot(xlim=xlims)
    plt.title(da.long_name, pad=5)
    plt.ylabel(da.name+' [' + da.units + ']')
    #plt.plot([da_ymean.year[0],da_ymean.year[-1]], [0,0], linewidth=1, color='gray')

# calculate and plot annual sum timeseries 
def plot_ysum(da):
    
    da_ymean = da.sum(dim=('lon','lat')).groupby('time.year').mean('time')
    xlims = (da_ymean.year[0].values,da_ymean.year[-1].values)
    da_tseries = da_ymean.plot(xlim=xlims)
    plt.title(da.long_name+' [' + da.units + ']' )


# plot timmean per selected region
def plot_yts_sel_regions(da_to_mask, selected_regions):
    mask = regionmask.defined_regions.srex.mask(da_to_mask)
    
    # annual means are already calculated
    if len(da_to_mask) < 50:
        da_mask_ts = da_to_mask.groupby(mask).mean('stacked_lat_lon').sel(region=selected_regions)
    else: 
        da_mask_ts = da_to_mask.groupby(mask).mean('stacked_lat_lon').groupby('time.year').mean().sel(region=selected_regions)

    # add abbreviations and names
    abbrevs = regionmask.defined_regions.srex[da_mask_ts.region.values].abbrevs
    names = regionmask.defined_regions.srex[da_mask_ts.region.values].names
    da_mask_ts.coords['abbrevs'] = ('region', abbrevs)
    da_mask_ts.coords['names'] = ('region', names)
    
    f, axes = plt.subplots(3, 2, figsize=(8,5))
    f.suptitle(da_to_mask.name, fontsize=14)

    low = da_mask_ts.min()
    high = da_mask_ts.max()
    for i in range(len(selected_regions)):
        (nx,ny) = axes.shape
        if i < nx : ax = axes[i,0]
        else      : ax = axes[i-nx,1]
            
        ts_region = da_mask_ts.isel(region=i)
        ts_region.plot(ax=ax)
        ax.set_title(da_mask_ts.isel(region=i).names.values)
        ax.set_ylim(low,high)
        ax.set_ylabel('('+da_to_mask.units+')')
        ax.set_xlim(da_mask_ts.year[0],da_mask_ts.year[-1])
        ax.plot([da_mask_ts.year[0],da_mask_ts.year[-1]], [0,0], linewidth=1, color='gray')
    
    plt.setp(axes, xlabel="")

    f.tight_layout()


# plot global map of difference 
def plot_delta_map(da_delta, plot_regions=False, vlims=False, calcsum=False, cmap='BrBG'):
    
    # calculate annual sum instead of mean (precip)
    if calcsum: 
        da_delta_ysum = da_delta.groupby('time.year').sum()
        da_delta_mean = da_delta_ysum.mean('year')
        da_delta_mean.attrs['units'] = 'mm/year'
    # only one value
    elif len(da_delta.dims) < 3: 
        da_delta_mean = da_delta
    # annual means already taken
    elif 'year' in da_delta.dims:
        da_delta_mean = da_delta.mean('year')
    else:
        da_delta_mean = da_delta.mean('time')
    
    plt.figure(figsize=(12,5))
    proj=ccrs.PlateCarree()
    ax = plt.subplot(111, projection=proj)
        
    # limiting values for plotting are given    
    if vlims==False: 
        da_delta_mean.plot(ax=ax, cmap=cmap, cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')', 'fraction': 0.02, 'pad': 0.04})
    else: 
        da_delta_mean.plot(ax=ax, cmap=cmap, vmin=vlims[0], vmax=vlims[1], extend='both',  cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')', 'fraction': 0.02, 'pad': 0.04}, add_labels=False)
        
    ax.set_title(da_delta.long_name, loc='right')
    ax.coastlines(color='dimgray', linewidth=0.5)
    # exclude Antactica from plot
    ax.set_extent((-180,180,-63,90), crs=proj) 
    ax.axis('off')
    
    if plot_regions: regionmask.defined_regions.srex.plot(ax=ax,add_ocean=False, coastlines=False, add_label=False) #label='abbrev'
    return ax



f, ax = plt.subplots()

# plot global map of difference 
def plot_delta_map_noax(ax, da_delta, plot_regions=False, vlims=False, calcsum=False, cmap='BrBG'):
    """plot difference maps without creating a figure within function"""
    # calculate annual sum instead of mean (precip)
    if calcsum: 
        da_delta_ysum = da_delta.groupby('time.year').sum()
        da_delta_mean = da_delta_ysum.mean('year')
        da_delta_mean.attrs['units'] = 'mm/year'
    # only one value
    elif len(da_delta.dims) < 3: 
        da_delta_mean = da_delta
    # annual means already taken
    elif 'year' in da_delta.dims:
        da_delta_mean = da_delta.mean('year')
    else:
        da_delta_mean = da_delta.mean('time')
    
    # limiting values for plotting are given    
    if vlims==False: 
        da_delta_mean.plot(ax=ax, cmap=cmap, cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')'})
    else: 
        da_delta_mean.plot(ax=ax, cmap=cmap, vmin=vlims[0], vmax=vlims[1], extend='both',  cbar_kwargs={'label': da_delta.name+' ('+da_delta.units+')'}, add_labels=False)
        
    ax.set_title(da_delta.long_name, loc='right')
    ax.coastlines(color='dimgray', linewidth=0.5)
    # exclude Antactica from plot
    ax.set_extent((-180,180,-63,90)) 

    if plot_regions: regionmask.defined_regions.srex.plot(ax=ax,add_ocean=False, coastlines=False, add_label=False) #label='abbrev'
    return ax


# calculate and plot difference map of variable
def calc_plot_delta(ds_delta, ds_ctl, variable,d_vlims): 
    
    
    if variable in d_vlims.keys(): 
        da_delta = ds_delta[variable]
        if 'levsoi' in da_delta.dims: 
            da_delta = da_delta.weighted(da_delta['levsoi']).mean(dim='levsoi')
        elif 'levgrnd' in da_delta.dims: 

            da_delta = da_delta.weighted(da_delta['levgrnd']).mean(dim='levgrnd')

        da_delta.attrs = ds_ctl[variable].attrs


        plot_delta_map(da_delta, plot_regions=False, vlims=d_vlims[variable], calcsum=False, cmap='BrBG')
        
    else: 
        print(variable +' not in limits')
              
# calculate and plot relative difference map of variable
def calc_plot_delta_rel(ds_delta, ds_ctl, variable, d_vlims): 
    
    if variable in d_vlims.keys(): 

        da_delta = ds_delta[variable] / abs(ds_ctl[variable]) *100
        if 'levsoi' in da_delta.dims: 
            da_delta = da_delta.weighted(da_delta['levsoi']).mean(dim='levsoi')
        elif 'levgrnd' in da_delta.dims: 

            da_delta = da_delta.weighted(da_delta['levgrnd']).mean(dim='levgrnd')

        da_delta.attrs = ds_ctl[variable].attrs
        da_delta.attrs['units'] = '%'


        ax = plot_delta_map(da_delta, plot_regions=False, vlims=d_vlims[variable], calcsum=False, cmap='BrBG')    
        return ax
    else: 
        print(variable +' not in limits')
        
        
# calculate and plot difference map of variable
def calc_plot_delta_da(da_delta, da_ctl, variable,d_vlims): 

    if variable in d_vlims.keys(): 
        if 'levsoi' in da_delta.dims: 
            da_delta = da_delta.weighted(da_delta['levsoi']).mean(dim='levsoi')
        elif 'levgrnd' in da_delta.dims: 

            da_delta = da_delta.weighted(da_delta['levgrnd']).mean(dim='levgrnd')

        da_delta.attrs = da_ctl.attrs


        plot_delta_map(da_delta, plot_regions=False, vlims=d_vlims[variable], calcsum=False, cmap='BrBG')
        
    else: 
        print(variable +' not in limits')

        
              
# calculate and plot relative difference map of variable
def calc_plot_delta_rel_da(da_delta, da_ctl, variable, d_vlims): 
    
    if variable in d_vlims.keys(): 

        da_delta = da_delta / da_ctl *100
        if 'levsoi' in da_delta.dims: 
            da_delta = da_delta.weighted(da_delta['levsoi']).mean(dim='levsoi')
        elif 'levgrnd' in da_delta.dims: 

            da_delta = da_delta.weighted(da_delta['levgrnd']).mean(dim='levgrnd')

        da_delta.attrs = da_ctl.attrs
        da_delta.attrs['units'] = '%'


        ax = plot_delta_map(da_delta, plot_regions=False, vlims=d_vlims[variable], calcsum=False, cmap='BrBG')    
        return ax
    else: 
        print(variable +' not in limits')
        