import os, shutil, subprocess, rasterio, time
os.environ['USE_PYGEOS'] = '0'
import numpy as np
from glob import glob
import pandas as pd
import geopandas as gpd
from sklearn.ensemble import RandomForestRegressor
from scipy.ndimage import uniform_filter
from joblib import load
import concurrent.futures

kernel_rad = 0  # radius of window for focal mean. 0 for no focal mean. 
out_nodata = np.nan
out_dtype = np.float32

# Mask must match the predictor set or the results are invalid
mask_path = r"D:\ECOFOR\boundaries\greaterkruger_utm36n.tif"
  
def buffer_window(window, big_window, overlap=1):
    window = rasterio.windows.Window(
                col_off=window.col_off - overlap,
                row_off=window.row_off - overlap,
                width=window.width + overlap * 2,
                height=window.height + overlap * 2)
    window = window.intersection(big_window)
    return window

def map_window(w, rf, paths_dict):
    """ Apply trained model to a window (w) of the rasters using Landsat data from lt_path.
    """

    np.seterr(all="ignore") # suppress warning when reading in NaNs for topo
    wshape = (w.height, w.width)
    feature_names = rf.feature_names_in_
    
    if kernel_rad>0:
        w_buf = buffer_window(w, big_window, overlap=kernel_rad)
        rstart = w.row_off-w_buf.row_off
        cstart = w.col_off-w_buf.col_off
    
    # load mask
    with rasterio.open(mask_path) as src:
        msk = src.read(1, window=w).astype(bool)
    
    if not msk.any():
        ypred = np.full(wshape, out_nodata, out_dtype)
        return ypred
    
    # prep each of the arrays to make a prediction with
    mask_list = []
    dfs = []
    for k, path_dict in paths_dict.items():
        # check if any bands are in the predictor layers
        path = path_dict['path']
        prefix = path_dict['band_prefix']
        with rasterio.open(path) as src:
            bdict = {i+1:prefix+b for i, b in enumerate(src.descriptions) if prefix+b in feature_names}
            bixs, bands = list(bdict.keys()), list(bdict.values())
            nodata_val = src.nodata
            
            if len(bdict)==0:
                print('None of the bands in', path,' are in the model predictor set. Quitting.')
                return None            
            
            if kernel_rad>0:
                parr = src.read(bixs, window=w_buf)

                # normalized convolution to deal with missing data
                mask = parr!=src.nodata
                parr = uniform_filter(parr * mask, size=[1,1+kernel_rad*2,1+kernel_rad*2], output=np.float32, mode='nearest')
                weights = uniform_filter(mask, size=[1,1+kernel_rad*2,1+kernel_rad*2], output=np.float32, mode='nearest')
                parr /= weights
                parr[~mask] = nodata_val

                # get unbuffered array
                parr = parr[:, rstart:rstart+w.height, cstart:cstart+w.width]
            else:
                parr = src.read(bixs, window=w)
                
        shape = parr.shape
        parr = parr.reshape(shape[0], shape[1]*shape[2]).T
        pdf = pd.DataFrame(parr, columns=bands)
        pdf.loc[~msk.ravel(),:] = nodata_val
        if np.isnan(nodata_val):
            pnodata = np.isnan(pdf).any(axis=1)
        else:
            pnodata = (pdf==nodata_val).any(axis=1)
    
        mask_list.append(pnodata)
        dfs.append(pdf)
    
    mask = ~(pd.concat(mask_list, axis=1).any(axis=1))
    if mask.sum()==0:
        ypred = np.full(wshape, out_nodata, out_dtype)
    else:
        nullix = mask[~mask].index
        nulls = pd.Series(out_nodata, index=nullix, dtype=out_dtype)

        # predict valid data
        X = pd.concat(dfs, axis=1)
        X = X[feature_names]                             # need to make df same order as training data
        Xsub = X[mask]
        ypred = pd.Series(rf.predict(Xsub), index = Xsub.index) # mean of trees
        ypred = ypred.astype(out_dtype)

        # add in nulls
        ypred = pd.concat([ypred,nulls])
        ypred.sort_index(inplace=True)

        # reshape and convert dtype
        ypred = ypred.values.reshape(wshape)
    
    return ypred


if __name__=='__main__':
    
    # Setup for all runs
    lt_dir = r"J:\projects\ECOFOR\lt" 
    palsar_dir = r"J:\projects\ECOFOR\palsar\lt_tiling_scheme"

    # Mask must match the predictor set or the results are invalid
    mask_path = r"D:\ECOFOR\boundaries\greaterkruger_utm36n.tif"
    
    # processpoolexecuter params
    num_workers = 18
    chunksize = 50
    
    runs = [
        {'model':r"D:\ECOFOR\gedi\models\v08\GEDI_2AB_2019to2023_leafon_sampy500m_all_lt-p-s-t_cover_v08.joblib",
         'outdir':r"D:\ECOFOR\gedi\maps\v08\lt-p-s-t\cover",
         'outbasename':'cover_'
        },
        {'model':r"D:\ECOFOR\gedi\models\v08\GEDI_2AB_2019to2023_leafon_sampy500m_all_lt-p-s-t_rh98_v08.joblib",
         'outdir':r"D:\ECOFOR\gedi\maps\v08\lt-p-s-t\rh98",
         'outbasename':'rh98_'
        },
        {'model':r"D:\ECOFOR\gedi\models\v08\GEDI_2AB_2019to2023_leafon_sampy500m_all_lt-p-s-t_fhd_normal_v08.joblib",
         'outdir':r"D:\ECOFOR\gedi\maps\v08\lt-p-s-t\fhd",
         'outbasename':'fhd_'
        },
           ]

    for run in runs:
        # Copy model to output directory
        os.makedirs(run['outdir'], exist_ok=True)
        mod_path = os.path.join(run['outdir'], os.path.basename(run['model']))
        shutil.copy(run['model'], mod_path)
        rf = load(run['model'])

        # run mapping for each year of the predictor time series
        years = list(range(2007, 2011)) + list(range(2015, 2023))        
        for year in years:
            start=time.time()
            
            paths_dict = {
                          'lt_dry':{'path':os.path.join(lt_dir, "dry", "lt_dry_"+str(year)+".vrt"),
                                    'band_prefix':'lt_dry_'},
                          'lt_wet':{'path':os.path.join(lt_dir, "wet", "lt_wet_"+str(year)+".vrt"),
                                    'band_prefix':'lt_wet_'},
                          'palsar':{'path':os.path.join(palsar_dir, "palsar_"+str(year)+".vrt"),
                                    'band_prefix':'palsar_'},
                          'topo':{'path': r"D:\ECOFOR\topo\topo_all.vrt", 
                                    'band_prefix':'topo_'},
                          'soil':{'path': r"D:\ECOFOR\soils\soil_all.vrt",
                                    'band_prefix':'soil_'},
                        }
            
            outpath = os.path.join(run['outdir'], run['outbasename']+str(year)+".tif")
            if os.path.exists(outpath):
                print(outpath, 'already exists. Skipping.')
                continue
            print(outpath)
            
            # Get windows to apply process to
            with rasterio.open(mask_path) as src:
                profile = src.profile
                windows = [window for ij, window in src.block_windows()] 
                big_window = rasterio.windows.Window(col_off=0, row_off=0, width=src.meta['width'], height=src.meta['height'])
                profile['nodata'] = out_nodata
                profile['dtype'] = out_dtype
                profile['count'] = 1
            
            path_iter = [paths_dict]*len(windows)
            rf_iter = [rf]*len(windows)
            args = (windows, rf_iter, path_iter)

            with rasterio.open(outpath, 'w', **profile) as dst:                                    
                with concurrent.futures.ProcessPoolExecutor(max_workers=num_workers) as executor:
                    for window, result in zip(windows, executor.map(map_window, *args, chunksize=chunksize)):
                        dst.write(result, 1, window=window)

            end = time.time()
            print(end-start)
        