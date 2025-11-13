# Kruger_GEDI_paper
Analysis for the paper on GEDI canopy metric mapping in Kruger

The general process for running the analysis was as follows:
1. Install miniconda and use ecofor_python_env.yml to create python environment. Install R and packages used in each R script.
2. Run gedi_download.R to get GEDI footprints in study area and time period
3. Run data_prep.ipynb to create most predictor layers. CCDC layers were created in through the javascript interface of Google Earth Engine with run_ccdc_tile.js.
4. Extract the predictors for the GEDI footprints with gedi_raster_extraction_bydate.r, gedi_raster_extraction_static.r, and ccdc_sample_extraction.js.
5. Prepare field data and create models with veg.ipynb (run up through the modeling section).
6. Use the chosen model to create maps of each GEDI metric with map_concurrent_multimodel.py.
7. Run the Maps section of veg.ipynb.
8. Create small area estimates with gedi_sae_emdi.r, and compare random effects structures with gedi_mixedmodels.r.
9. Create figures for paper by running Figures section of veg.ipynb.

Using Google Earth Engine in data_prep.ipynb and veg.ipynb requires the python GEE utility functions in <https://github.com/VogelerLab/pee>.
