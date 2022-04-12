# In Rainbows: Spatially-resolved SED fitting of local galaxies in Dustpedia

![](https://github.com/AndreaEnia/InRainbows/blob/672e25884e6a7d303ecd6efc0428195924f13e63/seds.mp4)

Upon which the results in [Enia et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.493.4107E/abstract), [Morselli et al., 2020](https://ui.adsabs.harvard.edu/abs/2020MNRAS.496.4606M/abstract), and [2021](https://ui.adsabs.harvard.edu/abs/2021MNRAS.502L..85M/abstract) are based.


## Introduction

To run a galaxy (i.e. NGCxxxx), there are three simple steps to follow:

- Create the galaxy folder, with the same name as the galaxy (NGCxxxx, mandatory)
- cut'n'paste inside the folder the notebook "GalaxyFlow.ipynb" found in scripts/
- follow the instructions inside the notebook

The notebook is a series of simple istructions, as it creates step by step the necessary folders with the reduced maps, results and plots. Check under for an exhaustive expalanation of what each function does, and how to customize them, in the meantime a general overview.

**Preliminary steps.**

1) Download from the DustPedia website bands and associated errors (if available) the user wants to work with. I suggest to pre-check in the site the bands and errors availability.

2) Data reduction, and for data reduction I mean: sky subtraction (if present), degradation to the worst map resolution, which as of today is assumed to always be SPIRE350. Yes, this thing could easily be generalized, but I don't tengo cazzi. This step is time consuming, so the reduced maps are saved in a folder, as in this way there is no need to always redo the whole procedure.

3)  Measure the maps rms (per pixel).

**Photometric steps.**

4) The user, depending on what results it desires, has to simply fill there four variables:
```python
* run_type = # The run name, i.e. pBp for pixel-by-pixel, or 1kpc for 1kpc apertures, o "Bh<oalf3-2Hhlad" if it feels creative.
* aperture_side = # How big the SIDE of the (square) apertures should be, either in parsec, kiloparsec or arcsec.
* reference_map_4grid = # The reference map for the apertures grid generation (coordinates.txt with values in ra e dec)
* reference_map_4check = # The reference map to check how the apertures were created 
* grid_extention_from_center = # How long the apertures grid should be extended from the galaxy center, either in pc, kpc or arcsec
* plot_radius = # How big the grid check plot should be. I suggest slightly bigger than the former grid_extention
```

5) Now it's the time for the big cell with all the photometric steps: the grid is generated, checked (eventually), the code performs photometry on that particular grid, builds the table that enters MAGPHYS as the input for SED fitting. This table is not built on all the apertures, but only on the ones over a certain threshold of positive photometric points (10 by default), given that if SNR goes below another threshold (1.5 by default), that point is excluded from photometry.

**SED fitting step and results.**

6) The table is copied within magphys folder, a bash script runs magphys, the user waits until it is onver, and then copies the results (.sed, .fit) inside the folder 'magphys'+_run_type. These are automatically read, saving the results in a comfy magphys_results.csv

7) One optional step, given the presence of sparse points disconnected from the other apertures given the thresholds reported earlier. This is a clustering removal step with DBScan, which tries to identify the main cluster of points belonging to the bigger group of apertures, discarding the isolated ones.

8) Finally, plots the results. These quick'n'dirty plots are quite ugly and preliminar, just to see what gets out.

Pedantic description of all the scripts.
-----------
All the described work is exectued by various scripts, contained in the folder scripts/ (duh!), with pretty self-explanatory titles:
* `GenericUsefulScripts` (GUS)
* `DataReduction` (DataRed)
* `PhotometryScripts` (PhotoScripts)
* `SedFittingScripts` (SEDfit)
* `ResultsScripts` (ResScripts)

Let's start over.
### GUS
GUS is a container of pretty much everything: scripts that download the maps (`GUS.download_data`, `GUS.download_data_errors`), classes that read galaxies physical properties (`GUS.GalaxyProperties`), to read the map (`GUS.FitsUtils`) and to convert the physical scales into angular and viceversa (`GUS.AngularPhysicalConv`), scripts to avoing ambiguity in observing band names (`GUS.band_disambiguation`) and to associate a colormap to each band (`GUS.associate_colormap`), and so on.

### DataRed
DataRed is the script containing all the data reduction processes (and to verify everything).
To reduce them, the syntax is `DataRed.data_reduction(galaxy_name, path_fits_input = 'Caapr/Maps')`. The input path depends on wheter the user wants CAAPR to run with star subtraction. I never run CAAPR with star subtraction, so I automatically inserted the path of the folder where the maps are downloaded. This generates a folder *_ReducedMaps* where the reduced maps are stored.

The check is performed with `DataRed.check_Dustpedia(galaxy_name, working_bands)`, and trivially measures the fluxes of the reduced maps within the same aperture used by the Dustpedia collaboration, and confronts it with their reported fluxes (all information that are within the folder *DustPedia_Tables*)

If a map reduction goes awry (i.e. fluxes completely different from the Dustpedia ones), or one feels the need to re-run the whole reductino process, the uses should not get frightened, as the script automatically skips the bands for which the reduced map is inside the folder. Therefore if that's the case, the user should just delete the bad map, and re-run the cell.

### PhotoScripts
PhotoScripts contain all those scripts that measure the maps rms, generate the grid of coordinates and perform photometry on those apertures, and finally build the photometric tables, including the ones (observations.dat and filters.dat) to give to MAGPHYS as input.

As for the **rms**, there are two versions (v1 and v2) within `PhotoScripts.evaluate_rms(galaxy_properties)`. Both versions share the majority of code, only the final part differs. For "majority of code" I mean: galaxy masking as prescribed by DustPedia tables, sigma clipping of the masked map, generation of 2500 (by default, N_ap = xxx to choose a number) random apertures of 10 arcsec radius (by default, ap_radius_physical = xxx to choose a value). Now the two methods diverge:
- v2 is adapted from CAAPR, with the rms measured on each of the 2500 apertures, and the map rms assumed as the peak of the distribution
- v1 takes the distribution of flux in all the 2500 apertures, and measures the rms as the standard deviation
In both cases, the rms is per pixel, so to associate an error to a measured flux within an aperture this rms is multiplied per the number of pixels in the aperture.

Now the actual photometric part, where the process is quite customizable. The user fills the variables cited before, and the scripts that run are:
- the one building the coordinates grid `PhotoScripts.generate_coordinate_grid_within_DustAp`. This is generated (in an artisanal way that could be greatly improved) within a circle of a given radius centered on the galaxy center. The user can also give an ellipticity, a position angle, i.e. for elongated galaxies. Within Dustpedia Aperture means that if one of angle, ellipticity or avoidance_semiaxes is False, it assumes the DustPedia aperture values, the one over which the collaboration measured the photometry;
- sometimes the grid is too big, or too small, so the validation script is: `PhotoScripts.check_coordinate_grid`;
- finally, to perform the photometry map-by-map on that grid of apertures, `PhotoScripts.do_photometry_build_table`. There are two thresholds that can be set: a SNR threshold (1.5 by default) under which all the photometric values with lower SNR are discarded and substituted with -99, and one with the number of negative photometric over which the entire aperture is discardes from the SED fitting process (10 by default), something that usually happens in the galaxies outskirts. This script generates three tables: **photometries_table.csv**, containing all the measued photometries even in the apertures which not satisfy the thresholds, **GOOD_photometries_table.csv**, that discards those points, **Final_for_MAGPHYS.csv**, the latter but formatted in the way that MAGPHYS likes.

### SEDfit
Here things get quite artisanal. The original MAGPHYS version wants as an input observations and filters file, but the output is generated directly on the /home/magphys directory without giving the chanche to insert a path, and as such ***it is impossible to automatically run multiple galaxies at the same time***. These last things should be hand fixed, and that's what those scripts here, `SEDfit.prepare_magphys`, `SEDfit.run_magphys`, `SEDfit.move_results` do.

**!!! CAVEATS !!!**

- `SEDfit.prepare_magphys` simply creates within /home/**whatever**/magphys a folder named as the galaxy, and a subfolder with the run name, inside which it puts  **filters.dat**, and partition the file **Final_for_MAGPHYS.csv** in many subtables **observations_xxx_xxx.dat** depending on how many parallel run you want to run: fill_cores, almost_fill_cores (all cores minus two), single (no parallelization), or a user inputted entries per file. If necessary, some filters might be excluded from SED fitting.
- `SEDfit.run_magphys` runs MAGPHYS. Now, this seems like a piece of cake, but its not. The original MAGPHYS was not thought to run on multiple resolved regions of the same galaxy, but on a sample of multiple galaxies at different redshift, and as such it follows this sequence: first it generates the sample redshift library, then the star+dust libraries for each associated redshift, and finally it fits the sample. But here I have up to 10k points at the same z, and MAGPHYS does not understand that it is the same galaxy, so it would generate 10k star+dust libraries overwriting the same file each time. As such, I have to do this: generate a obs_zero.dat with a single line, generate star+dust libraries on this line, **wait until these are created** with a brutal "stop for 300 seconds, then do the rest", then open a screen for every parallel run, modify the configuration file .magphys_tcshrc accordingly, and run fit_sample. All this is automatically done by the script.
- `SEDfit.run_magphys` checks that MAGPHYS SED fitting procedure is over, then move the results on the appropiate folders. How does it knows that MAGPHYS has stopped? As it runs, MAGPHYS generates 0 bytes size files, xxxx.fit, inside which the bestfit model is going to be wfitten. As long as 0 bytes files are in the path, MAGPHYS is still running, and that script will tell you to come back later.

### ResScripts
All those scripts which read the results from magphys, and plots the results.
