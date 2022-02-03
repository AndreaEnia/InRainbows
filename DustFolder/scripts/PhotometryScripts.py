#!/usr/bin/python
# -*- coding: latin-1 -*-

import os, sys, subprocess
import GenericUsefulScripts as GUS
import subprocess
import numpy as np
np.seterr(invalid='ignore')
import pandas as pd
DustPath = '/home/dustpedia/DustFolder/'
bandwavpath = DustPath+'scripts/bands_and_wvlghts.txt'
dict_df = pd.read_csv(bandwavpath, delimiter = '\t')
Wvlghts_dictionary = dict(zip(dict_df['name'], dict_df['lambda_eff']))
from astropy.io import fits, ascii
from astropy import units as u
from astropy.constants import c as v_lux
from astropy.cosmology import Planck15 as cosmo
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.wcs.utils import pixel_to_skycoord, skycoord_to_pixel
#from tqdm.autonotebook import trange, tqdm
from tqdm import trange, tqdm
space = '       '

def evaluate_rms(galaxy_properties, N_ap = 2500, ap_radius_physical = 10*u.arcsec, version = 'v1'):
    import random
    from matplotlib import pyplot as plt
    from matplotlib import gridspec
    import seaborn as sns
    import DataReduction as DataRed
    
    galaxy_name = galaxy_properties.galaxy_name
    subprocess.call('mkdir ../'+galaxy_name+'/maps_rms_evaluation', shell = True)
    
    if os.path.exists('maps_rms_evaluation/map_rms.txt'): return
    else: pass
    
    # ---------------------------------------------------------------------------
    # Galaxy Aperture Stuff, from Dustpedia (to mask and bkg evaluation purposes)
    DustPedia_Photom = pd.read_csv('../DustPedia_Tables/DustPedia_Aperture_Photometry_2.2.csv')
    subtable = DustPedia_Photom.loc[DustPedia_Photom['name'] == galaxy_name]
    ra, dec = subtable['ra'].values[0], subtable['dec'].values[0]
    ap_cen_coord = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    semimaj, axial_ratio, pos_angle = \
        subtable['semimaj_arcsec'].values[0], subtable['axial_ratio'].values[0], subtable['pos_angle'].values[0]
    # ---------------------------------------------------------------------------

    # ---------
    # READ MAPS
    print()
    print('Reading reduced maps...')
    list_data, list_band, list_wvl = [], [], []
    path_fits = '../'+galaxy_name+'/_ReducedMaps/'
    for file in os.listdir(path_fits):
        if not file.endswith('.fits'): continue
        elif file.startswith('InRGB.fits'): continue
        temp = GUS.FitsUtils(path_fits+'/'+file)
        list_data.append(temp)
        list_band.append(temp.bandname)
        list_wvl.append(temp.get_wavelength())
    
    list_wvl, list_band, list_data = (t for t in zip(*sorted(zip(list_wvl, list_band, list_data))))   
    print('...done!')
    print()
    # ---------

    list_map_rms = []
    for data in list_data:
                
        # 1) Read data, get pixel scale
        bandname = data.bandname
        
        if os.path.exists('maps_rms_evaluation/'+galaxy_name+'_'+bandname+'_rms_histogram.pdf'):
            print(space+'rms already evaluated for', data.bandname)
            continue

        print(space+'rms evaluation for', data.bandname)
        image = data.signal
        ok = np.where(image == 0)
        image[ok] = np.nan # This, to avoid that zeros bias the rms evaluation
                        
        pixel_scale = (data.get_pixel_scale()*u.deg).to('arcsec').value
        nx, ny = data.hdr['NAXIS1'], data.hdr['NAXIS2']
        range_x, range_y = [0 + nx/10, nx - nx/10], [0 + ny/10, ny - ny/10]
    
        # Galaxy Aperture Stuff, from Dustpedia (to mask galaxy emission for rms evaluation)
        centre_x, centre_y = ap_cen_coord.to_pixel(data.wcs)
        Gal_Ap_Stuff = centre_x, centre_y, semimaj/pixel_scale, axial_ratio, pos_angle      
                    
        # 2) Sigma clip the image
        clip_value = GUS.SigmaClip(image, tolerance=0.01, sigma_thresh=3.0, median=True)
        noise_value = clip_value[0]
        field_value = clip_value[1]
        cutoff_sigma = 2.0
        cutoff = field_value + ( cutoff_sigma * noise_value )
        
        # 3) Mask the image removing galaxy emission...
        image_masked = image.copy()
        centre_i, centre_j, mask_semimaj_pix, mask_axial_ratio, mask_angle = Gal_Ap_Stuff
        ellipse_mask = DataRed.EllipseMask(image, mask_semimaj_pix, mask_axial_ratio, mask_angle, centre_i, centre_j)
        image_masked[ np.where( ellipse_mask==1 ) ] = np.nan
        # ...and image pixels identified as having high SNR
        image_masked[ np.where( image_masked>cutoff ) ] = np.nan
        # ...and image pixels equal to zero, because Spitzer sucks
        image_masked[ np.where( image_masked==0 ) ] = np.nan
                
        # 4) Evaluate the rms in N_ap apertures, depending on ap_radius (in pixels) and image size (later)
        fig = plt.figure(figsize=(20,8))
        gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1]) 
        ax1 = plt.subplot(gs[0], projection = data.wcs)
        ax2 = plt.subplot(gs[1])
        ax1.set_title('Masked image')
        ax1.imshow(image_masked, origin = 'lower', interpolation = 'nearest', cmap = GUS.associate_colormap(bandname))

        if version == 'v1':
            ap_radius = (ap_radius_physical/pixel_scale).value
            rms_list = []
            for i in range(N_ap):
                i_centre = int(random.uniform(range_x[0], range_x[1]))
                j_centre = int(random.uniform(range_y[0], range_y[1]))
                try: check_val = np.nan_to_num(image_masked[i_centre, j_centre])
                except: continue
                if check_val == 0: continue
                #ax2.plot(i_centre, j_centre, 'r.', ms = 1.0)
                circle = plt.Circle((j_centre, i_centre), ap_radius, color='g', lw = 1.0, fill=False)
                ax1.add_artist(circle)
                ap_sum = DataRed.CircleSum(image_masked, i_centre, j_centre, ap_radius)
                pixels_in_ap = np.array(ap_sum[2])
                x = pixels_in_ap[~np.isnan(pixels_in_ap)]
                rms_list.append(np.sqrt(np.mean(x**2)))
            
            rms_list = np.array(rms_list)
            rms_list = rms_list[~np.isnan(rms_list)]
    
            n, b, patches = ax2.hist(rms_list, 50, density = True, histtype='stepfilled')
            map_rms = b[np.argmax(n)]
            list_map_rms.append(map_rms)
            ax2.set_title('Map rms = {0:.2E}'.format(map_rms))
            fig.savefig('maps_rms_evaluation/'+galaxy_name+'_'+bandname+'_rms_histogram.pdf', bbox_inches = 'tight')
            print(space+space+'rms is {0:.2E}'.format(map_rms))
            
        elif version == 'v2':
            ap_radius = (ap_radius_physical/pixel_scale).value
            fluxes_list = []
            for i in range(N_ap):
                i_centre = int(random.uniform(range_x[0], range_x[1]))
                j_centre = int(random.uniform(range_y[0], range_y[1]))
                try: check_val = np.nan_to_num(image_masked[i_centre, j_centre])
                except: continue
                if check_val == 0: continue
                circle = plt.Circle((j_centre, i_centre), ap_radius, color = 'g', lw = 1.0, fill=False)
                ax1.add_artist(circle)
                ap_sum = DataRed.CircleSum(image_masked, i_centre, j_centre, ap_radius)
                pixels_in_ap = np.array(ap_sum[2])
                x = pixels_in_ap[~np.isnan(pixels_in_ap)]
                fluxes_list.append(np.sum(x)/len(x))
            
            fluxes_list = np.array(fluxes_list)
            fluxes_list = fluxes_list[~np.isnan(fluxes_list)] 
            map_rms = np.std(fluxes_list)
            
            rms_hist = sns.distplot(fluxes_list, bins = 15, kde=True, color = 'tab:blue', ax = ax2)
            #rms_hist.set(xlabel=r'${\rm Fluxes}$', ylabel=' ')
            ax2.set_title('Map rms = {0:.2E}'.format(map_rms))
            fig.savefig('maps_rms_evaluation/'+galaxy_name+'_'+bandname+'_rms_histogram.pdf', bbox_inches = 'tight')
            print(space+space+'rms is {0:.2E}'.format(map_rms))
         
    ascii.write([list_band, list_map_rms], 'maps_rms_evaluation/map_rms.txt', names = ['Band', 'rms'], overwrite = True)
    
    return

def generate_coordinate_grid(galaxy_properties, fits_base_band, cell_side, run_type, avoidance_radius):
    
    galaxy_name = galaxy_properties.galaxy_name
    
    AngPhyConv = GUS.AngularPhysicalConv(galaxy_properties.z_source)
    
    subprocess.call('mkdir ../'+galaxy_name+'/Photometry_'+run_type, shell = True)
    folder_path = '../'+galaxy_name+'/Photometry_'+run_type
    
    fits_base = GUS.FitsUtils('../'+galaxy_name+'/_ReducedMaps/'+fits_base_band+'.fits')
    nx, ny = fits_base.hdr['NAXIS1'], fits_base.hdr['NAXIS2']
    pixel_scale = fits_base.get_pixel_scale()*u.deg.to('arcsec')
    
    cen_pos = SkyCoord(galaxy_properties.ra, galaxy_properties.dec, frame = 'icrs')
    cen_xpix, cen_ypix = skycoord_to_pixel(cen_pos, fits_base.wcs)
    
    try: unit = avoidance_radius.unit
    except: return 'You should give a unit to the avoidance radius.'
    if unit == 'deg' or unit == 'arcsec': avoidance_radius = avoidance_radius.to('arcsec')
    else: avoidance_radius = AngPhyConv.P2A(avoidance_radius).to('arcsec')
    
    avoidance_radius_pixel = avoidance_radius.value/pixel_scale
    x_avoid_min, x_avoid_max = cen_xpix-avoidance_radius_pixel, cen_xpix+avoidance_radius_pixel
    y_avoid_min, y_avoid_max = cen_ypix-avoidance_radius_pixel, cen_ypix+avoidance_radius_pixel
        
    if cell_side.unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_side.unit == 'arcmin': cell_side = cell_side.to('arcsec')
        
    try: unit = cell_side.unit
    except: return 'You should give a unit to the cell size.'
    
    if unit == 'arcsec' or unit == 'arcmin' or unit == 'deg':
        print('You want cells of {0:.3f} side'.format(cell_side))
        cell_in_kpc = AngPhyConv.A2P(cell_side).to('kpc')
        print('This means you will sample scales of {0:.2f}'.format(cell_in_kpc))
    if unit == 'kpc' or unit == 'pc':
        cell_in_kpc = cell_side.to('kpc')
        print('You want cells of '+str(cell_side)+' side')
        cell_side = AngPhyConv.P2A(cell_side).to('arcsec')
        print('At z = {0:.6f} this means cells of {1:.3f} side'.format(galaxy_properties.z_source, cell_side)) 
    if cell_in_kpc < 0.2*u.kpc:
        print("Warning! You're asking me to sample scales below the typical ones for GMCs (approx 200 pc).")
        print("Please, stop. Get some help.")
        
    cell_side_pixel = cell_side.value/pixel_scale
    cell_radius = cell_side/2
    cell_radius_pixel = cell_radius.value/pixel_scale
    number_of_squares = nx/cell_radius_pixel
    x_pos, y_pos = np.arange(0, nx, cell_side_pixel), np.arange(0, ny, cell_side_pixel)
    
    list_ra, list_dec = [], []
    for ix in tqdm(x_pos, ascii = True, desc = 'Generating coordinates grid'):
        range_x = int(cell_radius_pixel)
        if (ix+range_x) >= nx: continue
        elif (ix-range_x) <= 0: continue
        elif (ix < x_avoid_min): continue
        elif (ix > x_avoid_max): continue
        
        for iy in y_pos:
            range_y = int(cell_radius_pixel)      
            if (iy+range_y) >= nx: continue
            elif (iy-range_y) <= 0: continue
            elif (iy < y_avoid_min): continue
            elif (iy > y_avoid_max): continue
                    
            coord = pixel_to_skycoord(ix, iy, fits_base.wcs)
            list_ra.append(coord.ra.value), list_dec.append(coord.dec.value)

    ascii.write([list_ra, list_dec], folder_path+'/coordinates.txt', names = ['ra', 'dec'])
    print('Coordinates grid saved in '+folder_path)
    return

def generate_coordinate_grid_within_DustAp(galaxy_properties, fits_base_band, cell_side, run_type, avoidance_radius, \
                                           avoidance_semimaj = False, avoidance_angle = False, avoidance_axial_ratio = False):
    import matplotlib.pyplot as plt
    import DataReduction as DataRed

    galaxy_name = galaxy_properties.galaxy_name
    AngPhyConv = GUS.AngularPhysicalConv(galaxy_properties.z_source)
      
    subprocess.call('mkdir ../'+galaxy_name+'/Photometry_'+run_type, shell = True)
    folder_path = '../'+galaxy_name+'/Photometry_'+run_type

    fits_base = GUS.FitsUtils('../'+galaxy_name+'/_ReducedMaps/'+fits_base_band+'.fits')
    nx, ny = fits_base.hdr['NAXIS1'], fits_base.hdr['NAXIS2']
    pixel_scale = fits_base.get_pixel_scale()*u.deg.to('arcsec')
    
    # Galaxy Aperture Stuff, from Dustpedia
    DustPedia_Photom = pd.read_csv('../DustPedia_Tables/DustPedia_Aperture_Photometry_2.2.csv')
    subtable = DustPedia_Photom.loc[DustPedia_Photom['name'] == galaxy_name]
    ap_cen_coord = SkyCoord(subtable['ra'].values[0]*u.deg, subtable['dec'].values[0]*u.deg, frame = 'fk5')
    cen_xpix, cen_ypix = ap_cen_coord.to_pixel(fits_base.wcs)
    if avoidance_semimaj == False: mask_semimaj_pix, mask_axial_ratio, mask_angle = \
        subtable['semimaj_arcsec']/pixel_scale, subtable['axial_ratio'].values[0], subtable['pos_angle'].values[0]
    else: mask_semimaj_pix, mask_axial_ratio, mask_angle = \
        avoidance_semimaj/pixel_scale, subtable['axial_ratio'].values[0], subtable['pos_angle'].values[0]
    try: mask_semimaj_pix = mask_semimaj_pix.value
    except: pass
    if avoidance_angle == False: pass 
    else: mask_angle = np.copy(avoidance_angle)
    if avoidance_axial_ratio == False: pass 
    else: mask_axial_ratio = np.copy(avoidance_axial_ratio)
    
    image_masked = fits_base.signal.copy()
    ellipse_mask = DataRed.EllipseMask(fits_base.signal, mask_semimaj_pix, mask_axial_ratio, mask_angle, cen_xpix, cen_ypix)
    image_masked[ np.where( ellipse_mask != 1 ) ] = -99
        
    try: unit = avoidance_radius.unit
    except: return 'You should give a unit to the avoidance radius.'
    if unit == 'deg' or unit == 'arcsec': avoidance_radius = avoidance_radius.to('arcsec')
    else: avoidance_radius = AngPhyConv.P2A(avoidance_radius).to('arcsec')
    
    avoidance_radius_pixel = avoidance_radius.value/pixel_scale
    x_avoid_min, x_avoid_max = cen_xpix-avoidance_radius_pixel, cen_xpix+avoidance_radius_pixel
    y_avoid_min, y_avoid_max = cen_ypix-avoidance_radius_pixel, cen_ypix+avoidance_radius_pixel
        
    if cell_side.unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_side.unit == 'arcmin': cell_side = cell_side.to('arcsec')
    try: unit = cell_side.unit
    except: return 'You should give a unit to the cell size.'
    if unit == 'arcsec' or unit == 'arcmin' or unit == 'deg':
        print("You want cells of {0:.3f} side".format(cell_side))
        cell_in_kpc = AngPhyConv.A2P(cell_side).to('kpc')
        print("This means you will sample scales of {0:.2f}".format(cell_in_kpc))
    if unit == 'kpc' or unit == 'pc':
        cell_in_kpc = np.copy(cell_side)
        print("You want cells of "+str(cell_side)+" side")
        cell_side = AngPhyConv.P2A(cell_side).to('arcsec')
        print("At z = {0:.6f} this means cells of {1:.3f} side".format(galaxy_properties.z_source, cell_side))  
    if cell_in_kpc < 0.2*u.kpc:
        print("Warning! You're asking me to sample scales below the typical ones for GMCs (approx 200 pc).")
        print("Please, stop. Get some help.")
        
    cell_side_pixel = cell_side.value/pixel_scale
    cell_radius = cell_side/2
    cell_radius_pixel = cell_radius.value/pixel_scale
    number_of_squares = nx/cell_radius_pixel
    x_pos, y_pos = np.arange(0, nx, cell_side_pixel), np.arange(0, ny, cell_side_pixel)
            
    list_ra, list_dec = [], []
    for ix in tqdm(x_pos, ascii = True, desc = 'Generating coordinates grid'):
        range_x = int(cell_radius_pixel)
        if (ix+range_x) >= nx: continue
        elif (ix-range_x) <= 0: continue
        elif (ix < x_avoid_min): continue
        elif (ix > x_avoid_max): continue
        
        for iy in y_pos:
            if image_masked[int(ix), int(iy)] == -99: continue
            range_y = int(cell_radius_pixel)      
            if (iy+range_y) >= nx: continue
            elif (iy-range_y) <= 0: continue
            elif (iy < y_avoid_min): continue
            elif (iy > y_avoid_max): continue
                    
            coord = pixel_to_skycoord(ix, iy, fits_base.wcs)
            list_ra.append(coord.ra.value), list_dec.append(coord.dec.value)

    ascii.write([list_ra, list_dec], folder_path+'/coordinates.txt', names = ['ra', 'dec'])
    return

    
def check_coordinate_grid(galaxy_properties, fits_base_band, run_type, plot_size):
    import matplotlib
    matplotlib.rcParams['backend'] = "Qt4Agg"
    import matplotlib.pyplot as plt
    from matplotlib import cm, gridspec

    galaxy_name = galaxy_properties.galaxy_name
    AngPhyConv = GUS.AngularPhysicalConv(galaxy_properties.z_source)
    fits_base = GUS.FitsUtils('../'+galaxy_name+'/_ReducedMaps/'+fits_base_band+'.fits')    
    pixel_scale = fits_base.get_pixel_scale()*u.deg.to('arcsec')
    cen_pos = SkyCoord(galaxy_properties.ra, galaxy_properties.dec, frame = 'icrs')
    cen_xpix, cen_ypix = skycoord_to_pixel(cen_pos, fits_base.wcs)
    
    # ---------
    # Plot size
    try: unit = plot_size.unit
    except: raise Exception('You must give a unit to the plot size.')
    if unit == 'arcsec' or unit == 'arcmin' or unit == 'deg':
        plot_size_in_kpc = 2*AngPhyConv.A2P(plot_size).to('kpc')
        print('Plot side will be of {0:.2f}'.format(plot_size_in_kpc))
    if unit == 'kpc' or unit == 'pc':
        plot_size = AngPhyConv.P2A(plot_size).to('arcsec')
        print('At z = {0:.6f} plot side will be of {1:.1f}'.format(galaxy_properties.z_source, 2*plot_size)) 
    plot_size_pixel = (plot_size/pixel_scale).value
    x_min, x_max = cen_xpix-plot_size_pixel, cen_xpix+plot_size_pixel
    y_min, y_max = cen_ypix-plot_size_pixel, cen_ypix+plot_size_pixel
    nx, ny = fits_base.hdr['NAXIS1'], fits_base.hdr['NAXIS2']
    pixel_scale = fits_base.get_pixel_scale()*u.deg.to('arcsec')
    # ---------

    data = ascii.read('../'+galaxy_name+'/Photometry_'+run_type+'/coordinates.txt')
    ra_apertures, dec_apertures = data['ra'].data*u.deg, data['dec'].data*u.deg
    coords = SkyCoord(ra_apertures, dec_apertures, frame = 'icrs')
    
    fig = plt.figure(figsize=(15,15))
    ax = plt.subplot(1,1,1)
    ax.imshow(fits_base.signal, origin = 'lower', interpolation = 'nearest', cmap = GUS.associate_colormap(fits_base.bandname))
    for coord in tqdm(coords, ascii = True, desc = 'Plotting coordinates'):
        ra_x, dec_y = coord.to_pixel(fits_base.wcs)[0], coord.to_pixel(fits_base.wcs)[1]
        ax.plot(ra_x, dec_y, marker = '*', markersize = 5, color = 'indianred')
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    fig.savefig('../'+galaxy_name+'/Photometry_'+run_type+'/grid.pdf', bbox_inches = 'tight')
    return

def read_apertures(path, cell_side):
    from photutils import SkyRectangularAperture
    ra_apertures, dec_apertures = ascii.read(path)['ra'].data*u.deg, ascii.read(path)['dec'].data*u.deg
    cell_radius = cell_side/2
    coords = SkyCoord(ra_apertures, dec_apertures, frame = 'icrs')
    return SkyRectangularAperture(coords, w=cell_side, h=cell_side, theta = 0*u.deg)

def read_calibration_dictionary(bandname):
    '''
    Band calibration uncertainties.
    INPUT: the band name
    OUTPUT: the calibration uncertainty, as a percentage of the measured flux
    WRITTEN: A. Enia, Dec 2019
    '''
    calib_dictionary = dict(zip(ascii.read(DustPath+'scripts/calibration_uncertainties.txt')['BandName'].data, \
                          ascii.read(DustPath+'scripts/calibration_uncertainties.txt')['CalibUnc'].data))
    return calib_dictionary[GUS.band_disambiguation(bandname)]

def remove_weirdities(table_path, threshold = 1):

    print()
    print('You selected the remove weird data points option.')
    print('Congratulations!')
    print('This means that weird bumps in photometries will be removed by placing -99 to the bump value.')
    print()
    magphys_table = pd.read_csv(table_path, sep = '\t', dtype={'id': object})
    index_photometry = np.arange(start = 2, stop = 22, step = 2) # 22 perche' e' 2MASS_Ks
    
    count_modif = 0
    modified_table = []
    for i_row, row in magphys_table.iterrows():
        # Condizione 1: sulla fotometria UV-to-nearIR, se c'e' un singolo punto a ramengo rispetto a due attigui, va scartato.
        #               A ramengo vuol dire che la distanza in log tra i punti attigui è maggiore di uno.
        photom = row.values[index_photometry]
        photom_iter = range(1, 8)
        for i_ph in photom_iter:
            if photom[i_ph] == -99.0 or photom[i_ph-1] == -99.0 or photom[i_ph+1] == -99.0: continue   
            diff1 = np.abs(np.log10(photom[i_ph]) - np.log10(photom[i_ph-1]))
            diff2 = np.abs(np.log10(photom[i_ph+1]) - np.log10(photom[i_ph]))
            if diff1 > threshold and diff2 > threshold:
                #print('Datapoint', i_ph, 'of row', i_row, 'xe marso')
                magphys_table.iat[i_row,i_ph*2+2] = -99.0
                count_modif += 1
    
    print('Total number of modified values', count_modif)
    print()
    magphys_table.to_csv(table_path, sep = '\t', index = None)
    
    return

def do_photometry_build_table(galaxy_properties, working_bands, cell_side, run_type, \
                  SNR_threshold = 1.5, skip_value = 5, remove_weird_points = True):
    '''
    Perform photometry and build the MAGPHYS fluxes table, with errors coming from the measured root mean square of the maps.
    INPUT: galaxy properties, list of band, cell side (either in arcsec or kpc), run type.
    OUTPUT: it builds three files, 'photometries_table.csv', FINAL.csv' and 'Final_for_MAGPHYS.csv' in the photometry folder.
    WRITTEN: A. Enia, Dec 2019
    MODIFIED: A. Enia, Dec 2020
    '''
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    import matplotlib
    matplotlib.rcParams['backend'] = "Qt4Agg"
    import matplotlib.pyplot as plt
    plt.close('all') # Clear all the previously opened graphic windows
    from matplotlib import cm, gridspec
    import DataReduction as DataRed

    # =================
    # GALAXY PROPERTIES
    galaxy_name = galaxy_properties.galaxy_name
    AngPhyConv = GUS.AngularPhysicalConv(galaxy_properties.z_source)
    D_L = galaxy_properties.dist
    source_z = np.round(galaxy_properties.z_source, 6)
    ra_cen, dec_cen = galaxy_properties.ra, galaxy_properties.dec
    # =================
    
    # ================
    # SQUARE APERTURES
    photometry_folder_path = '../'+galaxy_name+'/Photometry_'+run_type+'/'
    coordinates_path = photometry_folder_path+'coordinates.txt'
    try: cell_unit = cell_side.unit
    except: raise ValueError('You must give a unit to cell side.')
    if cell_unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_unit == 'arcmin' or cell_unit == 'deg': cell_side = cell_side.to('arcsec')
    if cell_side.unit == 'arcsec':
        apertures = read_apertures(coordinates_path, cell_side)
    elif cell_side.unit == 'kpc':
        cell_side = AngPhyConv.P2A(cell_side).to('arcsec')
        apertures = read_apertures(coordinates_path, cell_side)
    # ================

    # =========
    # READ MAPS
    list_data, list_band, list_wvl, list_pixel_scales = [], [], [], []
    path_fits = '../'+galaxy_name+'/_ReducedMaps/'
    for file in os.listdir(path_fits):
        if not file.endswith('.fits'): continue
        elif file.startswith('InRGB.fits') or file.startswith('_'): continue
        temp = GUS.FitsUtils(path_fits+'/'+file)
        list_data.append(temp)
        list_band.append(temp.bandname)
        list_wvl.append(temp.get_wavelength())
        list_pixel_scales.append(temp.get_pixel_scale()*u.deg.to('arcsec'))
    
    list_wvl, list_band, list_data, list_pixel_scales = (t for t in zip(*sorted(zip(list_wvl, list_band, list_data, list_pixel_scales))))   
    # =========

    # ==================
    # PERFORM PHOTOMETRY
    # Read dictionaries:
    #   - galactic extinction correction
    GalCorr_path = '../'+galaxy_name+'/galactic_extinction_correction.txt'
    if os.path.exists(GalCorr_path): pass
    else: DataRed.GalExtCorr(galaxy_name, working_bands, galaxy_properties.ra, galaxy_properties.dec)
    GalExtCorr_Dict = dict(zip(ascii.read(GalCorr_path)['Band'].data, \
                      ascii.read(GalCorr_path)['Correction'].data))
    #   - rms table. rms is per pixel. Pixels in the aperture vary between bands, so multiply accordingly.
    try: rms_table = ascii.read('../'+galaxy_name+'/maps_rms_evaluation/map_rms.txt')
    except: raise ReadError('You must estimate the maps rms before building the MAGPHYS table.')
    pixels_per_aperture = ((cell_side/list_pixel_scales)**2).value
    RmsPerAperture_Dict = dict(zip(rms_table['Band'].data, pixels_per_aperture*GUS.round_arr(rms_table['rms'].data, 4)))            
    #   - calibration error dictionary
    Calib_Dict = dict(zip(list_band, np.array([read_calibration_dictionary(band) for band in list_band])))
    
    # Ok, now we can do photometry.
    subprocess.call('mkdir '+photometry_folder_path+'/fluxes', shell = True)
    magphys_table = pd.DataFrame()
    magphys_table['ra'] = ascii.read(coordinates_path)['ra']
    magphys_table['dec'] = ascii.read(coordinates_path)['dec']
    from photutils import aperture_photometry
    for Map in tqdm(list_data, ascii = True, desc = 'Performing photometry on maps'):
        phot_table = aperture_photometry(Map.signal, apertures, wcs = Map.wcs, method = 'subpixel', subpixels = 32)
        phot_table['aperture_sum'].info.format = '%.2g'
        phot = phot_table['aperture_sum'].data
        phot = np.nan_to_num(phot) # Eh...
        try:
            err_table = aperture_photometry(Map.errormap, apertures, wcs = Map.wcs, method = 'subpixel', subpixels = 32)
            phot_err = err_table['aperture_sum'].data
        except: 
            phot_err = RmsPerAperture_Dict[Map.bandname] + phot*(Calib_Dict[Map.bandname]/100)
        # SNR threshold
        BAD_SNR = phot/phot_err < SNR_threshold
        phot[BAD_SNR] = -99
        phot_err[BAD_SNR] = -99
        # Save to table
        magphys_table[Map.bandname] = phot
        magphys_table[Map.bandname+'_ERROR'] = phot_err

    magphys_table.to_csv(photometry_folder_path+'photometries_table.csv', index = False, sep = '\t')
    # ==================
    
    print('Photometry phase over.')
    print('Building MAGPHYS table')
        
    # ==================
    # Remove rows with more than skip_value negative fluxes           
    count_df = 0*magphys_table.filter(regex = '^((?!err).)*$')
    count_df[magphys_table < 0] = 1
    GOOD_rows = count_df.sum(axis = 1) <= skip_value
    BAD_rows = count_df.sum(axis = 1) > skip_value
    magphys_table_GOOD = magphys_table[GOOD_rows].reset_index().drop(columns = ['index'])
    # And put negative numbers to -99
    negval_cond = magphys_table_GOOD <= 0
    negval_cond['ra'] = False
    negval_cond['dec'] = False
    magphys_table_GOOD[negval_cond] = -99 ############## !!!!!!!!!!! ################# 
    # ==================
    
    # ================================
    good_apertures = apertures[GOOD_rows]
    number_of_good_apertures = len(good_apertures)
    print('Good apertures: ', number_of_good_apertures)
    z_fill_number = int(np.trunc(np.log10(number_of_good_apertures)))+2
    magphys_table_GOOD['id'] = [idx.zfill(z_fill_number) for idx in magphys_table_GOOD.index.to_numpy().astype(str)]
    magphys_table_GOOD.set_index('id')
    # ================================

    # ==================
    # BUILD THE TABLE(S)
    magphys_table_GOOD['z'] = source_z
    magphys_table_GOOD = magphys_table_GOOD[['z'] + magphys_table_GOOD.columns[:-1].tolist()]
    
    z_fill_number = int(np.trunc(np.log10(number_of_good_apertures)))+2
    magphys_table_GOOD['id'] = [idx.zfill(z_fill_number) for idx in magphys_table_GOOD.index.to_numpy().astype(str)]
    magphys_table_GOOD.set_index('id')
    magphys_table_GOOD = magphys_table_GOOD[['id'] + magphys_table_GOOD.columns[:-1].tolist()]
    magphys_table_GOOD.to_csv(photometry_folder_path+'/GOOD_photometries_table.csv', index=False, sep='\t')
    magphys_table_FINAL = magphys_table_GOOD.drop(columns = ['ra', 'dec'])
    magphys_table_FINAL.to_csv(photometry_folder_path+'/Final_for_MAGPHYS.csv', index=False, sep='\t')
    # ==================
    
    if remove_weird_points: remove_weirdities(photometry_folder_path+'/Final_for_MAGPHYS.csv')
    
    return

########################################################################################################################

def photometry_on_aperture_v0(galaxy_name, list_data, aperture, GalExtCorr_dict, skip_value):
    from photutils import SkyRectangularAperture, aperture_photometry
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) 

    list_fluxes, list_flux_errors = [], []
    i_check = 0
    for obs_data in list_data:
        # Perform photometry
        phot_table = aperture_photometry(obs_data.signal, aperture, wcs = obs_data.wcs, method = 'subpixel', subpixels = 32)
        phot_table['aperture_sum'].info.format = '%.4g'  
        # Put results in a single file
        phot = GUS.round_arr(phot_table['aperture_sum'].data, 2)
        phot = np.nan_to_num(phot) # Eh...
        try:
            err_table = aperture_photometry(obs_data.error, aperture, wcs = obs_data.wcs, method = 'subpixel', subpixels = 32)
            err_table['aperture_sum'].info.format = '%.4g'  
            phot_err = GUS.round_arr(err_table['aperture_sum'].data, 2)
        except: 
            phot_err = np.array([0.0])  
        # If a single phot value is < 0, assign -99, magphys will later ignore it    
        if phot <= 0:
            list_fluxes.append(np.array([-99]))
            list_flux_errors.append(np.array([-99]))
            i_check += 1
            continue      
        # If more than "skip_value" band values are shit, give -99 to them all
        if i_check >= skip_value: 
            list_fluxes = np.zeros(len(list_data)) - 99.0
            list_flux_errors = np.zeros(len(list_data)) - 99.0 
            break           
        # Galactic extintion correction
        try:
            phot *= GalExtCorr_dict[obs_data.band]
            phot_err *= GalExtCorr_dict[obs_data.band]
        except: pass
        list_fluxes.append(phot), list_flux_errors.append(phot_err)
    # Sort and save
    fluxes, errors = np.array(np.squeeze(list_fluxes)), np.array(np.squeeze(list_flux_errors))
    return fluxes, errors

def do_photometry_v0(galaxy_properties, working_bands, cell_side, run_type, skip_value = 5, erase = False):
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    import matplotlib
    matplotlib.rcParams['backend'] = "Qt4Agg"
    import matplotlib.pyplot as plt
    plt.close('all') # Clear all the previously opened graphic windows
    from matplotlib import cm, gridspec
    Wvlghts_dictionary = dict(zip(ascii.read('../sample_properties/bands_and_wvlghts.txt')['name'].data, \
                          ascii.read('../sample_properties/bands_and_wvlghts.txt')['lambda_eff'].data))
  
    # =================
    # GALAXY PROPERTIES
    galaxy_name = galaxy_properties.galaxy_name
    AngPhyConv = GUS.AngularPhysicalConv(galaxy_properties.z_source)
    D_L = galaxy_properties.dist
    source_z = galaxy_properties.dist*cosmo.H0/v_lux.to('km/s')
    ra_cen, dec_cen = galaxy_properties.ra, galaxy_properties.dec
    # =================
    
    # ================
    # SQUARE APERTURES
    photometry_folder_path = '../'+galaxy_name+'/Photometry_'+run_type+'/'
    coordinates_path = photometry_folder_path+'coordinates.txt'
    try: cell_unit = cell_side.unit
    except: raise ValueError('You must give a unit to cell side.')
    
    if cell_unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_unit == 'arcmin' or cell_unit == 'deg': cell_side = cell_side.to('arcsec')
        
    if cell_side.unit == 'arcsec':
        list_apertures = read_apertures(coordinates_path, cell_side)
    elif cell_side.unit == 'kpc':
        cell_side = AngPhyConv.P2A(cell_side).to('arcsec')
        list_apertures = read_apertures(coordinates_path, cell_side)
    # ================

    # =========
    # READ MAPS
    list_data, list_band, list_wvl = [], [], []
    path_fits = '../'+galaxy_name+'/_ReducedMaps/'
    for file in os.listdir(path_fits):
        if not file.endswith('.fits'): continue
        elif file.startswith('InRGB.fits'): continue
        temp = GUS.FitsUtils(path_fits+'/'+file)
        list_data.append(temp)
        list_band.append(temp.bandname)
        list_wvl.append(temp.get_wavelength())
    
    list_wvl, list_band, list_data = (t for t in zip(*sorted(zip(list_wvl, list_band, list_data))))   
    # =========

    # ==================
    # PERFORM PHOTOMETRY
    GalCorr_path = '../'+galaxy_name+'/galactic_extinction_correction.txt'
    if os.path.exists(GalCorr_path): pass
    else: GalExtCorr(galaxy_name, working_bands, ra, dec)
    GalCorrection_dictionary = dict(zip(ascii.read(GalCorr_path)['Band'].data, \
                      ascii.read(GalCorr_path)['Correction'].data))
    
    subprocess.call('mkdir '+photometry_folder_path+'/fluxes', shell = True)
    if erase:
        subprocess.call('rm -rf '+photometry_folder_path+'/fluxes/ap*', shell = True)
    print()
    for aperture, i_ap in tqdm(zip(list_apertures, range(len(list_apertures))), ascii = True, desc = 'Performing photometry on the apertures'):
        # Perform photometry
        list_fluxes, list_flux_errors = photometry_on_aperture(galaxy_name, list_data, aperture, GalExtCorr_dict, skip_value)
        # Sort and save
        fluxes, errors = np.array(np.squeeze(list_fluxes)), np.array(np.squeeze(list_flux_errors))
        ascii.write([list_band, GUS.round_arr(np.array(list_wvl),2), GUS.round_arr(fluxes, 2), GUS.round_arr(errors, 2)], \
                photometry_folder_path+'fluxes/ap_{0:004}_fluxes.txt'.format(i_ap), \
                names = ['Band', 'Wvl', 'Fluxes', 'Errors'], overwrite = True)
        
    print('Photometry phase over.')
    print()
    # ==================
    return

def build_table_rms_v0(galaxy_properties, list_band, cell_side, run_type, skip_value = 5, remove_weird_points = False):
    '''
    Building the MAGPHYS fluxes table, with errors coming from the measured root mean square of the maps.
    INPUT: galaxy properties, list of band, cell side (either in arcsec or kpc), run type. Optional: skip value (the one over which ignore the aperture)
            and the option to remove weird points in the SED.
    OUTPUT: it builds two files, 'FINAL.csv' and 'Final_for_MAGPHYS.csv' in the photometry folder.
    WRITTEN: A. Enia, Dec 2019
    '''
        
    # =================
    # GALAXY PROPERTIES
    galaxy_name = galaxy_properties.galaxy_name
    D_L = galaxy_properties.dist
    source_z = galaxy_properties.dist*cosmo.H0/v_lux.to('km/s')
    source_z = np.round(source_z, 6)
    AngPhyConv = GUS.AngularPhysicalConv(galaxy_properties.z_source)
    ra_cen, dec_cen = galaxy_properties.ra, galaxy_properties.dec

    list_wvl, list_pixel_scales = [], []
    path_fits = '../'+galaxy_name+'/_ReducedMaps/'
    for file in os.listdir(path_fits):
        if not file.endswith('.fits'): continue
        elif file.startswith('InRGB.fits'): continue
        temp = GUS.FitsUtils(path_fits+'/'+file)
        list_pixel_scales.append(temp.get_pixel_scale()*u.deg.to('arcsec'))
        list_wvl.append(temp.get_wavelength())    
    list_wvl, list_pixel_scales = (t for t in zip(*sorted(zip(list_wvl, list_pixel_scales))))   
    # =================

    # =================
    # Read the apertures
    try: cell_unit = cell_side.unit
    except: raise ValueError('You must give a unit to cell side.')
    if cell_unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_unit == 'arcmin' or cell_unit == 'deg': cell_side = cell_side.to('arcsec')
    photometry_folder_path = '../'+galaxy_name+'/Photometry_'+run_type+'/'
    coordinates_path = photometry_folder_path+'coordinates.txt'
    if cell_side.unit == 'kpc': cell_side = AngPhyConv.P2A(cell_side).to('arcsec')
    list_apertures = read_apertures(coordinates_path, cell_side)
    # =================
           
    # ==================
    try: rms_table = ascii.read('../'+galaxy_name+'/maps_rms_evaluation/map_rms.txt')
    except: raise ReadError('You must estimate the maps rms before building the MAGPHYS table.')
    rms_bands = rms_table['Band'].data
    rms_per_band = GUS.round_arr(rms_table['rms'].data, 4)
    rms_dictionary = dict(zip(rms_bands, rms_per_band))
    calibration_uncertainties = np.array([read_calibration_dictionary(band) for band in list_band])
            
    # RMS IS PER PIXEL. THE PIXELS IN THE APERTURE VARY WITH THE BANDS
    pixels_per_aperture = ((cell_side/list_pixel_scales)**2).value
    
    final_fluxes, final_errors, apertures_mask = [], [], []
    for i_ap in trange(len(list_apertures), ascii = True, desc = 'Reading photometries'):
        data = ascii.read(photometry_folder_path+'/fluxes/ap_{0:004}_fluxes.txt'.format(i_ap))
        try: single_ap_flux, single_ap_error = data['Fluxes'].data, data['Errors'].data
        except: single_ap_flux, single_ap_error = data['Fluxes'].data, 0*data['Fluxes'].data
        check = np.count_nonzero(single_ap_flux == -99.0)
        if check <= skip_value:
            bad_err = np.where(single_ap_error == 0)
            single_ap_error[bad_err] = (pixels_per_aperture*rms_per_band)[bad_err] + (single_ap_flux*(calibration_uncertainties/100))[bad_err]
            final_fluxes.append(single_ap_flux), final_errors.append(single_ap_error), apertures_mask.append(i_ap)
        else: pass
    
    errors_matrix = np.array(final_errors)
    # ==================
    
    # ================================
    good_apertures = [list_apertures[i_ap] for i_ap in apertures_mask]
    number_of_good_apertures = len(good_apertures)
    print()
    print('Good apertures: ', number_of_good_apertures)
    print()
    # ================================
    
    # ==================
    # BUILD THE TABLE(S)
    ra_apertures, dec_apertures = ascii.read(coordinates_path)['ra'], ascii.read(coordinates_path)['dec']
    sub_ra_apertures, sub_dec_apertures = ra_apertures[apertures_mask], dec_apertures[apertures_mask]
    df = pd.DataFrame()
    df['ra'], df['dec'] = sub_ra_apertures, sub_dec_apertures
    table = pd.DataFrame(final_fluxes, columns = list_band)
    table['ra'], table['dec'] = sub_ra_apertures, sub_dec_apertures
    table = table[['dec'] + table.columns[:-1].tolist()]
    table = table[['ra'] + table.columns[:-1].tolist()]
    table.to_csv(photometry_folder_path+'/FINAL.csv')#, sep='\t')
    
    # Genero la lista degli id. delle aperture per Magphys e la "lista" dei redshift (uno)
    list_id, list_z = [], []
    for a, b, idx in zip(sub_ra_apertures, sub_dec_apertures, range(len(sub_ra_apertures))):
        list_id.append(str(idx)), list_z.append(source_z.value)
        
    if np.log10(number_of_good_apertures) > 4: list_id = ['{:06d}'.format(int(idx)) for idx in list_id]
    else: list_id = ['{:05d}'.format(int(idx)) for idx in list_id]
    
    df = pd.DataFrame()
    table = pd.DataFrame(final_fluxes, columns = list_band)
    table['z'] = list_z
    table = table[['z'] + table.columns[:-1].tolist()]
    for band, i_band in zip(list_band, np.arange(len(list_band)*2, step = 2)):
        i_column = i_band + 2
        header_column = band+'_ERROR'
#       error_column = rms_dictionary[band]       
        error_column = errors_matrix[:,int(i_band/2)]
        table.insert(i_column, header_column, error_column)      
    table['id'] = list_id
    table = table[['id'] + table.columns[:-1].tolist()]
    table.to_csv(photometry_folder_path+'/Final_for_MAGPHYS.csv', index=False, sep='\t')
    # ==================
    
    if remove_weird_points: remove_weirdities(photometry_folder_path+'/Final_for_MAGPHYS.csv')
    
    return
            
def build_table_NSR(galaxy_properties, list_band, cell_side, run_type, skip_value = 5, remove_weird_points = False):
    '''
    Building the MAGPHYS fluxes table, with errors coming from the measured noise-to-signal ratio of Dustpedia photometry.
    INPUT: galaxy properties, list of band, cell side (either in arcsec or kpc), run type. Optional: skip value (the one over which ignore the aperture)
            and the option to remove weird points in the SED.
    OUTPUT: it builds two files, 'FINAL.csv' and 'Final_for_MAGPHYS.csv' in the photometry folder.
    WRITTEN: A. Enia, Dec 2019
    '''
    # =================
    # GALAXY PROPERTIES
    galaxy_name = galaxy_properties.galaxy_name
    D_L = galaxy_properties.dist
    source_z = galaxy_properties.dist*cosmo.H0/v_lux.to('km/s')
    source_z = np.round(source_z, 6)
    ra_cen, dec_cen = galaxy_properties.ra, galaxy_properties.dec
    # =================
    
    # =================
    # Read the apertures
    if cell_side.unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_side.unit == 'arcmin': cell_side = cell_side.to('arcsec')

    photometry_folder_path = '../'+galaxy_name+'/Photometry_'+run_type+'/'
    coordinates_path = photometry_folder_path+'coordinates.txt'
    if cell_side.unit == 'kpc': cell_side *= cosmo.arcsec_per_kpc_comoving(source_z)
    list_apertures = read_apertures(coordinates_path, cell_side)
    # =================

    # =================
    # Read DustPedia photometry, to extract the SNR
    DustPedia_data = pd.read_csv('../'+galaxy_name+'/Reduction/'+galaxy_name+'_photometry.dat', sep = '\t')
    fluxes_Dustpedia = DustPedia_data['Flux'].values
    errors_Dustpedia = DustPedia_data['Error'].values
    NSR_Dustpedia = errors_Dustpedia/fluxes_Dustpedia
    flx_per_err = []
    for i_ap in tqdm(range(len(list_apertures))):
        data = ascii.read(photometry_folder_path+'/fluxes/ap_{0:004}_fluxes.txt'.format(i_ap))
        f = data['Fluxes'].data
        check = np.count_nonzero(f == -99.0)
        if check <= skip_value:
            ok = np.where(f == -99.0)
            f[ok] = 'nan'
            flx_per_err.append(f)
        else: pass
    #mean_flux_band = np.nanmean(flx_per_err, axis = 0)
    #error_per_band = mean_flux_band*NSR_Dustpedia
    error_per_band = flx_per_err*NSR_Dustpedia

    # Reread the results
    final_fluxes, final_errors, apertures_mask = [], [], []
    for i_ap in tqdm(range(len(list_apertures))):
        data = ascii.read(photometry_folder_path+'/fluxes/ap_{0:004}_fluxes.txt'.format(i_ap))
        try: f, e = data['Fluxes'].data, data['Errors'].data
        except: f, e = data['Fluxes'].data, 0
        ok = np.where(f == 0)
        f[ok] = -99.0
        check = np.count_nonzero(f == -99.0)
        if check <= skip_value:
            final_fluxes.append(f), final_errors.append(e), apertures_mask.append(i_ap)
        else: pass
    
    # ================================
    # Attach errors from DustPedia SNR
    # Read the results
    final_errors = np.array(final_errors)
    GOOD_err = np.where(final_errors != 0.0)
    error_per_band = np.array(final_fluxes*NSR_Dustpedia) # Each flux error is obtained from the band NSR of Dustpedia
    #for a,b in zip(GOOD_err[0], GOOD_err[1]):
    #    error_per_band[a,b] = final_errors[a,b]
    # ================================
    
    # ================================
    good_apertures = [list_apertures[i_ap] for i_ap in apertures_mask]
    number_of_good_apertures = len(good_apertures)
    print()
    print('Good apertures: ', number_of_good_apertures)
    print()
    # ================================
    
    # ==================
    # BUILD THE TABLE(S)
    ra_apertures, dec_apertures = ascii.read(coordinates_path)['ra'], ascii.read(coordinates_path)['dec']
    sub_ra_apertures, sub_dec_apertures = ra_apertures[apertures_mask], dec_apertures[apertures_mask]
    df = pd.DataFrame()
    df['ra'], df['dec'] = sub_ra_apertures, sub_dec_apertures
    table = pd.DataFrame(final_fluxes, columns = list_band)
    table['ra'], table['dec'] = sub_ra_apertures, sub_dec_apertures
    table = table[['dec'] + table.columns[:-1].tolist()]
    table = table[['ra'] + table.columns[:-1].tolist()]
    table.to_csv(photometry_folder_path+'/FINAL.csv')#, sep='\t')
    
    # Genero la lista degli id. delle aperture per Magphys e la "lista" dei redshift (uno)
    list_id, list_z = [], []
    for a, b, idx in zip(sub_ra_apertures, sub_dec_apertures, range(len(sub_ra_apertures))):
        list_id.append(str(idx)), list_z.append(source_z)
    list_id = ['{:05d}'.format(int(idx)) for idx in list_id]
    
    df = pd.DataFrame()
    table = pd.DataFrame(final_fluxes, columns = list_band)
    table['z'] = list_z
    table = table[['z'] + table.columns[:-1].tolist()]
    for band, error, i_band in zip(list_band, error_per_band.transpose(), np.arange(len(list_band)*2, step = 2)):
        i_column = i_band + 2
        header_column = band+'_ERROR'
        error_column = np.copy(error)
        table.insert(i_column, header_column, error_column)      
    table['id'] = list_id
    table = table[['id'] + table.columns[:-1].tolist()]
    table.to_csv(photometry_folder_path+'/Final_for_MAGPHYS.csv', index=False, sep='\t')
    # ==================

    if remove_weird_points: remove_weirdities(photometry_folder_path+'/Final_for_MAGPHYS.csv')

    return

