#!/usr/bin/python
# -*- coding: latin-1 -*-

import os, subprocess
import numpy as np
import GenericUsefulScripts as GUS
from astropy import units as u
from astropy.io import ascii, fits
from astropy.convolution import convolve
from astropy.stats import SigmaClip
from astropy.coordinates import SkyCoord
from photutils.background import MedianBackground, Background2D
from skimage.transform import resize
import multiprocessing
import ChrisFuncs
import pandas as pd
space = '       '

def data_reduction(galaxy_name, path_fits_input = 'standard'):

    # ---------------------------------------------------------------------------
    # Galaxy Aperture Stuff, from Dustpedia (to mask and bkg evaluation purposes)
    DustPedia_Photom = pd.read_csv('../DustPedia_Tables/DustPedia_Aperture_Photometry_2.2.csv')
    subtable = DustPedia_Photom.loc[DustPedia_Photom['name'] == galaxy_name]
    ra, dec = subtable['ra'].values[0], subtable['dec'].values[0]
    ap_cen_coord = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    semimaj = subtable['semimaj_arcsec'].values[0]
    axial_ratio, pos_angle = subtable['axial_ratio'].values[0], subtable['pos_angle'].values[0]
    # ---------------------------------------------------------------------------
    
    subprocess.call('mkdir ../'+galaxy_name+'/_ReducedMaps/', shell = True)
        
    list_data = []
    if path_fits_input == 'standard': path_fits_input = '../'+galaxy_name+'/Caapr/Temp/Processed_Maps'
    else: path_fits_input = '../'+galaxy_name+'/'+path_fits_input
    header_fits = '../'+galaxy_name+'/Caapr/Maps/'

    print('Reading original maps...')
    filelist = [x for x in os.listdir('Caapr/Maps') if x.endswith('.fits')]
    for file in filelist:
        if file.endswith('Thumbnail.fits'): continue  # Don't work with thumbnails
        elif file.endswith('Error.fits'): continue      # Don't work with Errors
        signal_path = path_fits_input+'/'+file
        list_data.append(GUS.FitsUtils(signal_path))
        print(space+signal_path+' read')
    print('...done!')
    print()
    
    for data in list_data:
        if os.path.exists('../'+galaxy_name+'/_ReducedMaps/'+data.bandname+'.fits'):
            print(data.bandname+'.fits already reduced, skipping to next band')
            continue
        else: print('Processing band', data.bandname)
            
        # Galaxy Aperture Stuff, from Dustpedia (to mask and bkg evaluation purposes)
        centre_x, centre_y = ap_cen_coord.to_pixel(data.wcs)
        pixel_scale = (data.get_pixel_scale()*u.deg).to('arcsec').value
        Gal_Ap_Stuff = centre_x, centre_y, semimaj/pixel_scale, axial_ratio, pos_angle

        # Reduce band
        signal_reduced = reduce(data, Gal_Ap_Stuff)
        
        # Save fits
        hdu = fits.PrimaryHDU(signal_reduced)
        hdu.header = data.hdr
        hdu.writeto('../'+galaxy_name+'/_ReducedMaps/'+data.bandname+'.fits')
    
    print()
    print('Data reduction phase over.')
    print()
    return


def data_reduction_parallel(galaxy_name, processes = 5, path_fits_input = 'standard'):
    from itertools import repeat

    # ---------------------------------------------------------------------------
    # Galaxy Aperture Stuff, from Dustpedia (to mask and bkg evaluation purposes)
    DustPedia_Photom = pd.read_csv('../DustPedia_Tables/DustPedia_Aperture_Photometry_2.2.csv')
    subtable = DustPedia_Photom.loc[DustPedia_Photom['name'] == galaxy_name]
    ra, dec = subtable['ra'].values[0], subtable['dec'].values[0]
    ap_cen_coord = SkyCoord(ra*u.deg, dec*u.deg, frame = 'fk5')
    semimaj = subtable['semimaj_arcsec'].values[0]
    axial_ratio, pos_angle = subtable['axial_ratio'].values[0], subtable['pos_angle'].values[0]
    # ---------------------------------------------------------------------------
    
    subprocess.call('mkdir ../'+galaxy_name+'/_ReducedMaps/', shell = True)
        
    list_data = []
    if path_fits_input == 'standard': path_fits_input = '../'+galaxy_name+'/Caapr/Temp/Processed_Maps'
    else: path_fits_input = '../'+galaxy_name+'/'+path_fits_input
    header_fits = '../'+galaxy_name+'/Caapr/Maps/'

    print('Reading original maps...')
    filelist = [x for x in os.listdir('Caapr/Maps') if x.endswith('.fits')]
    for file in filelist:
        if file.endswith('Thumbnail.fits'): continue  # Don't work with thumbnails
        elif file.endswith('Error.fits'): continue      # Don't work with Errors
        signal_path = path_fits_input+'/'+file
        list_data.append(GUS.FitsUtils(signal_path))
        print(space+signal_path+' read')
    print('...done!')
    print()
    
    pool = multiprocessing.Pool()
    with multiprocessing.Pool(processes=processes) as pool:
        func = zip(list_data, repeat(galaxy_name), \
                   repeat(ap_cen_coord), repeat(semimaj), repeat(axial_ratio), repeat(pos_angle))
        pool.starmap(reduction_loop_parallel, func)
        
    print()
    print('Data reduction phase over.')
    print()
    return

def reduction_loop_parallel(data, galaxy_name, ap_cen_coord, semimaj, axial_ratio, pos_angle):
    if os.path.exists('../'+galaxy_name+'/_ReducedMaps/'+data.bandname+'.fits'):
        print(data.bandname+'.fits already reduced, skipping to next band')
        return
    else: print('Processing band', data.bandname)
        
    # Galaxy Aperture Stuff, from Dustpedia (to mask and bkg evaluation purposes)
    centre_x, centre_y = ap_cen_coord.to_pixel(data.wcs)
    pixel_scale = (data.get_pixel_scale()*u.deg).to('arcsec').value
    Gal_Ap_Stuff = centre_x, centre_y, semimaj/pixel_scale, axial_ratio, pos_angle

    # Reduce band
    signal_reduced = reduce(data, Gal_Ap_Stuff)
    
    # Save fits
    hdu = fits.PrimaryHDU(signal_reduced)
    hdu.header = data.hdr
    hdu.writeto('../'+galaxy_name+'/_ReducedMaps/'+data.bandname+'.fits')
    return

def reduce(data, Gal_Ap_Stuff, psf_degrad = True, sky_sub = True):
    
            #if data.bandname[:7] == 'Spitzer':
            #    print 
            #    print(space+"Spitzer bands usually have a problem with sky subtraction")
            #    print(space+"Evaluated background average is "+str(bkg_average)+". Perhaps it's too low.")
            #    print(space+"Do you want to insert the bkg average by hand? (insert value or n)")
            #    answer = raw_input()
            #    if answer == 'n': pass
            #    else: bkg_average = float(answer)
            #else: pass
        
    ok_nan = np.where(np.nan_to_num(data.signal_with_nans-1) == 0) # I know, can't do anything 'bout it
    if sky_sub:
        # Sky subtraction
        print(space+'Sky subtraction for '+data.bandname+' ...')
        # 1) Flatten the background
        signal_flat, check_sub = sky_flattening(data, Gal_Ap_Stuff)
        # 2) If check_sub is sub, the sky has already been flattened + removed
        #    if not, remove the average background
        if check_sub == 'sub':
            signal_skysub = signal_flat.copy()
        elif check_sub == 'unsub':
            bkg_average = evaluate_bkg_avg(signal_flat, Gal_Ap_Stuff)
            if bkg_average < 0:
                print(space+"Evaluated background average is lower than 0. Returning original map.")
                signal_skysub = signal_flat.copy()
            else:
                print(space+"Evaluated background average is {0:.2E}".format(bkg_average))
                signal_skysub = signal_flat - bkg_average
    else: 
        print(space+'No sky flattening + subtraction requested. Hey, whatever you want.')
        signal_skysub = data.signal.copy()

    if psf_degrad:
        print(space+'PSF degradation for '+data.bandname+' ...')
        if data.bandname == 'SPIRE_350':
            return signal_skysub
        else:
            try:
                kernel_path = '../_kernels/Kernel_LoRes_'+data.bandname+'_to_SPIRE_350.fits'
                kernel = fits.getdata(kernel_path)
                kernel_resized = resize(kernel, (101, 101), preserve_range = True)
                signal_conv = convolve(signal_skysub, kernel = kernel_resized, boundary = None, preserve_nan = True)
                signal_conv[ok_nan] = np.nan
            except: 
                print(space+'No LowResolution kernel, switching to (slower) HighResolution.')
                kernel_path = '../_kernels/Kernel_HiRes_'+data.bandname+'_to_SPIRE_350.fits'
                kernel = fits.getdata(kernel_path)
                kernel_resized = resize(kernel, (101, 101), preserve_range = True)
                signal_conv = convolve(signal_skysub, kernel = kernel_resized, boundary = None, preserve_nan = True)   
                signal_conv[ok_nan] = np.nan
            return signal_conv
    else: 
        print(space+'No PSF degradation requested. I beg you to reconsider.')
        signal_skysub[ok_nan] = np.nan
        return signal_skysub

def sky_flattening(data, Gal_Ap_Stuff):
    from astropy.modeling.polynomial import Polynomial2D
    from astropy.modeling.fitting import LevMarLSQFitter
    from scipy.ndimage.interpolation import zoom
    
    # 1) Read data, get pixel scale
    image = data.signal_with_nans
    pix_size = (data.get_pixel_scale()*u.deg).to('arcsec').value
    bandname = data.bandname
    
    # 2) If image has pixels smaller than some limit, downsample image to improve processing time
    pix_size_limit = 2.0
    if pix_size<pix_size_limit: downsample_factor = int(np.ceil(pix_size_limit/pix_size))
    else: downsample_factor = 1
    image_ds = GUS.Downsample(image, downsample_factor)
    
    # 3) Sigma clip the downsampled image
    clip_value = GUS.SigmaClip(image_ds, tolerance=0.01, sigma_thresh=3.0, median=True)
    noise_value = clip_value[0]
    field_value = clip_value[1]
    cutoff_sigma = 2.0
    cutoff = field_value + ( cutoff_sigma * noise_value )
    
    # 4) Mask the image removing galaxy emission...
    image_masked = image_ds.copy()
    centre_i, centre_j, mask_semimaj_pix, mask_axial_ratio, mask_angle = Gal_Ap_Stuff
    ellipse_mask = EllipseMask(image_ds, mask_semimaj_pix/downsample_factor, mask_axial_ratio, mask_angle, centre_i/downsample_factor, centre_j/downsample_factor)
    image_masked[ np.where( ellipse_mask==1 ) ] = np.nan
    
    # ...and image pixels identified as having high SNR
    image_masked[ np.where( image_masked>cutoff ) ] = np.nan
    
    # 5) Use astropy to set up 2-dimensional polynomial to the image
    image_masked[ np.where( np.isnan(image_masked)==True ) ] = field_value
    poly_model = Polynomial2D(degree=5)
    i_coords, j_coords = np.mgrid[:image_masked.shape[0], :image_masked.shape[1]]
    fitter = LevMarLSQFitter()
    i_coords = i_coords.flatten()
    j_coords = j_coords.flatten()
    image_flattened = image_masked.flatten()
    good = np.where(np.isnan(image_flattened)==False)
    i_coords = i_coords[good]
    j_coords = j_coords[good]
    
    # 6) Attempt polynomial fit; if insufficient data then skip onwards
    image_flattened = image_flattened[good]
    try:
        fit = fitter(poly_model, i_coords, j_coords, image_flattened)
    except:
        print(space+'Error fitting polinomial sky model. Returning unalterated image.')
        return image
        
    # 7) Create final polynomial filter (undoing downsampling using lorenzoriano GitHub script)
    i_coords, j_coords = np.mgrid[:image_ds.shape[0], :image_ds.shape[1]]
    poly_fit = fit(i_coords, j_coords)
    poly_full = zoom(poly_fit, [ float(image.shape[0])/float(poly_fit.shape[0]), \
                                float(image.shape[1])/float(poly_fit.shape[1]) ], mode='nearest') 
    
    # 8) Establish background variation before application of filter
    sigma_thresh = 3.0
    clip_in = GUS.SigmaClip(image, tolerance=0.005, median=True, sigma_thresh=sigma_thresh)
    bg_in = image[ np.where( image<clip_in[1] ) ]
    spread_in = np.mean( np.abs( bg_in - clip_in[1] ) )
    
    # 9) How much reduction in background variation there was due to application of the filter
    image_sub = image - poly_full
    clip_sub = GUS.SigmaClip(image_sub, tolerance=0.005, median=True, sigma_thresh=sigma_thresh)
    bg_sub = image_sub[ np.where( image_sub < clip_sub[1] ) ]
    spread_sub = np.mean( np.abs( bg_sub - clip_sub[1] ) )
    spread_diff = spread_in / spread_sub
    
    # If the filter made significant difference, apply to image and return it; otherwise, just return the unaltered map
    if spread_diff>1.1:
        print(space+bandname+' background is significantly variable; removing polynomial background fit.')
        return image_sub, 'sub'
    else:
        print(space+bandname+' background is not significantly variable; leaving image unaltered.')
        return image, 'unsub'

def evaluate_bkg_avg(image, Gal_Ap_Stuff):
    '''
    Function to evaluate the mean background in an elliptical annulus between 1.25 and 1.601 times the galaxy semimajor axis (from DustPedia photometric table).
    Args: Array, semi-major axis of inside edge of annulus (pix), width of annulus (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse
    Returns: Numpy array containing the mean background per pixel.
    '''
    centre_x, centre_y, semimaj_pix, axial_ratio, pos_angle = Gal_Ap_Stuff
    # =========
    # Evaluate pixels in background annulus
    bg_inner_semimaj_pix = semimaj_pix * 1.25
    bg_width = (semimaj_pix * 1.601) - bg_inner_semimaj_pix
    bg_calc = AnnulusSum(image, bg_inner_semimaj_pix, bg_width, axial_ratio, pos_angle, centre_x, centre_y)
    bg_clip = GUS.SigmaClip(bg_calc[2], median=False, sigma_thresh=3.0)
    # =========
    return bg_clip[1]

def check_Dustpedia(galaxy_name, working_bands):
    '''
    Function to check if DustPedia photometric flux and the one measured in the same apertures with our data reduction are compatible. 
    Args: Galaxy name, working bands, if wanted, perform Galactic Extinction Correction
    Returns: Nothing, generates a plot in Reduction folder.
    '''
    import os, subprocess
    from astropy.io import fits, ascii
    from astropy import units as u
    import pandas as pd
    import numpy as np
    from photutils import SkyEllipticalAperture, SkyEllipticalAnnulus, aperture_photometry
    from astropy.coordinates import SkyCoord
    from matplotlib import pyplot as plt

    subprocess.call('mkdir ../'+galaxy_name+'/Reduction/', shell = True)
    path_galaxy_photometry = '../'+galaxy_name+'/Reduction/'+galaxy_name+'_photometry.dat'

    # =========
    # Read DustPedia Photometric Table
    DustPedia_Photom = pd.read_csv('../DustPedia_Tables/DustPedia_Aperture_Photometry_2.2.csv')
    # Rearrange in order of increasing effective wavelenght
    right_order = [u'name', u'ra', u'dec', u'semimaj_arcsec', u'axial_ratio', u'pos_angle', u'global_flag',
    u'GALEX_FUV', u'GALEX_FUV_err', u'GALEX_FUV_flag', u'GALEX_NUV', u'GALEX_NUV_err', u'GALEX_NUV_flag',
    u'SDSS_u', u'SDSS_u_err', u'SDSS_u_flag', u'SDSS_g', u'SDSS_g_err', u'SDSS_g_flag',
    u'SDSS_r', u'SDSS_r_err', u'SDSS_r_flag', u'SDSS_i', u'SDSS_i_err', u'SDSS_i_flag',
    u'SDSS_z', u'SDSS_z_err', u'SDSS_z_flag',
    u'2MASS_J', u'2MASS_J_err', u'2MASS_J_flag', u'2MASS_H', u'2MASS_H_err', u'2MASS_H_flag',
    u'2MASS_Ks', u'2MASS_Ks_err', u'2MASS_Ks_flag',
    u'WISE_3.4', u'WISE_3.4_err', u'WISE_3.4_flag', u'Spitzer_3.6', u'Spitzer_3.6_err', u'Spitzer_3.6_flag',
    u'Spitzer_4.5', u'Spitzer_4.5_err', u'Spitzer_4.5_flag', u'WISE_4.6', u'WISE_4.6_err', u'WISE_4.6_flag',
    u'Spitzer_5.8', u'Spitzer_5.8_err', u'Spitzer_5.8_flag', u'Spitzer_8.0', u'Spitzer_8.0_err', u'Spitzer_8.0_flag',
    u'WISE_12', u'WISE_12_err', u'WISE_12_flag', u'WISE_22', u'WISE_22_err', u'WISE_22_flag',
    u'Spitzer_24', u'Spitzer_24_err', u'Spitzer_24_flag', u'Spitzer_70', u'Spitzer_70_err', u'Spitzer_70_flag',
    u'PACS_70', u'PACS_70_err', u'PACS_70_flag', u'PACS_100', u'PACS_100_err', u'PACS_100_flag',
    u'PACS_160', u'PACS_160_err', u'PACS_160_flag', u'Spitzer_160', u'Spitzer_160_err', u'Spitzer_160_flag',
    u'SPIRE_250', u'SPIRE_250_err', u'SPIRE_250_flag', u'SPIRE_350', u'SPIRE_350_err', u'SPIRE_350_flag',
    u'SPIRE_500', u'SPIRE_500_err', u'SPIRE_500_flag']
    DustPedia_Photom = DustPedia_Photom[right_order]
    gal_phot = DustPedia_Photom.loc[DustPedia_Photom['name'] == galaxy_name]
    # Fist, remove _flag columns
    to_remove = gal_phot.columns.str.contains('flag', case=False)
    gal_phot = gal_phot.loc[:,~to_remove]
    # Extract ra, dec, semimaj, axial ratio and pos_angle, then remove them
    ra, dec = gal_phot['ra'].values[0], gal_phot['dec'].values[0]
    semimaj, axial_ratio, pos_angle = gal_phot['semimaj_arcsec'].values[0], gal_phot['axial_ratio'].values[0], gal_phot['pos_angle'].values[0]
    to_remove = ['name', 'ra', 'dec', 'semimaj_arcsec', 'axial_ratio', 'pos_angle']
    gal_phot = gal_phot.drop(columns=to_remove)
    # And remove empy columns
    #gal_phot = gal_phot.dropna(axis='columns')
    # Extract working bands fluxes and errors
    gal_phot_flux = gal_phot[working_bands]
    gal_phot_flux = gal_phot_flux.transpose()
    working_bands_err = [t+'_err' for t in working_bands]
    gal_phot_err = gal_phot[working_bands_err]
    gal_phot_err = gal_phot_err.transpose()
    galaxy_photometry = pd.DataFrame(np.concatenate((gal_phot_flux.values, gal_phot_err.values), axis=1))
    galaxy_photometry.columns = ['Flux', 'Error']
    galaxy_photometry.index = working_bands
    galaxy_photometry = galaxy_photometry.fillna(0) # Fill NaN entries with zeroes
    # Save
    galaxy_photometry.index.names = ['Band'] # Rename the index column as "Band"
    galaxy_photometry.to_csv(path_galaxy_photometry, sep='\t', index = False)
    # =========

    # =========
    # APERTURES
    # Read the apertures + radii
    positions = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    DustPedia_aperture = SkyEllipticalAperture(positions, a=semimaj*u.arcsec, b=semimaj*u.arcsec/axial_ratio, theta=pos_angle*u.deg)
    DustPedia_annulus = SkyEllipticalAnnulus(positions, a_in=semimaj*u.arcsec*1.25, a_out=semimaj*u.arcsec*1.601, \
                                       b_out=semimaj*u.arcsec/axial_ratio, theta=pos_angle*u.deg)
    # =========

    # =========
    # Galactic Extinction Correction dictionary
    GalCorr_path = '../'+galaxy_name+'/galactic_extinction_correction.txt'
    if os.path.exists(GalCorr_path): pass
    else: GalExtCorr(galaxy_name, working_bands, ra, dec)
    GalCorrection_dictionary = dict(zip(ascii.read(GalCorr_path)['Band'].data, \
                      ascii.read(GalCorr_path)['Correction'].data))
    # =========

    # =========
    # Read reduced data and perform photometry
    path_fits = '../'+galaxy_name+'/_ReducedMaps/'
    list_data = []
    for file in os.listdir(path_fits):
        if not file.endswith('.fits'): continue
        elif file.startswith('In'): continue
        list_data.append(GUS.FitsUtils(path_fits+file))
    
    list_fluxes = []
    for data in list_data:
        # Perform photometry
        phot_table = aperture_photometry(data.signal, DustPedia_aperture, wcs = data.wcs)
        phot_table['aperture_sum'].info.format = '%.4g'  
        # Put results in a single file
        phot = GUS.round_arr(phot_table['aperture_sum'].data, 2)
        # Galactic extintion correction
        phot *= GalCorrection_dictionary[data.bandname]
        list_fluxes.append(abs(phot))
        
    fluxes = np.array(list_fluxes)
        
    # Sort w.r.t wavelengths
    list_wvl = (t.get_wavelength() for t in list_data)
    list_band = (t.bandname for t in list_data)
    wvl, fluxes, bandnames = (t for t in zip(*sorted(zip(list_wvl, fluxes, list_band))))
    wvl, fluxes = np.array(wvl), np.array(fluxes)[:,0]

    # Save the results
    ascii.write([bandnames, GUS.round_arr(wvl,2), GUS.round_arr(fluxes, 2)], '../'+galaxy_name+'/Reduction/'+galaxy_name+'_fluxes.txt', \
                names = ['Band', 'Wvl', 'Fluxes'], overwrite=True)
    # =========

    # =========
    # Re-read Dustpedia Photometry
    data_CAAPR = ascii.read(path_galaxy_photometry)
    fluxes_CAAPR, errors_CAAPR = data_CAAPR['Flux'].data, data_CAAPR['Error'].data
    compatibility = np.abs(np.array(fluxes_CAAPR) - np.array(fluxes))/np.sqrt(np.array(errors_CAAPR)**2)
    ascii.write([GUS.round_arr(compatibility,2)], '../'+galaxy_name+'/Reduction/'+galaxy_name+'_comp.txt', format='fixed_width_two_line', \
                names = ['Comp'], overwrite=True)
    # =========

    # =========
    # Plot    
    xmin, xmax = np.array(wvl).min(), np.array(wvl).max()

    DustpediaCheckPlot = plt.figure(figsize=(15,5))
    plt.subplot(2,1,1)
    plt.plot(np.array(wvl), np.array(fluxes_CAAPR), \
             linestyle = 'None', marker = '.', color = 'navy', label = 'CAAPR+Literature Photometry')
    plt.plot(wvl, fluxes, linestyle = 'None', marker = '.', color = 'red', label = 'My Photometry')
    plt.xscale('log'), plt.yscale('log')
    plt.ylabel(r'Flux (Jy)')
    plt.legend()
    
    plt.subplot(2,1,2)
    plt.axhline(5, color = 'r', linestyle = '-')
    plt.plot(wvl, compatibility,  ms = 10.0, linestyle = 'None', color = 'k', marker = '.')
    for i in range(len(wvl)):
        plt.text(wvl[i], 0.5, bandnames[i], rotation = 90)
    plt.xscale('log'), plt.yscale('log')
    plt.xlabel(r'Wavelength ($\mu$m)'), plt.ylabel(r'Compatibility $\lambda$')
    plt.subplots_adjust(hspace=0.,wspace=0.)
    DustpediaCheckPlot.savefig('../'+galaxy_name+'/Reduction/'+galaxy_name+'_SED.pdf', bbox_inches = 'tight')
    # =========
    
    return

def GalExtCorr(galaxy_name, list_band, ra, dec):
    list_correction = []
    for band in list_band:
        try:
            if band == 'Spitzer_3.6': band = 'IRAC1'
            elif band == 'Spitzer_4.5': band = 'IRAC2'
            elif band == 'Spitzer_5.8': band = 'IRAC3'
            elif band == 'Spitzer_8.0': band = 'IRAC4'
            elif band == 'WISE_3.4': band = 'WISE1'
            elif band == 'WISE_4.6': band = 'WISE2'
            correction = ChrisFuncs.ExtCorrrct(ra, dec, band, verbose = False)
            list_correction.append(correction)
        except: list_correction.append(1)
    ascii.write([list_band, list_correction], \
                '../'+galaxy_name+'/galactic_extinction_correction.txt', names = ['Band', 'Correction'])
    return 

##################################
# QUI COPIO BRUTALMENTE DA CLARK #
##################################
def AnnulusSum(array, rad_inner, width, axial_ratio, angle, i_centre, j_centre):
    '''
    Function to sum all elements in an annulus centred upon the middle of the given array
    Args: Array, semi-major axis of inside edge of annulus (pix), width of annulus (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse
    Returns: Numpy array containing the sum of the pixel values in the annulus, the total number of pixels counted, and an array containing the pixel values
    '''
    # Create slice of input array, containing only the region of interest
    i_cutout_min = int(np.floor(max([0, i_centre-(rad_inner+width)])))
    i_cutout_max = int(np.ceil(min([(array.shape)[0], i_centre+(rad_inner+width)])))
    j_cutout_min = int(np.floor(max([0, j_centre-(rad_inner+width)])))
    j_cutout_max = int(np.ceil(min([(array.shape)[1], j_centre+(rad_inner+width)])))
    array_slice = array[ int(round(i_cutout_min)):int(round(i_cutout_max))+1, int(round(j_cutout_min)):int(round(j_cutout_max))+1 ]
    i_centre_slice = i_centre - i_cutout_min
    j_centre_slice = j_centre - j_cutout_min
    if array[int(i_centre),int(j_centre)]!=array_slice[int(i_centre_slice),int(j_centre_slice)]:
        if np.isnan(array[int(i_centre),int(j_centre)]==False) and np.isnan(array_slice[int(i_centre_slice),int(j_centre_slice)]==False):
            print('SEVERE ERROR: AnnulusSum check failed.')
            pdb.set_trace()
    else:
        array = array_slice
        i_centre = i_centre_slice
        j_centre = j_centre_slice

    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj_inner = float(rad_inner)
    semi_min_inner = float(semi_maj_inner) / float(axial_ratio)
    semi_maj_outer = float(rad_inner) + float(width)
    semi_min_outer  = float(semi_maj_outer) / float(axial_ratio)
    angle = np.radians(float(angle))

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, array.shape[0]-1, array.shape[0])
    j_linespace = np.linspace(0, array.shape[1]-1, array.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create array identifying which coordinates lie within inner ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check_inner = (j_trans**2 / semi_maj_inner**2) + (i_trans**2 / semi_min_inner**2 )

    # Use meshgrids to create array identifying which coordinates lie within outer ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check_outer = (j_trans**2 / semi_maj_outer**2) + (i_trans**2 / semi_min_outer**2 )

    # Calculate flux & pixels in aperture, and store pixel values
    annulus_where = np.where( (ellipse_check_outer<=1) & (ellipse_check_inner>1) & (np.isnan(array)==False) )
    annulus_tot = sum( array[ annulus_where ] )
    annulus_count = annulus_where[0].shape[0]
    annulus_pix = array[ annulus_where ]
    annulus_nan = np.where( (ellipse_check_outer<=1) & (ellipse_check_inner>1) & (np.isnan(array)==True) )

    # Return results
    return [annulus_tot, annulus_count, annulus_pix, annulus_nan]

def EllipseMask(array, rad, axial_ratio, angle, i_centre, j_centre):
    '''
    Function to return a mask identifying all pixels within an ellipse of given parameters
    Args: Array, semi-major axis (pix), axial ratio, position angle (deg), i & j coords of centre of ellipse
    Returns: Mask array of same dimensions as input array where pixels that lie within ellipse have value 1
    '''
    # Define semi-major & semi-minor axes, then convert input angle to radians
    semi_maj = float(rad)
    semi_min = float(rad) / float(axial_ratio)
    if angle.dtype != 'float': angle = float(angle.value)
    
    try:
        if angle.unit == 'rad': pass
        else: angle = np.radians(angle) # Convert the angle in radians
    except: angle = np.radians(angle) # Vabbè, assumo che sia da convertire e sticazzi

    # Create meshgrids with which to access i & j coordinates for ellipse calculations
    i_linespace = np.linspace(0, array.shape[0]-1, array.shape[0])
    j_linespace = np.linspace(0, array.shape[1]-1, array.shape[1])
    i_grid, j_grid = np.meshgrid(i_linespace, j_linespace, indexing='ij')

    # Use meshgrids to create array identifying which coordinates lie within ellipse
    i_trans = -(j_grid-float(j_centre))*np.sin(angle) + (i_grid-float(i_centre))*np.cos(angle)
    j_trans = (j_grid-float(j_centre))*np.cos(angle) + (i_grid-float(i_centre))*np.sin(angle)
    ellipse_check = (j_trans**2 / semi_maj**2) + (i_trans**2 / semi_min**2 )

    # Create ellipse mask
    ellipse_mask = np.zeros([array.shape[0], array.shape[1]])
    ellipse_mask[ np.where( ellipse_check<=1 ) ] = 1.0

    # Return array
    return ellipse_mask

def CircleSum(fits, i_centre, j_centre, r):
    '''
    Function to sum all pixel elements inside a given circle... the old-fashioned way
    Args: Array to be used, i & j coordinates of centre of circle, radius of circle
    Returns: Sum of elements within circle, number of pixels within circle
    '''    
    i_centre, j_centre, r = int(i_centre), int(j_centre), int(r)
    ap_sum = 0.0
    ap_pix = 0.0
    ap_values = []
    for i in range(-r, r+1):
        for j in range(-r, r+1):
            if i**2.0 + j**2.0 <= r**2.0:
                try:
                    ap_sum += fits[i_centre+i, j_centre+j]
                    ap_pix += 1.0
                    ap_values.append(fits[i_centre+i, j_centre+j])
                except:
                    continue
    return [ap_sum, ap_pix, ap_values]
