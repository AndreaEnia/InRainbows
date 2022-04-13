#!/usr/bin/python
# -*- coding: latin-1 -*-

import pandas as pd
DustPath = '/home/dustpedia/DustFolder/'
bandwavpath = DustPath+'scripts/bands_and_wvlghts.txt'
dict_df = pd.read_csv(bandwavpath, delimiter = '\t')
Wvlghts_dictionary = dict(zip(dict_df['name'], dict_df['lambda_eff']))

class GalaxyProperties:
    def __init__(self, galaxy_name, table_path = DustPath+'/DustPedia_Tables/DustPedia_HyperLEDA_Herschel.csv', cosmology = 'Planck15'):
        import pandas as pd
        from astropy import units as u
        from astropy.constants import c as v_lux
        if cosmology == 'Planck15':
            from astropy.cosmology import Planck15 as cosmo
        elif cosmology == 'FlatLambdaCDM':
            from astropy.cosmology import FlatLambdaCDM
            cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)  
            
        self.galaxy_name = galaxy_name   
        df = pd.read_csv(table_path)
        self.properties_table = df.loc[df['objname'] == self.galaxy_name]
        self.ra = self.properties_table['ra2000'].values[0]*u.deg
        self.dec = self.properties_table['de2000'].values[0]*u.deg
        self.t = self.properties_table['t'].values[0]
        self.d25 = self.properties_table['d25'].values[0]*u.arcmin
        self.incl = self.properties_table['incl'].values[0]
        self.dist = self.properties_table['dist_best'].values[0]*u.Mpc
        self.z_source = cosmo.H0*self.dist/v_lux.to('km/s')
        self.d25_physical = cosmo.kpc_comoving_per_arcmin(self.z_source)*self.d25

class FitsUtils:
    def __init__(self, signal_path):
        import numpy as np
        from astropy.io import fits, ascii
        from astropy.wcs import WCS
        self.fits_path = signal_path
        self.signal = np.nan_to_num(fits.getdata(self.fits_path))
        self.signal_with_nans = fits.getdata(self.fits_path)
        self.hdr = fits.getheader(self.fits_path)
        self.wcs = WCS(self.hdr)
        
        if signal_path.rsplit('/',1)[-1][0:3] == 'NGC': self.bandname = signal_path.rsplit('/',1)[1][8:-5]
        else: self.bandname = signal_path.rsplit('/',1)[1][:-5]
    
        error_path = signal_path.rsplit('/',1)[0]+'/ERROR_MAPS/'+signal_path.rsplit('/',1)[-1]
        error_path = error_path.rsplit('.',1)[0]+'_Error.'+error_path.rsplit('.',1)[1]
        try:
            self.errormap = fits.getdata(error_path)
        except:
            self.errormap = 0
    
    def get_pixel_scale(self):
        import numpy as np
        if ('CDELT1' in self.hdr) & ('CDELT2' in self.hdr):
            pixel_scale_x=abs(self.hdr['CDELT1'])
            pixel_scale_y=abs(self.hdr['CDELT2'])
        elif ('CD1_1' in self.hdr) & ('CD1_2' in self.hdr) & \
                ('CD2_1' in self.hdr) & ('CD2_2' in self.hdr):
            _ = np.arctan(self.hdr['CD2_1']/self.hdr['CD1_1'])
            pixel_scale_x = abs(self.hdr['CD1_1']/np.cos(_))
            pixel_scale_y = abs(self.hdr['CD2_2']/np.cos(_))
        else:
            raise ValueError
        pixel_scale = np.sqrt(pixel_scale_x*pixel_scale_y)
        return pixel_scale
    
    def get_wavelength(self):
        return Wvlghts_dictionary[self.bandname]

class AngularPhysicalConv:
    
    def __init__(self, z_source, cosmology_name = 'Planck15'):
        if cosmology_name == 'Planck15':
            from astropy.cosmology import Planck15 as cosmo
        elif cosmology_name == 'FlatLambdaCDM':
            from astropy.cosmology import FlatLambdaCDM
            cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Tcmb0=2.725 * u.K, Om0=0.3)
        else: raise ValueError('You should give a valid cosmology name.')
            
        self.cosmo = cosmo
        
        from astropy import units as u
        self.z_source = z_source

    def P2A(self, physical):
        try: _ = physical.unit
        except: raise ValueError('Physical quantity should have a unit.')
            
        if physical.unit == 'arcsec' or physical.unit == 'arcmin' or physical.unit == 'deg':
            raise ValueError('Give physical quantity, not angular, mona.')

        if physical.unit == 'pc': physical = physical.to('kpc')
        elif physical.unit == 'Mpc': physical = physical.to('kpc')
        angular = physical*self.cosmo.arcsec_per_kpc_comoving(self.z_source)
        return angular
    
    def A2P(self, angular):
        try: _ = angular.unit
        except: raise ValueError('Angular quantity should have a unit.')
            
        if angular.unit == 'pc' or angular.unit == 'kpc':
            raise ValueError('Give angular quantity, not physical, mona.')
        
        if angular.unit == 'arcsec': angular = angular.to('arcmin')
        elif angular.unit == 'deg': angular = angular.to('arcmin')
        physical = angular*self.cosmo.kpc_comoving_per_arcmin(self.z_source)
        return physical

class ManageTable:

    def __init__(self, galaxy_name, properties_path, table_path, chisq_threshold = 25):
        import sys
        sys.path.insert(0, '../') # To import Dustpedia Scripts
        sys.path.insert(0, '../scripts/') # To import Dustpedia Scripts
        import DustpediaScripts as DS
        import pandas as pd
        import numpy as np
        from astropy.cosmology import Planck15 as cosmo
        from astropy import units as u
        from astropy.constants import c as v_lux

        self.properties_path = properties_path
        self.galaxy_name = galaxy_name   
        df = pd.read_csv(properties_path)
        self.properties_table = df.loc[df['objname'] == self.galaxy_name]
        self.ra = self.properties_table['ra2000'].values[0]*u.deg
        self.dec = self.properties_table['de2000'].values[0]*u.deg
        self.t = self.properties_table['t'].values[0]
        self.r25 = .5*self.properties_table['d25'].values[0]*u.arcmin
        self.incl = self.properties_table['incl'].values[0]
        self.dist = self.properties_table['dist_best'].values[0]*u.Mpc
        self.z_source = cosmo.H0*self.dist/v_lux.to('km/s')
        self.r25_physical = cosmo.kpc_comoving_per_arcmin(self.z_source)*self.r25

        # Read Pixels Properties
        self.table_path = table_path
        self.table = pd.read_csv(self.table_path, sep = '\,', dtype={'Id. Aperture': object})
        self.chisq_threshold = chisq_threshold
        self.GOOD_chisq = np.where(self.table.Bestfit_Chisq < self.chisq_threshold)[0]
        self.BAD_chisq = np.where(self.table.Bestfit_Chisq >= self.chisq_threshold)[0]

    def update_MSdistance(self, fit_m, fit_q):
        import sys
        sys.path.insert(0, '../') # To import Dustpedia Scripts
        import DustpediaScripts as DS
        import numpy as np
        x, y = np.log10(self.table.SigmaMstar), np.log10(self.table.SigmaSFR)
        self.table['MS_distance'] = DS.shortest_distance(x, y, fit_m, fit_q)
        return

    def update_MSdistance(self, fit_m, fit_q):
        import sys
        sys.path.insert(0, '../') # To import Dustpedia Scripts
        import DustpediaScripts as DS
        import numpy as np
        x, y = np.log10(self.table.SigmaMstar), np.log10(self.table.SigmaSFR)
        self.table['MS_distance'] = DS.shortest_distance(x, y, fit_m, fit_q)
        return
    
    def update_MSdistance_ODR(self, fit_m, fit_q):
        import sys
        sys.path.insert(0, '../') # To import Dustpedia Scripts
        import DustpediaScripts as DS
        import numpy as np
        x, y = np.log10(self.table.SigmaMstar), np.log10(self.table.SigmaSFR)
        self.table['MS_distance_ODR'] = DS.shortest_distance(x, y, fit_m, fit_q)
        return
    
    def update_autoMSdistance(self, logSFR_limit = False):
        import sys
        sys.path.insert(0, '../') # To import Dustpedia Scripts
        import DustpediaScripts as DS
        import numpy as np
        from scipy import odr
        linear = odr.Model(DS.linear_function)
        x, y = np.log10(self.table.SigmaMstar[self.GOOD_chisq]), np.log10(self.table.SigmaSFR[self.GOOD_chisq])
        if logSFR_limit != False:
            ok = np.where(y >= logSFR_limit)
            x, y = x[ok], y[ok]

        MYodr = odr.ODR(odr.RealData(x, y), linear, beta0=[1,1])
        MYresults = MYodr.run()
        x, y = np.log10(self.table.SigmaMstar), np.log10(self.table.SigmaSFR)
        self.table['autoMS_distance'] = DS.shortest_distance(x, y, MYresults.beta[0], MYresults.beta[1])
        return MYresults.beta[0], MYresults.beta[1]
            
    def save_table(self, save_table_path, update = False):
        import pandas as pd
        import numpy as np
        if update == True: self.table.to_csv(self.table_path, index = False, sep = ',')
        else: self.table.to_csv(save_table_path, index = False, sep = ',')
        return
    
    def apply_M51b_mask(self):
        import numpy as np
        from astropy import units as u
        BAD_5194 = np.where(self.dec > 47.24052*u.deg)
        self.table.bestfit_chisq[BAD_5194] = 50
        self.GOOD_chisq = np.where(self.table.Bestfit_Chisq < self.chisq_threshold)
        self.BAD_chisq = np.where(self.table.Bestfit_Chisq >= self.chisq_threshold)
        return

    def apply_radius_threshold(self, radius_threshold):
        import numpy as np
        from astropy import units as u
        BAD = np.where(self.table.cendist_physical > radius_threshold.value)
        self.table.bestfit_chisq[BAD] = 100
        self.GOOD_chisq = np.where(self.table.Bestfit_Chisq < self.chisq_threshold)
        self.BAD_chisq = np.where(self.table.Bestfit_Chisq >= self.chisq_threshold)
        return
        
import numpy as np
import os, subprocess

def band_disambiguation(band):
    if band == 'GALEX_FUV' or band == 'GALEX FUV': return 'GALEXFUV'
    if band == 'GALEX_NUV' or band == 'GALEX NUV': return 'GALEXNUV'
    if band == 'SDSS u' or band == 'Sloan u': return 'SDSS_u'
    if band == 'SDSS g' or band == 'Sloan g': return 'SDSS_g'
    if band == 'SDSS r' or band == 'Sloan r': return 'SDSS_r'
    if band == 'SDSS i' or band == 'Sloan i': return 'SDSS_i'
    if band == 'SDSS z' or band == 'Sloan z': return 'SDSS_z'
    if band == '2MASS_J' or band == '2MASS J': return 'J'
    if band == '2MASS_H' or band == '2MASS H': return 'H'
    if band == '2MASS_Ks' or band == '2MASS Ks' or band == '2MASS K': return 'Ks'
    if band == 'WISE_3.4' or band == 'WISE 3.4': return 'WISEW1'
    if band == 'WISE_4.6' or band == 'WISE 4.6': return 'WISEW2'
    if band == 'WISE_12' or band == 'WISE 12': return 'WISEW3'
    if band == 'WISE_22' or band == 'WISE 22': return 'WISEW4'
    if band == 'Spitzer_3.6' or band == 'Spitzer 3.6': return 'IRAC1'
    if band == 'Spitzer_4.5' or band == 'Spitzer 4.5': return 'IRAC2'
    if band == 'Spitzer_5.8' or band == 'Spitzer 5.8': return 'IRAC3'
    if band == 'Spitzer_8.0' or band == 'Spitzer 8.0': return 'IRAC4'
    if band == 'Spitzer_24' or band == 'Spitzer 24': return 'MIPS24'
    if band == 'Spitzer_70' or band == 'Spitzer 70': return 'MIPS70'
    if band == 'Spitzer_160' or band == 'Spitzer 160': return 'MIPS160'
    if band == 'PACS_70' or band == 'PACS 70': return 'PACS70'
    if band == 'PACS_100' or band == 'PACS 100': return 'PACS100'
    if band == 'PACS_160' or band == 'PACS 160': return 'PACS160'
    if band == 'SPIRE_250' or band == 'SPIRE 250': return 'SPIRE250'
    if band == 'SPIRE_350' or band == 'SPIRE 350': return 'SPIRE350'
    if band == 'SPIRE_500' or band == 'SPIRE 500': return 'SPIRE500'
    if band == 'SMA' or band == 'SMA_800': return 'SMA800'
    return band

def associate_colormap(band):
    from matplotlib import cm
    band = band_disambiguation(band)
    if band == 'GALEXFUV' or band == 'GALEXNUV': return cm.bone
    if band == 'SDSS_u' or band == 'SDSS_g' or band == 'SDSS_r' \
        or band == 'SDSS_i' or band == 'SDSS_z': return cm.gray
    if band == 'J' or band == 'H' or band == 'Ks': return cm.gray
    if band == 'WISEW1' or band == 'WISEW2' \
        or band == 'IRAC1' or band == 'IRAC2' \
        or band == 'IRAC3' or band == 'IRAC4': return cm.pink
    if band == 'WISEW3' or band == 'WISEW4' \
        or band == 'MIPS24': return cm.Reds
    if band == 'MIPS70' or band == 'PACS70' or band == 'PACS100' \
        or band == 'PACS160' or band == 'MIPS160' or band == 'SPIRE250' \
        or band == 'SPIRE350' or band == 'SPIRE500': return cm.gist_heat

def round_sig(x, sig = 2):
    try: return round(x, sig-((np.floor(np.log10(abs(x))))-1).astype('int'))
    except: return np.round(x, sig)

def round_arr(x, sig = 2):
    return np.array([round_sig(el, sig) for el in x])

def linear_function(B, x):
    '''
    Linear function y = m*x + b
    '''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]

def pairwise(iterable):
    from itertools import izip
    "s -> (s0, s1), (s2, s3), (s4, s5), ..."
    a = iter(iterable)
    return izip(a, a)

def download_data(galaxy_name, bands_to_download, download_directory, processes = 10):
    from itertools import repeat
    import multiprocessing
    subprocess.call('mkdir '+download_directory.split('/')[0], shell = True)
    subprocess.call('mkdir '+download_directory.split('/')[0]+'/'+download_directory.split('/')[1], shell = True)
    actual_path = os.getcwd()
    os.chdir(download_directory)
    pool = multiprocessing.Pool()
    with multiprocessing.Pool(processes=processes) as pool:
        func = zip(bands_to_download, repeat(galaxy_name))
        pool.starmap(parallel_download, func)
    os.chdir(actual_path)
    return 'All maps have been stored in ' + download_directory

def parallel_download(band, galaxy_name):
    instrument = str(band.split('_', 1)[0])
    if instrument == 'Spitzer':
        link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'.fits.gz&instrument='+instrument+'"'
        command = 'wget -cO -'+link_path+' > '+galaxy_name+'_'+band+'.fits.gz'
        subprocess.call(command, shell = True)
        subprocess.call('gunzip '+galaxy_name+'_'+band+'.fits.gz', shell = True)
    else:
        link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'.fits&instrument='+instrument+'"'
        command = 'wget -cO -'+link_path+' > '+galaxy_name+'_'+band+'.fits'
        subprocess.call(command, shell = True)
    print('Downloaded '+galaxy_name+'_'+band)
    return

def download_data_errors(galaxy_name, bands_to_download, download_directory):
    # Sono poche, pesano una sega, non c'è bisogno di ||.
    subprocess.call('mkdir '+download_directory.split('/')[0], shell = True)
    subprocess.call('mkdir '+download_directory.split('/')[0]+'/'+download_directory.split('/')[1], shell = True)
    actual_path = os.getcwd()
    os.chdir(download_directory)
    for band in bands_to_download:
        instrument = str(band.split('_', 1)[0])
        if instrument == 'Spitzer':
            link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'_Error.fits.gz&instrument='+instrument+'"'
            command = 'wget -cO -'+link_path+' > '+band+'_Error.fits.gz'
        else:
            link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'_Error.fits&instrument='+instrument+'"'
            command = 'wget -cO -'+link_path+' > '+band+'_Error.fits'
        print('Downloading '+galaxy_name+'_'+band+' error map.')
        subprocess.call(command, shell = True)
    os.chdir(actual_path)
    return 'All error maps have been stored in '+download_directory
    
def run_CAAPR(galaxy_name):
    actual_path = os.getcwd()
    os.chdir('Caapr')
    subprocess.call('python Caapr_'+galaxy_name[3:]+'.py', shell = True)
    os.chdir(actual_path)
    return 'Done!'
    
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    QUI STO COPIANDO DA CHRIS CLARK
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

def Downsample(myarr, factor, estimator=np.nanmean):
    '''
    Keflavich function to downsample an array
    Args: Array to downsample, downsampling factor, and estiamtor
    Returns: Downsampled array
    '''
    ys,xs = myarr.shape
    crarr = myarr[:ys-(ys % int(factor)),:xs-(xs % int(factor))]
    dsarr = estimator( np.concatenate([[crarr[i::factor,j::factor]
        for i in range(factor)]
        for j in range(factor)]), axis=0)
    return dsarr

def SigmaClip(values, tolerance=0.001, median=False, sigma_thresh=3.0, no_zeros=False):
    '''
    Function to perform a sigma clip upon a set of values
    Args: Array of values, convergence tolerance, state if median instead of mean should be used for clip centrepoint,
          clipping threshold, boolean for whether sigma of zero can be accepted
    Returns: List containing the clipped standard deviation, the average, and the values themselves

    '''
    # Remove NaNs from input values
    values = np.array(values)
    values = values[ np.where(np.isnan(values)==False) ]
    values_original = np.copy(values)

    # Continue loop until result converges
    diff = 10E10
    while diff>tolerance:

        # Assess current input iteration
        if median == False:
            average = np.mean(values)
        elif median == True:
            average = np.median(values)
        sigma_old = np.std(values)

        # Mask those pixels that lie more than 3 stdev away from mean
        check = np.zeros([len(values)])
        check[ np.where( values>(average+(sigma_thresh*sigma_old)) ) ] = 1
        check[ np.where( values<(average-(sigma_thresh*sigma_old)) ) ] = 1
        values = values[ np.where(check<1) ]

        # Re-measure sigma and test for convergence
        sigma_new = np.std(values)
        diff = abs(sigma_old-sigma_new) / sigma_old

    # Perform final mask
    check = np.zeros([len(values)])
    check[ np.where( values>(average+(sigma_thresh*sigma_old)) ) ] = 1
    check[ np.where( values<(average-(sigma_thresh*sigma_old)) ) ] = 1
    values = values[ np.where(check<1) ]

    # If required, check if calculated sigma is zero
    if no_zeros==True:
        if sigma_new==0.0:
            sigma_new = np.std(values_original)
            if median==False:
                average = np.mean(values)
            elif median==True:
                average = np.median(values)

    # Return results
    return [sigma_new, average, values]

def LogError(value, error):
    '''
    # New function to convert an uncertainty to log space
    # Args: Value, uncertainty
    # Returns: Logarithmic uncertainty (up, down, and mean bw up and down)
    '''
    value, error = np.array(value), np.array(error)
    frac = 1.0 + (error/value)
    error_up = value * frac
    error_down = value / frac
    log_error_up = np.abs( np.log10(error_up) - np.log10(value) )
    log_error_down = np.abs( np.log10(value) - np.log10(error_down) )
    return log_error_up, log_error_down, 0.5*(log_error_up+log_error_down)

def Most_Common(lst):
    from collections import Counter
    data = Counter(lst)
    return data.most_common(1)[0][0]

######

def old_download_data(galaxy_name, bands_to_download, download_directory):
    # Ancora non è parallelizzato.
    subprocess.call('mkdir '+download_directory.split('/')[0], shell = True)
    subprocess.call('mkdir '+download_directory.split('/')[0]+'/'+download_directory.split('/')[1], shell = True)
    actual_path = os.getcwd()
    os.chdir(download_directory)
    for band in bands_to_download:
        instrument = str(band.split('_', 1)[0])
        if instrument == 'Spitzer':
            link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'.fits.gz&instrument='+instrument+'"'
            command = 'wget -cO -'+link_path+' > '+galaxy_name+'_'+band+'.fits.gz'
            print('Downloading '+galaxy_name+'_'+band)
            subprocess.call(command, shell = True)
            subprocess.call('gunzip '+galaxy_name+'_'+band+'.fits.gz', shell = True)
        else:
            link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'.fits&instrument='+instrument+'"'
            command = 'wget -cO -'+link_path+' > '+galaxy_name+'_'+band+'.fits'
            print('Downloading '+galaxy_name+'_'+band)
            subprocess.call(command, shell = True)
    os.chdir(actual_path)
    return 'Done!'
