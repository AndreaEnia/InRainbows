#!/usr/bin/python
# -*- coding: latin-1 -*-

# Function to find distance 
def shortest_distance(xp, yp, coeff_m, coeff_q):
    '''
    Approach: The distance (i.e shortest distance) from a given point to a line is the perpendicular distance from that point to the given line.
    The equation of a line in the plane is given by the equation ax + by + c = 0, where a, b and c are real constants.
    -by = ax + c
    -1y = mx + q
    The coordinate of the point is (xp, yp)
    The formula for distance between a point and a line in 2-D is given by:
    Distance = (| a*x1 + b*y1 + c |) / (sqrt( a*a + b*b))
    '''
    import math 
    return -((coeff_m*xp - yp + coeff_q)) / (math.sqrt(coeff_m**2 + 1)) 

def renormalize(array, a, b):
    return (b-a)*((array - array.min())/array.ptp()) + a

def density_scatter(m1, m2):
    from scipy.stats import gaussian_kde
    xy = np.vstack([m1,m2])
    z = gaussian_kde(xy)(xy)
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    return m1[idx], m2[idx], z[idx]

def density_estimation(m1, m2):
    from scipy.stats import gaussian_kde
    import numpy as np
    xmin = m1.min()
    xmax = m1.max()
    ymin = m2.min()
    ymax = m2.max()
    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]                                                     
    positions = np.vstack([X.ravel(), Y.ravel()])                                                       
    values = np.vstack([m1, m2])                                                                        
    kernel = gaussian_kde(values)                                                                 
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z

def read_data(signal_path):
    '''
    Read .fits data
    INPUT:
        path of the desired signal
    OUTPUT:
        (in order): Signal, header, pixel scale and world coordinate system from header
    '''
    import numpy as np
    from astropy.io import fits
    from astropy import units as u
    from image_tools import array_positions
    from astropy.wcs import WCS
    signal = np.nan_to_num(fits.getdata(signal_path))
    hdr = fits.getheader(signal_path)
    try: pixel_scale = round(np.abs(hdr['CD1_1'])*3600, 5)*u.arcsec
    except: pixel_scale = round(np.abs(hdr['CDELT1'])*3600, 5)*u.arcsec
    wcs = WCS(hdr)
    return signal, hdr, pixel_scale, wcs

def read_data_full(signal_path, noise_path):
    '''
    Read .fits data
    INPUT:
        path of the desired signal, path of the eventual noise
    OUTPUT:
        (in order): signal, noise, SNR, hdr, wcs, nx, ny, pixel_scale, extent_im, index_pixel_x, index_pixel_y, x_arcsec, y_arcsec, ra_ref, dec_ref
    '''
    from astropy.io import fits
    import numpy as np
    from image_tools import array_positions
    from astropy.wcs import WCS

    signal = np.nan_to_num(fits.getdata(signal_path))
    try: noise  = np.nan_to_num(fits.getdata(noise_path))
    except: noise = 1.0
    SNR = signal/noise
    hdr = fits.getheader(signal_path)
    try: pixel_scale = round(abs(hdr['CDELT1'])*3600, 5)*u.arcsec
    except: pixel_scale = round(abs(hdr['CD1_1'])*3600, 5)*u.arcsec
    nx, ny = hdr['naxis1'], hdr['naxis2']
    extent_im = [-nx*pixel_scale.value/2, nx*pixel_scale.value/2, -ny*pixel_scale.value/2, ny*pixel_scale.value/2]
    index_pixel_x, index_pixel_y, x_arcsec, y_arcsec = array_positions(nx, ny, pixel_scale)
    x_cen_arcsec, y_cen_arcsec = nx*pixel_scale/2.0, ny*pixel_scale/2.0
    ra_ref, dec_ref = hdr['crval1']*u.deg, hdr['crval2']*u.deg
    wcs = WCS(hdr)
    return signal, noise, SNR, hdr, wcs, nx, ny, pixel_scale, extent_im, index_pixel_x, index_pixel_y, x_arcsec, y_arcsec, ra_ref, dec_ref

def find_nearest(array, value):
    '''
    Find the argument in [array] with the value nearest to [value]
    INPUT:
        Array, value
    OUTPUT:
        Most similiar argument
    '''
    import numpy as np
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()

def evaluate_SFRs_wErrors(temp_res, lam_observed):
    '''
    Evaluate SFRs from magphys .fit and .sed results.
    INPUT:
        temp_res, what the magphys_read module gives back
        lam_observed, the wavelenghts of the observed fluxes
    OUTPUT:
        SFR_IR: Dust-obscured SFR from 8-1000 micrometer integration of the SED, SFR_IR = 2.8E-44*L_IR(8-1000) [erg/s], Kennicutt (1998) 
        SFR_UV: Direct UV SFR from 2800 Angstrom observation, SFR_UV = 1.0E-28*Lnu(2800) [erg/s/Hz], Bell et al. (2005)
                Both rescaled for Chabrier (2003) IMF
        SFR_TOTAL  Simple sum of the two terms
        eSFR_IR: error on the IR-SFR, evaluated by integrating the SED obtained summing to the best-fut one the error from the observed flux SNR.
        eSFR_UV: error on the UV-SFR, evaluated from the nearest 2800 Angstrom measurement SNR (that is, GALEX NUV)
        eSFR_TOTAL: 
    Written: A. Enia (May19)
    '''
    import numpy as np
    from astropy import units as u
    from astropy import constants as const
    lambda_model = (temp_res.sed_model_logwaves).to('Angstrom')  # Wavelenghts of the SED model, in A
    nu_model = lambda_model.to(u.Hz, equivalencies=u.spectral()) # Frequencies of the SED model, in Hz
    # IR
    SED_attenuated = 10**(temp_res.sed_model[:,1])*u.L_sun/u.Angstrom # Attenuated SED in original magphys units (Lsol/A)
    y = (SED_attenuated*(lambda_model**2/const.c.to('Angstrom/s'))).value*u.L_sun/u.Hz # Attenuated SED in Lsol/Hz
    y = y.to('erg/s/Hz') # Attenuated SED in erg/s/Hz
    ok_integr = np.logical_and(lambda_model.to('micron').value >= 8.0, lambda_model.to('micron').value <= 1000.0)
    Lir_integr = (np.abs(np.trapz(nu_model[ok_integr],y[ok_integr]))) # Integrated 8-1000 IR luminosity in erg/s
    #SFR_IR = 2.8E-44*Lir_integr.value # IR-SFR in Msun/yr
    SFR_IR = 2.64E-44*Lir_integr.value # IR-SFR in Msun/yr
    # UV
    ok_280 = find_nearest(lambda_model.value, 2800)
    #SFR_UV = 1.0E-28*y[ok_280].value # UV-SFR in Msun/yr
    SFR_UV = 0.82E-28*y[ok_280].value # UV-SFR in Msun/yr

    # IR error, I integrate the SED obtained summing the same error of the observed signal to the best fit SED.
    yerr = np.copy(y)
    for lam, i_lam in zip(lambda_model, range(len(lambda_model))):
        ok_nearest = find_nearest(lam.value, lam_observed.to('Angstrom').value)
        SNR = temp_res.obs_flux[ok_nearest]/temp_res.obs_flux_err[ok_nearest]
        if SNR < 1: SNR = (0.1)**(-1) # Realistic, 10% error
        error = (y[i_lam]/SNR).value
        yerr[i_lam] = y[i_lam].value + error
        
    eSFR_IR = 2.64E-44*(np.abs(np.trapz(nu_model[ok_integr],yerr[ok_integr]))).value - SFR_IR # IR-SFR in Msun/yr
    # UV error, I take the same error associated to the nearest 2800 Angstrom measure, that is GALEX NUV [1].
    ok_nearest = find_nearest(lam.value, 2800)
    eSFR_UV = 0.82E-28*y[ok_280].value*(temp_res.obs_flux_err[ok_nearest]/temp_res.obs_flux[ok_nearest])

    if eSFR_IR == 0: eSFR_IR = 0.1*SFR_IR # Fiduciary 10% error of the measure
    if eSFR_UV == 0: eSFR_UV = 0.1*SFR_UV # Fiduciary 10% error of the measure
    return SFR_IR, SFR_UV, SFR_IR+SFR_UV, eSFR_IR, eSFR_UV, np.sqrt(eSFR_IR**2 + eSFR_UV**2)

def evaluate_SFRs(temp_res, UV = 'NUV'):
    '''
    Evaluate SFRs from magphys .fit and .sed results.
    INPUT:
        temp_res, what the magphys_read module gives back
        lam_observed, the wavelenghts of the observed fluxes
    OUTPUT:
        SFR_IR: Dust-obscured SFR from 8-1000 micrometer integration of the SED, SFR_IR = 2.8E-44*L_IR(8-1000) [erg/s], Kennicutt (1998) 
        SFR_UV: Direct UV SFR from 2800 Angstrom observation, SFR_UV = 1.0E-28*Lnu(2800) [erg/s/Hz], Bell et al. (2005)
                Both rescaled for Chabrier (2003) IMF
        SFR_TOTAL  Simple sum of the two terms
        eSFR_IR: error on the IR-SFR, evaluated by integrating the SED obtained summing to the best-fut one the error from the observed flux SNR.
        eSFR_UV: error on the UV-SFR, evaluated from the nearest 2800 Angstrom measurement SNR (that is, GALEX NUV)
        eSFR_TOTAL: 
    Written: A. Enia (May19)
    '''
    import numpy as np
    from astropy import units as u
    from astropy import constants as const
    lambda_model = (temp_res.sed_model_logwaves).to('Angstrom')  # Wavelenghts of the SED model, in A
    nu_model = lambda_model.to(u.Hz, equivalencies=u.spectral()) # Frequencies of the SED model, in Hz
    SED_attenuated = 10**(temp_res.sed_model[:,1])*u.L_sun/u.Angstrom # Attenuated SED in original magphys units (Lsol/A)
    y = (SED_attenuated*(lambda_model**2/const.c.to('Angstrom/s'))).value*u.L_sun/u.Hz # Attenuated SED in Lsol/Hz
    y = y.to('erg/s/Hz') # Attenuated SED in erg/s/Hz
    # IR
    ok_integr = np.logical_and(lambda_model.to('micron').value >= 8.0, lambda_model.to('micron').value <= 1000.0)
    Lir_integr = (np.abs(np.trapz(nu_model[ok_integr],y[ok_integr]))) # Integrated 8-1000 IR luminosity in erg/s
    #SFR_IR = 2.8E-44*Lir_integr.value # IR-SFR in Msun/yr
    SFR_IR = 2.64E-44*Lir_integr.value # IR-SFR in Msun/yr
    ## UV
    if UV == 'NUV':
        ok_280 = find_nearest(lambda_model.value, 2800)
        SFR_UV = 0.82E-28*y[ok_280].value # UV-SFR in Msun/yr
        #SFR_UV = 1.0E-28*y[ok_280].value # UV-SFR in Msun/yr
    elif UV == 'FUV':
        SFR_UV = SFRs_UV_1500(lambda_model, y)
    
    return SFR_IR, SFR_UV, SFR_IR+SFR_UV

def evaluate_SFRs_Bell2005(temp_res):
    '''
    Evaluate SFRs from magphys .fit and .sed results.
    INPUT:
        temp_res, what the magphys_read module gives back
        lam_observed, the wavelenghts of the observed fluxes
    OUTPUT:
        SFR_IR: Dust-obscured SFR from L_IR, the 8-1000 micrometer integration of the SED
        SFR_UV: Direct UV SFR from 2800 Angstrom observation, L_UV = 1.5*nuLnu,2800
        SFR_TOTAL: from Bell et al., 2005, and Whitaker et al., 2014, SFR [Msol/yr] = 1.09E-10(L_IR+2.2L_UV)
                        All rescaled for Chabrier (2003) IMF
    Written: A. Enia (Dec19)
    '''
    import numpy as np
    from astropy import units as u
    from astropy import constants as const
    lambda_model = (temp_res.sed_model_logwaves).to('Angstrom')  # Wavelenghts of the SED model, in A
    nu_model = lambda_model.to(u.Hz, equivalencies=u.spectral()) # Frequencies of the SED model, in Hz
    # IR
    SED_attenuated = 10**(temp_res.sed_model[:,1])*u.L_sun/u.Angstrom # Attenuated SED in original magphys units (Lsol/A)
    y = (SED_attenuated*(lambda_model**2/const.c.to('Angstrom/s'))).value*u.L_sun/u.Hz # Attenuated SED in Lsol/Hz
    y = y.to('erg/s/Hz') # Attenuated SED in erg/s/Hz
    ok_integr = np.logical_and(lambda_model.to('micron').value >= 8.0, lambda_model.to('micron').value <= 1000.0)
    Lir_integr = (np.abs(np.trapz(nu_model[ok_integr],y[ok_integr]))) # Integrated 8-1000 IR luminosity in erg/s
    Lir_integr = Lir_integr.to('Lsun')
    SFR_IR = 1.09E-10*Lir_integr.value
    # UV
    ok_2800 = find_nearest(lambda_model.value, 2800)
    nu_2800 = (2800 * u.Angstrom).to(u.Hz, equivalencies=u.spectral())
    L_UV = (1.5*nu_2800*y[ok_2800]).to('Lsun')
    SFR_UV = 2.2*1.09E-10*L_UV.value
    return SFR_IR, SFR_UV, SFR_IR+SFR_UV

def SFRs_UV_1500(lambda_model, y):
    '''
    Evaluate UV unobscured SFR from magphys .fit and .sed results.
    INPUT:
        temp_res, what the magphys_read module gives back
    OUTPUT:
        SFR_UV: Direct UV SFR from 1500 Angstrom observation, rescaled in Chabrier: L_UV = 0.63*1.13E-28*nuLnu,1500
        o 0.63*1.4. Non si capisce un cazzo dio maiale.
    Written: A. Enia (Dec19)
    '''
    import numpy as np
    from astropy import units as u
    from astropy import constants as const
    ok_1500 = find_nearest(lambda_model.value, 1500)
    nu_1500 = (1500 * u.Angstrom).to(u.Hz, equivalencies=u.spectral())
    return 0.63*1.4E-28*y[ok_1500].value 

def read_tables(path, chisq_limit = 25):
    '''
    Read tables with physical properties.
    INPUT:
        path
    OUTPUT:
        Id_Ap, bestfit_chisq, GOOD_chisq, SigmaSFR, SigmaMstar, eSigmaSFR, A_IRX, T_cold, cendist, cosmo_factor, \
        SigmaSFR_integrated, SigmaMstar_integrated

    Written: A. Enia (May19)
    Modified: A. Enia (Oct19)
    '''

    import pandas as pd
    import numpy as np

    print('Reading '+path)
    table = pd.read_csv(path, delimiter = ',', dtype={'Id. Aperture': object})
    Id_Ap, cendist = table['Id. Aperture'].values, table['cendist_physical'].values
    ra, dec = table['ra_apertures'].values, table['dec_apertures'].values
    try: bestfit_chisq = table['Bestfit_Chisq'].values
    except: bestfit_chisq = table['Bestfit Chisq'].values
    GOOD_chisq = np.where(bestfit_chisq <= chisq_limit)
    SFR_total, Mstar = table['SFR_total'].values, table['Mstar'].values
    SigmaSFR, SigmaMstar = table['SigmaSFR'].values, table['SigmaMstar'].values
    A_IRX, T_C_ISM, T_W_BC = table['A_IRX'].values, table['T_C_ISM'].values, table['T_W_BC'].values
    return Id_Ap, ra, dec, cendist, bestfit_chisq, GOOD_chisq, SFR_total, Mstar, SigmaSFR, SigmaMstar, A_IRX, T_C_ISM, T_W_BC

def read_integrated_tables(path):
    '''
    Read tables with physical properties from INTEGRATED run
    INPUT:
        path
    OUTPUT:
        Id_Ap, bestfit_chisq, physical_scale, Mstar, SFR_magphys, SFR_IR, SFR_UV, SFR_total, SigmaSFR, SigmaMstar,
        SigmaSFR_magphys, A_IRX, Ldust, Mdust, T_C_ISM, T_W_BC
    Written: A. Enia (May19)
    Modified: A. Enia (Oct19)
    '''

    import pandas as pd
    import numpy as np

    print('Reading '+path)
    table = pd.read_csv(path, delimiter = '\t', dtype={'id_ap': object})
    Id_Ap, bestfit_chisq = table['id_ap'].values, table['bestfit_chisq'].values
    physical_scale = table['physical_scale'].values
    Mstar = table['Mstar'].values    
    SFR_magphys = table['SFR'].values
    SFR_IR, SFR_UV, SFR_total = table['SFR_IR'].values, table['SFR_UV'].values, table['SFR_total'].values
    try: eSFR = table['eSFR_total'].values
    except: eSFR = 0
    SigmaSFR, SigmaMstar = table['SigmaSFR'].values, table['SigmaMstar'].values
    SigmaSFR_magphys = table['SigmaSFR_magphys'].values
    A_IRX = table['SFR_IR'].values/table['SFR_total'].values
    Ldust, Mdust = table['Ldust'].values, table['Mdust'].values
    T_C_ISM, T_W_BC = table['T_C_ISM'].values, table['T_W_BC'].values
    return Id_Ap, bestfit_chisq, physical_scale, Mstar, SFR_magphys, SFR_IR, SFR_UV, SFR_total, eSFR, \
            SigmaSFR, SigmaMstar, SigmaSFR_magphys, A_IRX, Ldust, Mdust, T_C_ISM, T_W_BC

def round_sig(x, sig=2):
    import numpy as np
    return round(x, sig-((np.floor(np.log10(abs(x))))-1).astype('int'))

def round_arr(x, sig=2):
    import numpy as np
    return np.array([round_sig(el, sig) for el in x])

def linear_function(B, x):
    '''Linear function y = m*x + b'''
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

def download_data(galaxy_name, working_bands, download_directory):
    # Ancora non è parallelizzato.
    import subprocess
    import os
    actual_path = os.getcwd()
    os.chdir(download_directory)
    for band in working_bands:
        instrument = str(band.split('_', 1)[0])
        if instrument == 'Spitzer':
            link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'.fits.gz&instrument='+instrument+'"'
            command = 'wget -cO -'+link_path+' > '+galaxy_name+'_'+band+'.fits.gz'
        else:
            link_path = ' "http://dustpedia.astro.noa.gr/Data/GetImage?imageName='+galaxy_name+'_'+band+'.fits&instrument='+instrument+'"'
            command = 'wget -cO -'+link_path+' > '+galaxy_name+'_'+band+'.fits'
        print('Downloading '+galaxy_name+'_'+band)
        subprocess.call(command, shell = True)
    os.chdir(actual_path)
    return 'Done!'

def run_CAAPR(galaxy_name):
    import os
    import subprocess
    actual_path = os.getcwd()
    os.chdir('Caapr')
    subprocess.call('python Caapr_'+galaxy_name[3:]+'.py', shell = True)
    os.chdir(actual_path)
    return 'Done!'
    