import sys, os
import pandas as pd
from astropy.io import fits
from tqdm import tqdm, trange
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm, gridspec
from astropy.visualization import (MinMaxInterval, LogStretch, ImageNormalize)
from photutils import SkyRectangularAperture    
from astropy.coordinates import SkyCoord

import GenericUsefulScripts as GUS

fits_base = GUS.FitsUtils('/data/dustpedia/DustFolder/NGC4535/_ReducedMaps/SPIRE_350.fits')
phot_df = pd.read_csv('GOOD_photometries_table.csv', sep = '\t')
photcols = ['GALEX_FUV', 'GALEX_NUV', 'SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', \
            '2MASS_J', '2MASS_H', '2MASS_Ks', 'WISE_3.4', 'Spitzer_3.6',  'Spitzer_4.5', \
            'WISE_4.6', 'WISE_12', 'WISE_22', 'PACS_70', 'PACS_100', 'PACS_160', 'SPIRE_250', 'SPIRE_350']
for band in photcols:
    phot_df[band+'_SNR'] = phot_df[band]/phot_df[band+'_ERROR']

lam_observed = [0.1528, 0.2271, 0.3551, 0.4686, 0.6166, 0.7480, 0.8932, \
           1.25, 1.65, 2.16, 3.40, 3.56, 4.51, 4.60, \
           12.0, 22.0, 70.0, 100., 160., 250., 350.]

 aperture_side = 8*u.arcsec

for row in tqdm(phot_df.iterrows()):
    plt.close('all')
    if os.path.exists('seds/{0:04}_sed.png'.format(row[0])): continue

    SED_observed = row[1][photcols]
    SED_snr= row[1][[photcol + '_SNR' for photcol in photcols]]

    coord = SkyCoord(row[1].ra*u.deg, row[1].dec*u.deg, frame = 'icrs')
    aperture = SkyRectangularAperture(coord, w=aperture_side, h=aperture_side, theta = 0*u.deg)

    SED = plt.figure(figsize=(15,5))
    
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.8, 1]) 
    ax1 = plt.subplot(gs[0]) # SED with observed datapoints
    ax2 = plt.subplot(gs[1], projection = fits_base.wcs) # Galaxy with selected aperture

    colors = ['tab:red' if t < 2 else \
              ('orange' if t > 2 and t < 3 else \
                ('gold' if t > 3 and t < 5 else \
                ('lightgreen' if t > 5 and t < 8 else \
                ('tab:green' if t > 8 and t < 10 else 'darkgreen')))) for t in SED_snr]

    ax1.scatter(lam_observed, SED_observed, s = 40, color = colors, linestyle = 'None', marker = 'o')
    ax1.set_xscale('log'), ax1.set_yscale('log')
    ax1.set_xlabel(r'Wavelength ($\mu$m)'), ax1.set_ylabel(r'Flux (Jy)')
    ymin, ymax = 1E-7, 1E2 # Fixed y axis
    ax1.set_xlim(1E-1, 7E2), ax1.set_ylim(ymin, ymax)
    
    ax1.plot(0, 0, color = 'tab:red', linewidth = 0, ms = 10, marker = 'o', label = 'SNR < 2')
    ax1.plot(0, 0, color = 'orange', linewidth = 0, ms = 10, marker = 'o', label = '2 < SNR < 3')
    ax1.plot(0, 0, color = 'gold', linewidth = 0, ms = 10, marker = 'o', label = '3 < SNR < 5')
    ax1.plot(0, 0, color = 'lightgreen', linewidth = 0, ms = 10, marker = 'o', label = '5 < SNR < 8')
    ax1.plot(0, 0, color = 'tab:green', linewidth = 0, ms = 10, marker = 'o', label = '8 < SNR < 10')
    ax1.plot(0, 0, color = 'darkgreen', linewidth = 0, ms = 10, marker = 'o', label = 'SNR > 10')
    ax1.legend(loc = 'upper left', frameon=False)
    
    # Create an ImageNormalize object
    norm = ImageNormalize(fits_base.signal, interval=MinMaxInterval(), stretch=LogStretch())
    ax2.imshow(fits_base.signal, origin = 'lower', interpolation = 'nearest', cmap = cm.inferno, norm = norm)
    
    ra_x, dec_y = aperture.positions.to_pixel(fits_base.wcs)[0], aperture.positions.to_pixel(fits_base.wcs)[1]
    pixel_scale = fits_base.get_pixel_scale()*u.deg
    term = (.5*aperture_side/pixel_scale.to('arcsec')).value
    rectangle = plt.Rectangle((ra_x-term, dec_y+term),
                              width = aperture_side/pixel_scale, height = aperture_side/pixel_scale,
                              color='w', linewidth = 3.0, fill = True)
    
    ax2.add_artist(rectangle)
    
    ax2.set_xlim(175, 275), ax2.set_ylim(180, 280)
    
    plt.tight_layout()
    SED.savefig('seds/{0:04}_sed.png'.format(row[0]), bbox_inches = 'tight')

# Generate .gif
import imageio
images = [imageio.imread('seds/'+filename) for filename in sorted(os.listdir('seds/')) if not filename.startswith('.')]
imageio.mimsave('seds.gif', images, duration = 0.1)

import moviepy.editor as mp
clip = mp.VideoFileClip("seds.gif")
clip.write_videofile("seds.mp4")