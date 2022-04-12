#!/usr/bin/python
# -*- coding: latin-1 -*-
import sys, os, subprocess
import GenericUsefulScripts as GUS
import PhotometryScripts as PhotoScripts
sys.path.insert(0, '../')
import magphys_read
import DustpediaScripts as DS
import MainSequences as MS
import numpy as np
import pandas as pd
from tqdm import tqdm, trange
from astropy import units as u
from astropy.cosmology import Planck15 as cosmo
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import ascii, fits
from astropy.constants import c as v_lux
from photutils import SkyRectangularAperture    

def read_magphys_v0(GalProp, run_type, cell_side, print_totals = True):
    
    # =================
    # Galaxy properties
    arcsec_per_kpc = cosmo.arcsec_per_kpc_comoving(GalProp.z_source)
    kpc_per_arcmin = cosmo.kpc_comoving_per_arcmin(GalProp.z_source)
    magphys_results_folder = '../'+GalProp.galaxy_name+'/magphys_'+run_type+'/run_'+run_type
    
    # =====================
    # Paths for .sed and .fit
    sed_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.sed')])
    fit_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.fit')])
    # =====================
    
    # =================
    # Read the apertures
    if cell_side.unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_side.unit == 'arcmin': cell_side = cell_side.to('arcsec')
    photometry_folder_path = '../'+GalProp.galaxy_name+'/Photometry_'+run_type+'/'
    data = pd.read_csv(photometry_folder_path+'FINAL.csv')
    ra_apertures, dec_apertures = data['ra'].values*u.deg, data['dec'].values*u.deg   
    coords = SkyCoord(ra_apertures, dec_apertures, frame = 'icrs')
    if cell_side.unit == 'kpc': cell_side *= arcsec_per_kpc        
    list_apertures = [SkyRectangularAperture(coords[i_ap], w=cell_side, h=cell_side, theta = 0*u.deg) \
                      for i_ap in trange(len(ra_apertures), ascii = True, desc = 'Reading apertures')]
    # =================

    # ==============
    # Read results from Magphys, then save everything in a table.
    id_ap, bestfit_chisq = [], []
    fmuIR, fmuSFH = [], []
    Ldust, Mdust, = [], []
    Mstar, Mstar_m1s, Mstar_p1s = [], [], []  # m1s = Minus 1 sigma, p1s = plus 1 sigma
    T_C_ISM, T_W_BC = [], []
    SFR_magphys, sSFR_magphys, tauv, tvism = [], [], [], []
    SFR_IR, SFR_UV, SFR_total = [], [], []
    eSFR_IR, eSFR_UV, eSFR_total = [], [], []
    xi_C_tot, xi_MIR_tot, xi_PAH_tot, xi_W_tot = [], [], [], []
    tform, age_wm = [], []
    
    for sed, fit in tqdm(zip(sed_path, fit_path), ascii = True, desc = 'Reading MAGPHYS results'):
        i_file = sed[-8:-4] # Ricordati di cambiare da sed[14:17] a sed[14:-4] quando arrivano i nuovi risultati.
        id_ap.append(i_file)
        try: results = magphys_read.MagphysOutput(fit, sed)
        except:
            print('There has been an error in reading this results:')
            print(sed, fit)
            emergency_path = '00'+str(int(i_file)-1)
            results = magphys_read.MagphysOutput(fit_path[0][:-9]+emergency_path+'.fit', sed_path[0][:-9]+emergency_path+'.sed')
            print('So, I am reading the next one, just to go on with the code. Check what went wrong.')
        # Evaluate SFRs
        SFRs = DS.evaluate_SFRs(results, UV = 'FUV')
        SFR_IR.append(SFRs[0]), SFR_UV.append(SFRs[1]), SFR_total.append(SFRs[2])
        # Append the rest
        bestfit_chisq.append(results.bestfit_chi2)
        fmuIR.append(results.fmuIR), fmuSFH.append(results.fmuSFH)
        Ldust.append(results.Ldust), Mdust.append(results.Mdust)
        Mstar.append(results.Mstar)
        sigm1, sigp1 = 10**results.marginal_percentiles['Mstars'][1], 10**results.marginal_percentiles['Mstars'][3]
        Mstar_m1s.append(sigm1), Mstar_p1s.append(sigp1)
        T_C_ISM.append(results.T_C_ISM), T_W_BC.append(results.T_W_BC)
        SFR_magphys.append(results.SFR), sSFR_magphys.append(results.sSFR)
        tauv.append(results.tauv), tvism.append(results.tvism)
        xi_C_tot.append(results.xi_C_tot), xi_MIR_tot.append(results.xi_MIR_tot)
        xi_PAH_tot.append(results.xi_PAH_tot), xi_W_tot.append(results.xi_W_tot)
        tform.append(results.tform), age_wm.append(results.age_wm)
        
    # Distances from galaxy centre.
    cen_coord = SkyCoord(GalProp.ra, GalProp.dec, frame = 'icrs')
    sep = coords.separation(cen_coord)
    cendist_physical = sep.arcminute*u.arcmin*kpc_per_arcmin # In kpc
    
    # From list to array
    Mstar = np.array(Mstar)
    SFR_total = np.array(SFR_total)
    SFR_UV = np.array(SFR_UV)
    SFR_magphys = np.array(SFR_magphys)
    
    # Resolved properties
    aperture_side_in_kpc = (cell_side/arcsec_per_kpc).value
    SigmaSFR_magphys = SFR_magphys/aperture_side_in_kpc**2
    SigmaSFR = SFR_total/aperture_side_in_kpc**2
    try: eSigmaSFR = np.array(eSFR_total)/aperture_side_in_kpc**2
    except: pass
    SigmaMstar = np.array(Mstar)/aperture_side_in_kpc**2
    
    # Distance from MS
    answer = input('new or old pipeline?')
    if answer == 'old':
        REM_results = ascii.read('../Enia+20_oldpipe/MSfit_stuff/REM_sequence_params.txt')
    else:
        REM_results = ascii.read('../Enia+20_newpipe/MSfit_stuff/MSfit_pBp_params.txt')

    REM_m, REM_q = REM_results['m'][0], REM_results['q'][0]
    x, y = np.log10(SigmaMstar), np.log10(SigmaSFR)
    MS_distance = DS.shortest_distance(x, y, REM_m, REM_q)        
    
    table = pd.DataFrame()
    table['Id. Aperture'] = id_ap
    table['ra_apertures'] = ra_apertures
    table['dec_apertures'] = dec_apertures
    table['cendist_physical'] = cendist_physical                # Distance from center, in units of kpc
    table['Bestfit_Chisq'] = np.nan_to_num(bestfit_chisq)
    table['MS_distance'] = np.nan_to_num(MS_distance)           # Distance from the (best-fit) MS
    table['Mstar'] = np.nan_to_num(Mstar)                       # Stellar mass in units of solar masses
    table['SigmaMstar'] = np.nan_to_num(SigmaMstar)             # Mass density, in units of solar masses per kpc^2
    table['SFR_magphys'] = np.nan_to_num(SFR_magphys)
    table['SigmaSFR_magphys'] = SigmaSFR_magphys
    table['SFR_IR'] = np.nan_to_num(SFR_IR)                     # Star formation rate in units of solar masses per year
    try: table['eSFR_IR'] = np.nan_to_num(eSFR_IR)
    except: pass
    table['SFR_UV'] = np.nan_to_num(SFR_UV)                     # Star formation rate in units of solar masses per year
    try: table['eSFR_UV'] = np.nan_to_num(eSFR_UV)
    except: pass
    table['SFR_total'] = np.nan_to_num(SFR_total)               # Star formation rate in units of solar masses per year
    try: table['eSFR_total'] = np.nan_to_num(eSFR_total)
    except: pass
    table['SigmaSFR'] = np.nan_to_num(SigmaSFR)                 # SFR density, in units of solar masses per year per kpc^2
    try: table['eSigmaSFR'] = eSigmaSFR
    except: pass
    table['A_IRX'] = 2.5*np.log10(SFR_total/SFR_UV)             # Continuum attenuation (Puglisi+16).
    table['sSFR'] = np.nan_to_num(Mstar/SFR_total)              # Specific star formation rate in units of year^-1
    table['sSFR_magphys'] = np.nan_to_num(sSFR_magphys)         # Specific star formation rate in units of year^-1
    table['Mdust'] = np.nan_to_num(Mdust)                       # Dust mass in units of solar masses
    table['Ldust'] = np.nan_to_num(Ldust)                       # Dust luminosity in units of solar luminosities
    table['T_C_ISM'] = np.nan_to_num(T_C_ISM)                   # Equilibrium temperature of cold dust in the diffuse ISM
    table['T_W_BC'] = np.nan_to_num(T_W_BC)                     # Equilibrium temperature of warm dust in birth clouds
    table['tauv'] = np.nan_to_num(tauv)
    table['tvism'] = np.nan_to_num(tvism)
    table['fmuIR'] = np.nan_to_num(fmuIR)
    table['fmuSFH'] = np.nan_to_num(fmuSFH)
    table['xi_C_tot'] = np.nan_to_num(xi_C_tot)
    table['xi_MIR_tot'] = np.nan_to_num(xi_MIR_tot)
    table['xi_PAH_tot'] = np.nan_to_num(xi_PAH_tot)
    table['xi_W_tot'] = np.nan_to_num(xi_W_tot)
    table['tform'] = np.nan_to_num(tform)
    table['age_wm'] = np.nan_to_num(age_wm)
    table.to_csv(magphys_results_folder+'/'+GalProp.galaxy_name+'_results_'+run_type+'.dat', index = False, sep = ',')    
    print
    print('MAGPHYS results saved in '+magphys_results_folder+'/'+GalProp.galaxy_name+'_results_'+run_type+'.dat')
    print
    # ==============

    if print_totals:
        print('Total SFR:', np.sum(SFR_total))
        print('Total log stellar mass:', np.log10((np.sum(table['Mstar']))))
    return

def read_magphys(GalProp, run_type, cell_side, print_totals = True):
    
    # =================
    # Galaxy properties
    arcsec_per_kpc = cosmo.arcsec_per_kpc_comoving(GalProp.z_source)
    kpc_per_arcmin = cosmo.kpc_comoving_per_arcmin(GalProp.z_source)
    magphys_results_folder = '../'+GalProp.galaxy_name+'/magphys_'+run_type+'/run_'+run_type
    # =================
    
    # =====================
    # Paths for .sed and .fit
    sed_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.sed')])
    fit_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.fit')])
    # =====================
    
    # =================
    # Read the apertures
    if cell_side.unit == 'pc': cell_side = cell_side.to('kpc')
    if cell_side.unit == 'arcmin': cell_side = cell_side.to('arcsec')
    photometry_folder_path = '../'+GalProp.galaxy_name+'/Photometry_'+run_type+'/'
    data = pd.read_csv(photometry_folder_path+'GOOD_photometries_table.csv', sep = '\t')
    ra_apertures, dec_apertures = data['ra'].values*u.deg, data['dec'].values*u.deg   
    coords = SkyCoord(ra_apertures, dec_apertures, frame = 'icrs')
    if cell_side.unit == 'kpc': cell_side *= arcsec_per_kpc        
    list_apertures = [SkyRectangularAperture(coords[i_ap], w=cell_side, h=cell_side, theta = 0*u.deg) \
                      for i_ap in trange(len(ra_apertures), ascii = True, desc = 'Reading apertures')]
    # =================

    # ==============
    # Read results from Magphys, then save everything in a table.
    malloppone = []
    id_ap, bestfit_chisq = [], []
    SFR_IR, SFR_UV, SFR_total = [], [], []
    list_varnames = []
    for sed, fit in tqdm(zip(sed_path, fit_path), ascii = True, desc = 'Reading MAGPHYS results'):
        i_file = sed.split('/')[-1][:-4]
        try: results = magphys_read.MagphysOutput(fit, sed)
        except:
            print('There has been an error in reading this cell:')
            print(sed, fit)
            results = magphys_read.MagphysOutput('../scripts/emergency.fit', '../scripts/emergency.sed')
            print('So, I am reading the emergency one, just to go on with the code. Check what went wrong.')
        id_ap.append(i_file)
        bestfit_chisq.append(results.bestfit_chi2)
        malloppino = []
        for varname, varvalue in vars(results).items():
            if varname in results.bestfitparams_names:
                list_varnames.append(varname)
                malloppino.append(varvalue)
        # Evaluate SFRs
        SFRs = DS.evaluate_SFRs(results, UV = 'FUV')
        SFR_IR.append(SFRs[0]), SFR_UV.append(SFRs[1]), SFR_total.append(SFRs[2])
        malloppone.append(malloppino)
    malloppone = np.array(malloppone).T

    # Distances from galaxy centre.
    cen_coord = SkyCoord(GalProp.ra, GalProp.dec, frame = 'icrs')
    sep = coords.separation(cen_coord)
    cendist_physical = sep.arcminute*u.arcmin*kpc_per_arcmin # In kpc

    table = pd.DataFrame()
    table['Id. Aperture'] = id_ap
    table['ra_apertures'] = ra_apertures
    table['dec_apertures'] = dec_apertures
    table['cendist_physical'] = cendist_physical                # Distance from center, in units of kpc
    table['Bestfit_Chisq'] = np.nan_to_num(bestfit_chisq)
    table['SFR_IR'] = np.nan_to_num(SFR_IR)                     # Star formation rate in units of solar masses per year
    try: table['eSFR_IR'] = np.nan_to_num(eSFR_IR)
    except: pass
    table['SFR_UV'] = np.nan_to_num(SFR_UV)                     # Star formation rate in units of solar masses per year
    try: table['eSFR_UV'] = np.nan_to_num(eSFR_UV)
    except: pass
    table['SFR_total'] = np.nan_to_num(SFR_total)               # Star formation rate in units of solar masses per year
    try: table['eSFR_total'] = np.nan_to_num(eSFR_total)
    except: pass
    table['A_IRX'] = 2.5*np.log10(table['SFR_total']/table['SFR_UV'])             # Continuum attenuation (Puglisi+16).
    for column, column_name in zip(malloppone, list_varnames):
        if column_name == 'SFR': column_name = 'SFR_magphys'
        if column_name == 'sSFR': column_name = 'sSFR_magphys'
        table[column_name] = column
        
    # Resolved properties
    aperture_side_in_kpc = (cell_side/arcsec_per_kpc).value
    table['SigmaSFR_magphys'] = table['SFR_magphys']/aperture_side_in_kpc**2
    table['SigmaSFR'] = np.nan_to_num(table['SFR_total']/aperture_side_in_kpc**2) # SFR density, in units of solar masses per year per kpc^2
    try: table['eSigmaSFR'] = np.array(eSFR_total)/aperture_side_in_kpc**2
    except: pass

    table['SigmaMstar'] = table['Mstar']/aperture_side_in_kpc**2
    
    # Distance from MS
    x, y = np.log10(table['SigmaMstar']), np.log10(table['SigmaSFR'])
    answer = input('MS distance from E20 or M20? Or you do not want to measure it? (E/M/n)')
    if answer == 'n': pass
    elif answer == 'E':
        REM_m, REM_q = 0.82, -8.69
        table['MS_distance'] = np.nan_to_num(DS.shortest_distance(x, y, REM_m, REM_q))           # Distance from the (best-fit) MS
    elif answer == 'M':
        REM_m, REM_q = 0.73, -7.98
        table['MS_distance'] = np.nan_to_num(DS.shortest_distance(x, y, REM_m, REM_q))           # Distance from the (best-fit) MS
    table.to_csv(magphys_results_folder+'/'+GalProp.galaxy_name+'_results_'+run_type+'.dat', index = False, sep = ',')    
    print('************')
    print('MAGPHYS results saved in '+magphys_results_folder+'/'+GalProp.galaxy_name+'_results_'+run_type+'.dat')
    print('************')
    # =====================
    if print_totals:
        print('Total SFR:', np.sum(SFR_total))
        print('Total log Mstar:', np.log10((np.sum(table['Mstar']))))
    return


def clustering_removal(GalProp, run_type, s_size, eps_value = 0.09):
    from sklearn.cluster import DBSCAN
    from sklearn import metrics
    from sklearn.preprocessing import StandardScaler
    from sklearn.datasets import load_digits
    from sklearn.preprocessing import scale
    from sklearn.cluster import KMeans
    import matplotlib.pyplot as plt
    plt.close('all') # Clear all the previously opened graphic windows
    from matplotlib import cm

    galaxy_name = GalProp.galaxy_name
    ra_cen, dec_cen = GalProp.ra, GalProp.dec
    magphys_results_folder = '../'+galaxy_name+'/magphys_'+run_type+'/run_'+run_type
    
    sample_properties_path = '../DustPedia_Tables/DustPedia_HyperLEDA_Herschel.csv'
    results = GUS.ManageTable(galaxy_name, sample_properties_path, magphys_results_folder+'/'+galaxy_name+'_results_'+run_type+'.dat')
    min_mst, max_mst = np.log10(results.table.SigmaMstar).min(), np.log10(results.table.SigmaMstar).max()
    min_sfr, max_sfr = np.log10(results.table.SigmaSFR).min(), np.log10(results.table.SigmaSFR).max()

    ra_apertures = results.table.ra_apertures.values
    dec_apertures = results.table.dec_apertures.values
    RaDec_arr = np.zeros((np.shape(ra_apertures)[0], 2))
    RaDec_arr[:,0], RaDec_arr[:,1] = ra_apertures, dec_apertures
    
    print('Defaul eps value is '+str(eps_value))
    
    check = 'n'
    while not check == 'y':
        # Compute DBSCAN
        X = StandardScaler().fit_transform(RaDec_arr)
        db = DBSCAN(eps=eps_value, min_samples=1).fit(X)
        core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        labels = db.labels_
        # Number of clusters in labels, ignoring noise if present.
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        # Black removed and is used for noise instead.
        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each)
                  for each in np.linspace(0, 1, len(unique_labels))]
        
        fig = plt.figure(figsize=(12,12))
        plt.subplot(2,2,1)
        #plt.scatter(table.ra, table.dec, s = s_size, c = np.log10(table.SFR), cmap = cm.magma, marker = 's')
        for k, col in zip(unique_labels, colors):
            if k == -1:
                # Black used for noise.
                col = [0, 0, 0, 1]
            class_member_mask = (labels == k)
            xy = RaDec_arr[class_member_mask & core_samples_mask]
            plt.scatter(xy[:, 0], xy[:, 1], marker = '.', color = tuple(col), s = s_size/2)
        plt.gca().invert_xaxis() 
        plt.title('Estimated number of clusters: %d' % n_clusters_)
        
        MainCl = np.where(labels == GUS.Most_Common(labels))[0]
        OtherCl = np.where(labels != GUS.Most_Common(labels))[0]
        
        plt.subplot(2,2,2)
        plt.scatter(ra_apertures[MainCl], dec_apertures[MainCl], s = s_size, c = np.log10(results.table.SFR_total[MainCl]), \
                    cmap = cm.magma, marker = 's')
        plt.gca().invert_xaxis()        
        plt.subplot(2,2,3)
        plt.scatter(np.log10(results.table.SigmaMstar), np.log10(results.table.SigmaSFR), s = 10, c = np.log10(results.table.cendist_physical), \
                    cmap = cm.viridis, marker = '.')
        plt.xlim(min_mst, max_mst), plt.ylim(min_sfr, max_sfr)
        plt.subplot(2,2,4)
        plt.scatter(np.log10(results.table.SigmaMstar[MainCl]), np.log10(results.table.SigmaSFR[MainCl]), s = 10, c = np.log10(results.table.cendist_physical[MainCl]), \
                    cmap = cm.viridis, marker = '.')
        plt.xlim(min_mst, max_mst), plt.ylim(min_sfr, max_sfr)
        plt.show()   
        
        answer = input('Do you like what you are seeing? (y/n/break)')
        if answer == 'y': break
        elif answer == 'break': return
        else:
            eps_value = float(input( 'Give me a value for eps (roughly between 0.06, meaning more clusters, and 0.5, meaning less clusters)'))
    
    fig.savefig(magphys_results_folder+'/'+galaxy_name+'_ClusteringRemoval.pdf', bbox_inches = 'tight')
    
    print('Assigning a chisq of 100 to cells not belonging to the main cluster.')
    results.table.Bestfit_Chisq[OtherCl] = 100
    results.save_table(magphys_results_folder+'/'+galaxy_name+'_results_'+run_type+'_ClRm.dat')
    return

def plot_results(GalProp, run_type, aperture_side, fits_base_band, plot_size, s_size, chisq_threshold = 25, radius_threshold = False):
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    rc('text', usetex=True)
    from matplotlib import cm, gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    import matplotlib.axes as maxes

    # =================
    # GALAXY PROPERTIES
    galaxy_name = GalProp.galaxy_name
    D_L = GalProp.dist
    source_z = GalProp.dist*cosmo.H0/v_lux.to('km/s')
    arcsec_per_kpc = cosmo.arcsec_per_kpc_comoving(source_z)
    kpc_per_arcmin = cosmo.kpc_comoving_per_arcmin(source_z)
    ra_cen, dec_cen = GalProp.ra, GalProp.dec
    # =================

    magphys_results_folder = '../'+galaxy_name+'/magphys_'+run_type+'/run_'+run_type
    subprocess.call('mkdir '+magphys_results_folder+'/phys_plots', shell = True)
        
    try: table = pd.read_csv(magphys_results_folder+'/'+galaxy_name+'_results_'+run_type+'_ClRm.dat', \
                        sep = ',', dtype={'Id. Aperture': object})
    except: table = pd.read_csv(magphys_results_folder+'/'+galaxy_name+'_results_'+run_type+'.dat', \
                        sep = ',', dtype={'Id. Aperture': object})

    id_ap = table['Id. Aperture'].values
    ra = table['ra_apertures'].values
    dec = table['dec_apertures'].values
    cendist_physical = table['cendist_physical'].values
    bestfit_chisq = table['Bestfit_Chisq'].values
    MS_distance = table['MS_distance'].values
    Mstar = table['Mstar'].values
    SigmaMstar = table['SigmaMstar'].values
    SFR_magphys = table['SFR_magphys'].values
    sSFR_magphys = table['sSFR_magphys'].values
    SigmaSFR_magphys = table['SigmaSFR_magphys'].values
    SFR_IR = table['SFR_IR'].values
    try: eSFR_IR = table['eSFR_IR'].values
    except: pass
    SFR_UV = table['SFR_UV'].values
    try: eSFR_UV = table['eSFR_UV'].values
    except: pass
    SFR = table['SFR_total'].values
    try: eSFR = table['eSFR_total'].values
    except: pass
    SigmaSFR = table['SigmaSFR'].values
    try: eSigmaSFR = table['eSigmaSFR'].values
    except: pass
    #sSFR = table['sSFR'].values
    A_IRX = table['A_IRX'].values
    Mdust = table['Mdust'].values
    Ldust = table['Ldust'].values
    T_C_ISM = table['T_C_ISM'].values
    T_W_BC = table['T_W_BC'].values
    tauv = table['tauv'].values
    tvism = table['tvism'].values
    fmuIR = table['fmuIR'].values
    fmuSFH = table['fmuSFH'].values
    xi_C_tot = table['xi_C_tot'].values
    xi_MIR_tot = table['xi_MIR_tot'].values
    xi_PAH_tot = table['xi_PAH_tot'].values
    xi_W_tot = table['xi_W_tot'].values
    tform = table['tform'].values
    age_wm = table['age_wm'].values
    
    if radius_threshold:
        print('Assigning a chisq of 100 to every pixel over the radius threshold.')
        try: radius_threshold_unit = radius_threshold.unit
        except: raise Exception('Radius threshold should have a unit!')
        if radius_threshold.unit == 'arcsec': radius_threshold = radius_threshold.to('arcmin')
        if radius_threshold.unit == 'pc': radius_threshold = radius_threshold.to('kpc')
        if radius_threshold.unit == 'arcmin':
            radius_threshold = kpc_per_arcmin*radius_threshold
        over_threshold = np.where(table.cendist_physical > radius_threshold.value)
        table.bestfit_chisq[over_threshold] = 100
    
    GOOD_chisq = np.where(bestfit_chisq <= chisq_threshold)
    BAD_chisq = np.where(bestfit_chisq > chisq_threshold)
    
    MS_plot = plt.figure(figsize=(8,8))
    plt.scatter(np.log10(SigmaMstar[GOOD_chisq]), np.log10(SigmaSFR[GOOD_chisq]), c = cendist_physical[GOOD_chisq], \
             marker = '*', cmap = cm.viridis, linestyle = 'None', label = 'SFR from SED')
    plt.xlabel(r'$\log \Sigma_* [{\rm M}_{\odot} {\rm kpc}^{-2}]$', fontsize = 20)
    plt.ylabel(r'$\log \Sigma_{\rm SFR} [{\rm M}_{\odot} {\rm yr}^{-1} {\rm kpc}^{-2}]$', fontsize = 20)
    plt.legend()
    plt.colorbar()
    MS_plot.savefig(magphys_results_folder+'/phys_plots/logSigmaSFR-logMstar.pdf', bbox_inches = 'tight')
    
    SFR_comparison = plt.figure(figsize=(8,8))
    plt.plot(np.log10(SigmaSFR[GOOD_chisq]), np.log10(SigmaSFR_magphys[GOOD_chisq]), marker = '*', color = 'indianred', linestyle = 'None')
    plt.plot([-8,0], [-8, 0], linestyle = '-', color = 'k', linewidth = 3.0)
    plt.xlim(-5, -1)
    plt.xlabel(r'$\log \Sigma_{\rm SFRKennicutt} [{\rm M}_{\odot} {\rm yr}^{-1} {\rm kpc}^{-2}]$', fontsize = 20)
    plt.ylabel(r'$\log \Sigma_{\rm SFRmagphys} [{\rm M}_{\odot} {\rm yr}^{-1} {\rm kpc}^{-2}]$', fontsize = 20)
    SFR_comparison.savefig(magphys_results_folder+'/phys_plots/SFR_comparison.pdf', bbox_inches = 'tight')
    
    fits_base = GUS.FitsUtils('../'+galaxy_name+'/_ReducedMaps/'+fits_base_band+'.fits')    
    pixel_scale = fits_base.get_pixel_scale()*u.deg.to('arcsec')
    cen_pos = SkyCoord(GalProp.ra, GalProp.dec, frame = 'icrs')
    cen_xpix, cen_ypix = skycoord_to_pixel(cen_pos, fits_base.wcs)
    
    coords = SkyCoord(ra*u.deg, dec*u.deg, frame = 'icrs')
    ra_x, dec_y = coords.to_pixel(fits_base.wcs) # In pixels.

    # ---------
    # Plot size
    try: unit = plot_size.unit
    except: raise Exception('You must give a unit to plot size.')
    if plot_size.unit == 'pc': plot_size = plot_size.to('kpc')
    if plot_size.unit == 'arcmin': plot_size = plot_size.to('arcsec')    
    if plot_size.unit == 'arcsec':
        plot_size_in_kpc = (2*plot_size)/arcsec_per_kpc
        print('Plot side will be of {0:.2f}'.format(plot_size_in_kpc))
    if plot_size.unit == 'kpc':
        plot_size *= arcsec_per_kpc
        print('At z = {0:.6f} plot side will be of {1:.1f}'.format(GalProp.z_source, 2*plot_size))
    plot_size_pixel = (plot_size/pixel_scale).value
    xmin, xmax = cen_xpix-plot_size_pixel, cen_xpix+plot_size_pixel
    ymin, ymax = cen_ypix-plot_size_pixel, cen_ypix+plot_size_pixel
    # ---------

    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(1,1,1, projection = fits_base.wcs)
    ax.coords[0].set_axislabel('RA [deg]', size = 20), ax.coords[1].set_axislabel('DEC [deg]', size = 20)
    ax.coords[0].set_ticklabel(size=15), ax.coords[1].set_ticklabel(size=15)
    c_BAD = np.log10(SFR[BAD_chisq])
    c_GOOD = np.log10(SFR[GOOD_chisq])
    plt.scatter(ra_x[BAD_chisq], dec_y[BAD_chisq], c = c_BAD, s = s_size, marker = '.', alpha = 0.3, cmap = cm.inferno)
    sca = ax.scatter(ra_x[GOOD_chisq], dec_y[GOOD_chisq], c = c_GOOD, s = s_size, marker = 's', alpha = 1.0, cmap = cm.inferno)
    ax.axis('equal')
    if xmin != 0: ax.set_xlim(xmin, xmax), ax.set_ylim(ymin, ymax)
    else: pass
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", axes_class=maxes.Axes, pad=0.0)
    cax.tick_params(direction='in') 
    cbar = fig.colorbar(sca, cax=cax, orientation='horizontal')
    cbar.outline.set_edgecolor('black')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='both', which='major', length = 4.0, labelsize=15)
    cax.set_title(r'$\log {\rm SFR} [M_{\odot}/yr]$', fontsize = 25, pad = 30)
    fig.savefig(magphys_results_folder+'/phys_plots/resolved_SFR.pdf', bbox_inches = 'tight')
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(1,1,1, projection = fits_base.wcs)
    ax.coords[0].set_axislabel('RA [deg]', size = 20), ax.coords[1].set_axislabel('DEC [deg]', size = 20)
    ax.coords[0].set_ticklabel(size=15), ax.coords[1].set_ticklabel(size=15)
    c_BAD = np.log10(Mstar[BAD_chisq])
    c_GOOD = np.log10(Mstar[GOOD_chisq])
    plt.scatter(ra_x[BAD_chisq], dec_y[BAD_chisq], c = c_BAD, s = s_size, marker = '.', alpha = 0.3, cmap = cm.magma)
    sca = ax.scatter(ra_x[GOOD_chisq], dec_y[GOOD_chisq], c = c_GOOD, s = s_size, marker = 's', alpha = 1.0, cmap = cm.magma)
    ax.axis('equal')
    if xmin != 0: ax.set_xlim(xmin, xmax), ax.set_ylim(ymin, ymax)
    else: pass
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", axes_class=maxes.Axes, pad=0.0)
    cax.tick_params(direction='in') 
    cbar = fig.colorbar(sca, cax=cax, orientation='horizontal')
    cbar.outline.set_edgecolor('black')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='both', which='major', length = 4.0, labelsize=15)
    cax.set_title(r'$\log {\rm M}_* [M_{\odot}]$', fontsize = 25, pad = 30)
    fig.savefig(magphys_results_folder+'/phys_plots/resolved_Mstar.pdf', bbox_inches = 'tight')
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(1,1,1, projection = fits_base.wcs)
    ax.coords[0].set_axislabel('RA [deg]', size = 20), ax.coords[1].set_axislabel('DEC [deg]', size = 20)
    ax.coords[0].set_ticklabel(size=15), ax.coords[1].set_ticklabel(size=15)
    plt.scatter(ra_x[BAD_chisq], dec_y[BAD_chisq], c = MS_distance[BAD_chisq], s = s_size, marker = '.', alpha = 0.3, cmap = cm.RdBu)
    plt.clim(-1,1)
    sca = ax.scatter(ra_x[GOOD_chisq], dec_y[GOOD_chisq], c = MS_distance[GOOD_chisq], s = s_size, marker = 's', alpha = 1.0, cmap = cm.RdBu)
    sca.set_clim(-1,1)
    ax.axis('equal')
    if xmin != 0: ax.set_xlim(xmin, xmax), ax.set_ylim(ymin, ymax)
    else: pass
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", axes_class=maxes.Axes, pad=0.0)
    cax.tick_params(direction='in') 
    cbar = fig.colorbar(sca, cax=cax, orientation='horizontal')
    cbar.outline.set_edgecolor('black')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='both', which='major', length = 4.0, labelsize=15)
    cax.set_title(r'${\rm Distance}$ ${\rm from}$ ${\rm MS}$', fontsize = 25, pad = 30)
    fig.savefig(magphys_results_folder+'/phys_plots/resolved_MSdistance.pdf', bbox_inches = 'tight')
    
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(1,1,1, projection = fits_base.wcs)
    ax.coords[0].set_axislabel('RA [deg]', size = 20), ax.coords[1].set_axislabel('DEC [deg]', size = 20)
    ax.coords[0].set_ticklabel(size=15), ax.coords[1].set_ticklabel(size=15)
    plt.scatter(ra_x[BAD_chisq], dec_y[BAD_chisq], c = bestfit_chisq[BAD_chisq], s = s_size, marker = '.', alpha = 0.3, cmap = cm.RdBu_r)
    plt.clim(0,2*chisq_threshold)
    sca = ax.scatter(ra_x[GOOD_chisq], dec_y[GOOD_chisq], c = bestfit_chisq[GOOD_chisq], s = s_size, marker = 's', alpha = 1.0, cmap = cm.RdBu_r)
    sca.set_clim(0,2*chisq_threshold)
    ax.axis('equal')
    if xmin != 0: ax.set_xlim(xmin, xmax), ax.set_ylim(ymin, ymax)
    else: pass
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", axes_class=maxes.Axes, pad=0.0)
    cax.tick_params(direction='in') 
    cbar = fig.colorbar(sca, cax=cax, orientation='horizontal')
    cbar.outline.set_edgecolor('black')
    cax.xaxis.set_ticks_position('top')
    cax.tick_params(axis='both', which='major', length = 4.0, labelsize=15)
    cax.set_title(r'$\chi^2$', fontsize = 25, pad = 30)
    fig.savefig(magphys_results_folder+'/phys_plots/resolved_chisquare.pdf', bbox_inches = 'tight')
    
    return

def plot_seds(GalProp, run_type, aperture_side, fits_base_band, plot_size, \
              vmin = False, vmax = False, ymin = 1E-7, ymax = 1E2, tile_color = 'w', fill_tile = False, \
              map_title = False, overwrite = False, make_m4a = True, sort_by_cendist = True):
    
    from astropy.cosmology import Planck15 as cosmo
    from astropy.constants import L_sun
    import matplotlib.pyplot as plt
    from matplotlib import cm, gridspec
    from astropy.visualization import (MinMaxInterval, LogStretch, ImageNormalize)
    from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator, LogLocator)
    
    # =================
    # GALAXY PROPERTIES
    galaxy_name = GalProp.galaxy_name
    D_L = GalProp.dist
    z_source = GalProp.z_source
    DL = cosmo.luminosity_distance(z_source).to('m') # Sta roba è da risolvere.
    arcsec_per_kpc = cosmo.arcsec_per_kpc_comoving(z_source)
    ra_cen, dec_cen = GalProp.ra, GalProp.dec
    cen_pos = SkyCoord(ra_cen, dec_cen, frame = 'icrs')
    # =================

    # =================
    magphys_results_folder = '../'+galaxy_name+'/magphys_'+run_type+'/run_'+run_type
    # =================

    # =====================
    # Paths for .sed and .fit
    sed_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.sed')])
    fit_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.fit')])
    # =====================
    
    # ==============
    # Read some physical properties
    table = pd.read_csv(magphys_results_folder+'/'+galaxy_name+'_results_'+run_type+'.dat', \
                        sep = ',', dtype={'Id. Aperture': object})
    id_ap = table['Id. Aperture'].values
    ra_apertures = table['ra_apertures'].values*u.deg
    dec_apertures = table['dec_apertures'].values*u.deg
    cendist_physical = table['cendist_physical'].values
    bestfit_chisq = table['Bestfit_Chisq'].values
    
    if aperture_side.unit == 'pc': aperture_side = aperture_side.to('kpc')
    if aperture_side.unit == 'arcmin': aperture_side = aperture_side.to('arcsec')
    if aperture_side.unit == 'kpc': aperture_side *= arcsec_per_kpc
    else: pass
    coords = SkyCoord(ra_apertures, dec_apertures, frame = 'icrs')
    list_apertures = []
    for i_ap in trange(len(ra_apertures), ascii = True, desc = 'Reading apertures'):
        aperture = SkyRectangularAperture(coords[i_ap], w=aperture_side, h=aperture_side, theta = 0*u.deg)
        list_apertures.append(aperture)
    # ==============

    # ==============
    # Now start with the plot
    fits_base = GUS.FitsUtils('../'+galaxy_name+'/_ReducedMaps/'+fits_base_band+'.fits') 
    pixel_scale = (fits_base.get_pixel_scale()*u.deg).to('arcsec')
    cen_pos = SkyCoord(GalProp.ra, GalProp.dec, frame = 'icrs')
    cen_xpix, cen_ypix = skycoord_to_pixel(cen_pos, fits_base.wcs)
    try: unit = plot_size.unit
    except: raise Exception('You must give a unit to plot size.')
    if unit == 'pc': plot_size = plot_size.to('kpc')
    if unit == 'arcmin': plot_size = plot_size.to('arcsec')    
    if unit == 'arcsec':
        plot_size_in_kpc = 2*plot_size/arcsec_per_kpc
        print('Plot side will be of {0:.2f}'.format(plot_size_in_kpc))
    if unit == 'kpc':
        plot_size *= arcsec_per_kpc
        print('At z = {0:.6f} plot side will be of {1:.1f}'.format(GalProp.z_source, 2*plot_size))
    plot_size_pixel = (plot_size/pixel_scale).value
    x_min, x_max = cen_xpix-plot_size_pixel, cen_xpix+plot_size_pixel
    y_min, y_max = cen_ypix-plot_size_pixel, cen_ypix+plot_size_pixel
    # ==============

    # ==============
    subprocess.call('mkdir '+magphys_results_folder+'/seds', shell = True)
    for fit, sed, aperture, chi_sq in tqdm(zip(fit_path, sed_path, list_apertures, bestfit_chisq)):
        
        plt.close('all') # Mandatory for survival
        i_file = int(sed[-8:-4])
        
        if overwrite: pass
        else:
            if os.path.exists(magphys_results_folder+'/seds/{0:04}_sed.png'.format(i_file)): continue
    
        SED = plt.figure(figsize=(15,5))
    
        gs = gridspec.GridSpec(1, 2, width_ratios=[1.8, 1]) 
        ax1 = plt.subplot(gs[0]) # SED with observed datapoints + magphys fit
        ax2 = plt.subplot(gs[1], projection = fits_base.wcs) # Galaxy with selected aperture
    
        temp_res = magphys_read.MagphysOutput(fit, sed)
        filters = temp_res.obs_filters
        lam_observed = [Wvlghts_dictionary[i] for i in filters]*u.micrometer
    
        SED_observed = 1E26*(temp_res.obs_flux.value*L_sun)/(4*np.pi*DL**2)             # From LSOL HZ^-1 To Jansky
        SED_observed_errors = 1E26*(temp_res.obs_flux_err.value*L_sun)/(4*np.pi*DL**2)  # From LSOL HZ^-1 To Jansky
        SED_predicted = 1E26*(temp_res.obs_predict.value*L_sun)/(4*np.pi*DL**2)         # From LSOL HZ^-1 To Jansky
        SED_snr= SED_observed/SED_observed_errors
        
        lambda_model = (temp_res.sed_model_logwaves).to('Angstrom')
        y = 10**temp_res.sed_model[:,1]*(lambda_model**2/v_lux.to('Angstrom/s'))
        fsed = (y*(1+z_source)/(4*np.pi*DL**2))*1E26 # Ricordati di cambiare a 1E23 quando arriva blabla
        fsed = fsed*L_sun
        
        cond_colors = ['tab:red' if t < 2 else \
                      ('orange' if t > 2 and t < 3 else \
                      ('gold' if t > 3 and t < 5 else \
                      ('lightgreen' if t > 5 and t < 8 else \
                      ('tab:green' if t > 8 and t < 10 else 'darkgreen')))) for t in SED_snr]
                  
        ax1.plot(lambda_model.to('micrometer'), fsed, color = 'k', linewidth = 3.0, zorder = 0)
        ax1.scatter(lam_observed, SED_observed, s = 80, color = cond_colors, linestyle = 'None', marker = 'o', zorder = 1)
        
        SEDpoint_labels = ['SNR < 2', '2 < SNR < 3', '3 < SNR < 5', '5 < SNR < 8', '8 < SNR < 10', 'SNR > 10']
        SEDpoint_colors = ['tab:red', 'orange', 'gold', 'lightgreen', 'tab:green', 'darkgreen']
        [ax1.plot(0, 0, color = c, linewidth = 0, ms = 10, marker = 'o', label = l) for c, l in zip(SEDpoint_colors, SEDpoint_labels)]
        ax1.legend(loc = 'upper left', frameon=False)
    
        cond_color = 'tab:red' if chi_sq > 25 else ('gold' if chi_sq > 15 and chi_sq < 25 else 'tab:green')
        ax1.text(1E2, 2E1, r'$\chi^2$  = '+str(np.round(chi_sq, 1)), color = cond_color, fontsize = 15)
        ax1.set_xscale('log'), ax1.set_yscale('log')
        ax1.set_xlabel(r'$\log_{10} (\lambda/\mu{\rm m})$', size = 15), ax1.set_ylabel(r'$\log_{10} (f_{\lambda}\,/\,{\rm Jy})$', size = 15)
        ax1.tick_params(axis='both', which='major', direction = 'in', labelsize=15)
        ax1.tick_params(axis='both', which='minor', direction = 'in', labelsize=8)
        ax1.yaxis.set_minor_locator(LogLocator(base=10))
    
        #ymin, ymax = 10**(np.log10(SED_observed.min().value)-0.5), 10**(np.log10(SED_observed.max().value)+0.5)
        ax1.set_xlim(1E-1, 7E2), ax1.set_ylim(ymin, ymax)
    
        # Create an ImageNormalize object
        if vmin: norm = ImageNormalize(fits_base.signal, interval=MinMaxInterval(), stretch=LogStretch(), vmin = vmin, vmax = vmax)
        else: norm = ImageNormalize(fits_base.signal, interval=MinMaxInterval(), stretch=LogStretch())
        ax2.imshow(fits_base.signal, origin = 'lower', interpolation = 'nearest', cmap = GUS.associate_colormap(fits_base.bandname), norm = norm)
        if map_title: ax2.set_title(fits_base_band)
    
        ra_x, dec_y = aperture.positions.to_pixel(fits_base.wcs)[0], aperture.positions.to_pixel(fits_base.wcs)[1]
        term = (.5*aperture_side/pixel_scale).value
        rectangle = plt.Rectangle((ra_x-term, dec_y+term),
                                  width = aperture_side/pixel_scale, height = aperture_side/pixel_scale,
                                  color = tile_color, linewidth = 3.0, fill = fill_tile)
        ax2.add_artist(rectangle)
    
        ax2.set_xlim(x_min, x_max), ax2.set_ylim(y_min, y_max)
        ax2.coords[1].set_axislabel_position('r'), ax2.coords[1].set_ticks_position('r'), ax2.coords[1].set_ticklabel_position('r')
        ax2.coords[0].set_axislabel('RA [deg]', size = 15), ax2.coords[1].set_axislabel('DEC [deg]', size = 15)
        ax2.coords[0].set_ticklabel(size=10), ax2.coords[1].set_ticklabel(size=10)
        
        plt.tight_layout()
        SED.savefig(magphys_results_folder+'/seds/{0:04}_sed.png'.format(i_file), bbox_inches = 'tight')
        
    if make_m4a:
        import imageio
        import moviepy.editor as mp
        # Generate .gif
        path = magphys_results_folder+'/seds/'
        fullpath = [path+filename for filename in sorted(os.listdir(path)) if not filename.startswith('.')]
        if sort_by_cendist:
            print('.gif and .mp4 will be sorted by distance from the center.')
            fullpath = [x for _, x in sorted(zip(cendist_physical, fullpath), key=lambda pair: pair[0])]
        else: print('.gif and .mp4 will be sorted by aperture number (top-left to bottom-right).')
        images = [imageio.imread(t) for t in fullpath]
        imageio.mimsave(magphys_results_folder+'/seds.gif', images, duration = 0.1)
        # Generate .mp4
        clip = mp.VideoFileClip(magphys_results_folder+'/seds.gif')
        clip.write_videofile(magphys_results_folder+'/seds.mp4')
        
    return

def plot_seds_old(GalProp, run_type, aperture_side, fits_base_band, plot_size):
    from astropy.constants import L_sun
    import matplotlib.pyplot as plt
    from matplotlib import cm, gridspec
    from astropy.visualization import (MinMaxInterval, LogStretch, ImageNormalize)

    # =================
    # GALAXY PROPERTIES
    galaxy_name = GalProp.galaxy_name
    D_L = GalProp.dist
    z_source = GalProp.z_source
    DL = cosmo.luminosity_distance(z_source).to('m') # Sta roba è da risolvere.
    arcsec_per_kpc = cosmo.arcsec_per_kpc_comoving(z_source)
    ra_cen, dec_cen = GalProp.ra, GalProp.dec
    cen_pos = SkyCoord(ra_cen, dec_cen, frame = 'icrs')
    # =================
    
    magphys_results_folder = '../'+galaxy_name+'/magphys_'+run_type+'/run_'+run_type

    # =====================
    # Paths for .sed and .fit
    sed_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.sed')])
    fit_path = sorted([magphys_results_folder+'/final_results/'+file \
                   for file in os.listdir(magphys_results_folder+'/final_results/') if file.endswith('.fit')])
    # =====================

    # ==============
    # Read some physical properties
    table = pd.read_csv(magphys_results_folder+'/'+galaxy_name+'_results_'+run_type+'.dat', \
                        sep = ',', dtype={'Id. Aperture': object})
    id_ap = table['Id. Aperture'].values
    ra_apertures = table['ra_apertures'].values*u.deg
    dec_apertures = table['dec_apertures'].values*u.deg
    cendist_physical = table['cendist_physical'].values
    bestfit_chisq = table['Bestfit_Chisq'].values
    SFRs = table['SFR_total'].values
    Mstars = table['Mstar'].values

    if aperture_side.unit == 'pc': aperture_side = aperture_side.to('kpc')
    if aperture_side.unit == 'arcmin': aperture_side = aperture_side.to('arcsec')
    if aperture_side.unit == 'kpc': aperture_side *= arcsec_per_kpc
    else: pass
    coords = SkyCoord(ra_apertures, dec_apertures, frame = 'icrs')
    list_apertures = []
    for i_ap in trange(len(ra_apertures), ascii = True, desc = 'Reading apertures'):
        aperture = SkyRectangularAperture(coords[i_ap], w=aperture_side, h=aperture_side, theta = 0*u.deg)
        list_apertures.append(aperture)
    # ==============

    # ==============
    # Sort w.r.t cendist
    _, apertures_sorted, sed_path_sorted, fit_path_sorted, bestfit_chisq_sorted, SFRs_sorted, Mstars_sorted = \
         (t for t in zip(*sorted(zip(cendist_physical, list_apertures, sed_path, fit_path, bestfit_chisq, SFRs, Mstars))))
    
    # Colorscale
    cm_subsection = np.linspace(0, 1, len(apertures_sorted)) 
    colors = [cm.inferno(x) for x in cm_subsection]
    # ==============

    # Now start with the plot
    fits_base = GUS.FitsUtils('../'+galaxy_name+'/_ReducedMaps/'+fits_base_band+'.fits') 
    pixel_scale = (fits_base.get_pixel_scale()*u.deg).to('arcsec')
    cen_pos = SkyCoord(GalProp.ra, GalProp.dec, frame = 'icrs')
    cen_xpix, cen_ypix = skycoord_to_pixel(cen_pos, fits_base.wcs)
    try: unit = plot_size.unit
    except: raise Exception('You must give a unit to plot size.')
    if unit == 'pc': plot_size = plot_size.to('kpc')
    if unit == 'arcmin': plot_size = plot_size.to('arcsec')    
    if unit == 'arcsec':
        plot_size_in_kpc = 2*plot_size/arcsec_per_kpc
        print('Plot side will be of {0:.2f}'.format(plot_size_in_kpc))
    if unit == 'kpc':
        plot_size *= arcsec_per_kpc
        print('At z = {0:.6f} plot side will be of {1:.1f}'.format(GalProp.z_source, 2*plot_size))
    plot_size_pixel = (plot_size/pixel_scale).value
    x_min, x_max = cen_xpix-plot_size_pixel, cen_xpix+plot_size_pixel
    y_min, y_max = cen_ypix-plot_size_pixel, cen_ypix+plot_size_pixel

    subprocess.call('mkdir '+magphys_results_folder+'/seds', shell = True)
    for fit, sed, aperture, chi_sq, SFR, Mstar, i_ap \
    in tqdm(zip(fit_path_sorted, sed_path_sorted, apertures_sorted, bestfit_chisq_sorted, SFRs_sorted, Mstars_sorted, range(len(apertures_sorted)))):
        
        plt.close('all') # Mandatory for survival
        i_file = int(sed[-8:-4])
        
        #if os.path.exists(folder+'/seds/{0:04}_sed.png'.format(i_file)): continue
    
        SED = plt.figure(figsize=(15,5))
    
        gs = gridspec.GridSpec(1, 2, width_ratios=[1.8, 1]) 
        ax1 = plt.subplot(gs[0]) # SED with observed datapoints + magphys fit
        ax2 = plt.subplot(gs[1], projection = fits_base.wcs) # Galaxy with selected aperture
    
        temp_res = magphys_read.MagphysOutput(fit, sed)
        filters = temp_res.obs_filters
        lam_observed = [Wvlghts_dictionary[i] for i in filters]*u.micrometer
    
        SED_observed = 1E26*(temp_res.obs_flux.value*L_sun)/(4*np.pi*DL**2)             # From LSOL HZ^-1 To Jansky
        SED_observed_errors = 1E26*(temp_res.obs_flux_err.value*L_sun)/(4*np.pi*DL**2)  # From LSOL HZ^-1 To Jansky
        SED_predicted = 1E26*(temp_res.obs_predict.value*L_sun)/(4*np.pi*DL**2)         # From LSOL HZ^-1 To Jansky
            
        lambda_model = (temp_res.sed_model_logwaves).to('Angstrom')
        y = 10**temp_res.sed_model[:,1]*(lambda_model**2/v_lux.to('Angstrom/s'))
        fsed = (y*(1+z_source)/(4*np.pi*DL**2))*1E26 # Ricordati di cambiare a 1E23 quando arriva blabla
        fsed = fsed*L_sun
    
        ax1.plot(lambda_model.to('micrometer'), fsed, color = colors[i_ap], linewidth = 3.0, label = 'Magphys fit')
        ax1.plot(lam_observed, SED_observed, linestyle = 'None', marker = 's', \
                 markersize = 12, color = colors[i_ap], label = 'Observed')
        ax1.text(1E2, 1E1, r'$\chi^2$  = '+str(np.round(chi_sq, 1)), fontsize = 13)
        ax1.text(1E1, 2E1, r'$\log$ M$_*$  = {0:.2f}'.format(np.log10(Mstar)), fontsize = 13)
        ax1.text(1E1, 8E0, r'$\log$ SFR  = {0:.2f}'.format(np.log10(SFR)), fontsize = 13)
        ax1.set_xscale('log'), ax1.set_yscale('log')
        ax1.set_xlabel(r'Wavelength ($\mu$m)'), ax1.set_ylabel(r'Flux (Jy)')
        #ymin, ymax = 10**(np.log10(SED_observed.min().value)-0.5), 10**(np.log10(SED_observed.max().value)+0.5)
        ymin, ymax = 1E-7, 1E2 # Fixed y axis
        ax1.set_xlim(1E-1, 7E2), ax1.set_ylim(ymin, ymax)
        #ax1.legend(loc='upper left')

        # Create an ImageNormalize object
        norm = ImageNormalize(fits_base.signal, interval=MinMaxInterval(), stretch=LogStretch())
        ax2.imshow(fits_base.signal, origin = 'lower', interpolation = 'nearest', cmap = cm.viridis, norm = norm)
        #ax2.set_title(fits_base_band)
    
        ra_x, dec_y = aperture.positions.to_pixel(fits_base.wcs)[0], aperture.positions.to_pixel(fits_base.wcs)[1]
        term = (.5*aperture_side/pixel_scale).value
        rectangle = plt.Rectangle((ra_x-term, dec_y+term),
                                  width = aperture_side/pixel_scale, height = aperture_side/pixel_scale,
                                  color='w', linewidth = 3.0, fill = False)
        ax2.add_artist(rectangle)
    
        ax2.set_xlim(x_min, x_max), ax2.set_ylim(y_min, y_max)
    
        plt.tight_layout()
        SED.savefig(magphys_results_folder+'/seds/{0:04}_sed.png'.format(i_file), bbox_inches = 'tight')

    return

def map_gridding(phys_prop, ra_ap, dec_ap, pixel_scale):
    pixel_scale = (8*u.arcsec).to('deg')
    xxx = np.array(ra_ap*u.deg/pixel_scale)
    yyy = np.array(dec_ap*u.deg/pixel_scale)
    xxx -= xxx.min()
    yyy -= yyy.min()
    
    int_xxx = [int(np.round(t,1)) for t in xxx]
    int_yyy = [int(np.round(t,1)) for t in yyy]
    
    x_, x_idx = np.unique(np.ravel(int_yyy), return_inverse=True)
    y_, y_idx = np.unique(np.ravel(int_xxx), return_inverse=True)
    newArray = np.zeros((len(x_), len(y_)), dtype=phys_prop.dtype)
    newArray[x_idx, y_idx] = np.ravel(phys_prop)
    newArray[np.where(newArray == 0)] = np.nan
    return newArray, x_, y_

def save_fits(newArray, fits_path, pixel_scale, ref_ra, ref_dec, nx, ny, ref_pixel_nx, ref_pixel_ny):
    hdu = fits.PrimaryHDU(newArray)
    header = hdu.header
    header['NAXIS1'] = nx
    header['NAXIS2'] = ny
    header['CRVAL1'] = ref_ra
    header['CRVAL2'] = ref_dec
    header['CRPIX1'] = ref_pixel_nx
    header['CRPIX2'] = ref_pixel_ny
    header['CDELT1'] = pixel_scale.value
    header['CDELT2'] = pixel_scale.value
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CTYPE1'] = 'RA---TAN'
    header['CTYPE2'] = 'DEC--TAN'
    hdu.header = header
    hdu.writeto(fits_path+'.fits', overwrite = True)

def generate_fits(GalProp, run_type, pixel_scale, phys_prop):
    
    # GALAXY PROPERTIES
    galaxy_name = GalProp.galaxy_name
    physical_scale = cosmo.arcsec_per_kpc_comoving(z_source)
    ra_cen, dec_cen = GalProp.ra, GalProp.dec

    table = pd.read_csv('run_pixBYpix_newproc/NGC0628_results_pixBYpix.dat')
    ra_ap, dec_ap = table['ra_apertures'].values, table['dec_apertures'].values
    
    if pixel_scale.unit == 'pc': pixel_scale.to('kpc')
    elif pixel_scale.unit == 'arcmin': pixel_scale.to('arcsec')
    if unit == 'arcsec':
        plot_size_in_kpc = 2*plot_size/physical_scale
        print('Plot side will be of {0:.2f}'.format(plot_size_in_kpc))
    if unit == 'kpc':
        plot_size *= physical_scale
        print('At z = {0:.6f} plot side will be of {1:.1f}'.format(GalProp.z_source, 2*plot_size))

    fits_path = '../'+galaxy_name+'/magphys_'+run_type+'/run_'+run_type
    newArray, x_, y_ = generate_fits(table[phys_prop].values, ra_ap, dec_ap, pixel_scale)
    save_fits(newArray, fits_path, pixel_scale, ra_cen, dec_cen, len(x_), len(y_), len(x_)/2, len(y_)/2)
    return