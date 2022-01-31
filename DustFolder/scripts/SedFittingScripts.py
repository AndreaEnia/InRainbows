#!/usr/bin/python
# -*- coding: latin-1 -*-

import os, subprocess
import multiprocessing
import pandas as pd
import numpy as np
import math
import time 
import GenericUsefulScripts as GUS

def create_filter_file(working_bands, magphys_folder, photometry_folder, excluded_band = False):
    magphys_filter_dictionary = pd.read_csv('../scripts/magphys_filters.csv', sep = '\t')
    disam_bands = [GUS.band_disambiguation(t) for t in working_bands]
    # Create filter, save it to magphys_folder
    filters = magphys_filter_dictionary[magphys_filter_dictionary['name'].isin(disam_bands)]
    filters['fit?'] = [1 for t in working_bands]
    if excluded_band:
        disam_excluded_bands = [GUS.band_disambiguation(t) for t in excluded_band]
        filters.loc[filters['name'].isin(disam_excluded_bands), 'fit?'] = 0
    filters.to_csv(magphys_folder+'/filters.dat', index = False, sep = '\t')
    # And also in Photometry/observations
    filters.to_csv(photometry_folder+'/observations/filters.dat', index = False, sep = '\t')
    return

def prepare_magphys(galaxy_name, run_type, working_bands, entries_per_file, exclude_filters = False):
    
    magphys_folder = '~/magphys/'+galaxy_name+'/run_'+run_type+'/'
    photometry_folder = '../'+galaxy_name+'/'+'Photometry_'+run_type+'/'
    subprocess.call('mkdir '+photometry_folder+'observations', shell = True)
    subprocess.call('mkdir ~/magphys/'+galaxy_name, shell = True)
    subprocess.call('mkdir ~/magphys/'+galaxy_name+'/run_'+run_type, shell = True)
    
    create_filter_file(working_bands, magphys_folder, photometry_folder, excluded_band = exclude_filters)
    
    # Remove previous obs. files, to avoid confusion
    subprocess.call(['rm '+photometry_folder+'observations/obs*'], shell = True)
    
    # Split the big table into multiple files.
    Magphys_table = pd.read_csv(photometry_folder+'Final_for_MAGPHYS.csv', sep = '\t', dtype={'id': object})
    number_of_apertures = len(Magphys_table.index)
    subset_zero = Magphys_table.loc[0:0]
    subset_zero.to_csv(photometry_folder+'observations/obs_zero.dat', index = False, sep = '\t')
    z_fill_number = int(np.trunc(np.log10(number_of_apertures)))+2
    aperture_idxs = [str(idx).zfill(z_fill_number) for idx in range(number_of_apertures)]
    #aperture_idxs = [idx.zfill(z_fill_number) for idx in number_of_apertures.astype(str)]
    magphys_subtable = Magphys_table[Magphys_table['id'].isin(aperture_idxs)]
        
    if type(entries_per_file) == float or type(entries_per_file) == int:
        entries_per_file = int(entries_per_file)
        occupied_cores = int(math.ceil(number_of_apertures/entries_per_file))
        if occupied_cores > multiprocessing.cpu_count():
            raise ValueError('WARNING, you want to occupy more cores ({0:2d}) than available on this machine ({1:2d})'.format(occupied_cores, multiprocessing.cpu_count()))
        print('Each file will have {0:3d} entries. This will occupy {1:2d} cores'.format(entries_per_file, occupied_cores))
        est_time_min = time.strftime("%H:%M", time.gmtime(60*3.0*entries_per_file))
        est_time_max = time.strftime("%H:%M", time.gmtime(60*7.0*entries_per_file))
        print('Est. time between '+est_time_min+' hours and '+est_time_max+' hours')
        print('(based on a badly rough esteem between 3 and 7 minutes per aperture)')

    if entries_per_file == 'fill_cores':
        occupied_cores = multiprocessing.cpu_count()
        entries_per_file = int(math.ceil(number_of_apertures/occupied_cores))
        print('You want to fill the cores ({0:2d}) so each file will have {1:3d} entries.'.format(occupied_cores, entries_per_file))
        est_time_min = time.strftime("%H:%M", time.gmtime(60*3.0*entries_per_file))
        est_time_max = time.strftime("%H:%M", time.gmtime(60*7.0*entries_per_file))
        print('Est. time between '+est_time_min+' hours and '+est_time_max+' hours')
        print('(based on a rough esteem between 3 and 7 minutes per entry)')
    elif entries_per_file == 'almost_fill_cores':
        occupied_cores = multiprocessing.cpu_count()-3
        entries_per_file = int(math.ceil(number_of_apertures/occupied_cores))
        print('You almost want to fill the cores ({0:2d}) so each file will have {1:3d} entries.'.format(occupied_cores, entries_per_file))
        est_time_min = time.strftime("%H:%M", time.gmtime(60*3.0*entries_per_file))
        est_time_max = time.strftime("%H:%M", time.gmtime(60*7.0*entries_per_file))
        print('Est. time between '+est_time_min+' hours and '+est_time_max+' hours')
        print('(based on a rough esteem between 3 and 7 minutes per entry)')
    elif entries_per_file == 'single':
        occupied_cores = 1
        print('You choose to run everything on a single non-parallel run.')
        print('You weirdo.')
        est_time = time.strftime("%H:%M", time.gmtime(60*3.0*number_of_apertures))
        print('Est. time between '+est_time_min+' hours and '+est_time_max+' hours')
        print('(based on a badly rough esteem between 3 and 7 minutes per aperture)')
        print('I beg you to reconsider')
        entries_per_file = np.copy(number_of_apertures)
    elif entries_per_file > number_of_apertures:
        raise Exception('Entries per file should not exceed the number of data points, '+str(number_of_apertures))
    else: pass
    
    list_of_subsets = np.array_split(magphys_subtable, occupied_cores)
    for subset in list_of_subsets:
        start_val, end_val = str(subset['id'].values[0]), str(subset['id'].values[-1])
        subset.to_csv(photometry_folder+'observations/observations_'+start_val+'_'+end_val+'.dat', index = False, sep = '\t')
 
    # Purtroppo va fatto a mano.
    z_source = np.round(Magphys_table['z'][0], 4)
    make_zlibs(str(z_source))
        
    # Copy and paste the splitted observations into the magphys folder.
    subprocess.call(['rm '+magphys_folder+'obs*'], shell = True)
    subprocess.call(['cp -rp '+photometry_folder+'observations/obs* '+magphys_folder+'.'], shell = True)
    return

def make_zlibs(z_source):
    f = open('/home/dustpedia/magphys/zlibs.dat', 'w')
    f.write('     1  '+z_source)
    f.close()
    return
    
def run_magphys(galaxy_name, run_type):
    # Run the Jewels
    os.chdir('/home/dustpedia/magphys/')
    subprocess.call('sh run_SED_fitting.sh '+galaxy_name+' '+run_type, shell = True)
    print('MAGPHYS is running. Come back later for the results. Bye.')
    return

def run_magphys_leftovers(galaxy_name, run_type):
    # Run the Jewels
    os.chdir('/home/dustpedia/magphys/')
    subprocess.call('sh run_SED_fitting_leftovers.sh '+galaxy_name+' '+run_type, shell = True)
    print('MAGPHYS leftovers are running. Come back later for the results. Bye.')
    return

def move_results(galaxy_name, run_type, original_path):
    try:
        zero_size_files = subprocess.check_output('find ~/magphys/*.fit -type f -size 0c', shell = True)
        zero_size_files = [t.split('/')[-1].split('.')[0] for t in zero_size_files.decode().split('\n')[:-1]]
        print('MAGPHYS is now processing these apertures:', zero_size_files)
    except:
        zero_size_files = 'No zero size files'
        return 'over'
    
    if len(subprocess.check_output('find ~/magphys/*.fit -type f -size 0c', shell = True)) == 0:
        os.chdir(original_path)
        magphys_results_folder = original_path+'/magphys_'+run_type+'/run_'+run_type+'/final_results/'
        magphys_folder = '~/magphys/'+galaxy_name+'/run_'+run_type
        print('Results will be moved in folders:')
        print(magphys_results_folder)
        print(magphys_folder)
        answer = input('Are you ok with that? (y/n)')
        if answer == 'y': pass
        else: magphys_results_folder = input('Give me an alternative path:')         
        subprocess.call('mkdir magphys_'+run_type, shell = True)
        subprocess.call('mkdir magphys_'+run_type+'/run_'+run_type, shell = True)
        subprocess.call('mkdir '+magphys_results_folder, shell = True)
        subprocess.call('cp ~/magphys/*.sed '+magphys_results_folder, shell = True)
        subprocess.call('cp ~/magphys/*.fit '+magphys_results_folder, shell = True)
        subprocess.call('mv ~/magphys/*.sed '+magphys_folder, shell = True)
        subprocess.call('mv ~/magphys/*.fit '+magphys_folder, shell = True)
        return 'over'
    else:
        print('MAGPHYS is still running. Come back later for the results. Bye.')
        return 'running'
    
def check_results(galaxy_name, run_type, galaxy_folder_path):
    from os import listdir
    from os.path import isfile, join
    mypath = galaxy_folder_path+'/magphys_'+run_type+'/run_'+run_type+'/final_results/'
    onlyfiles = sorted([f.split('.')[0] for f in listdir(mypath) if isfile(join(mypath, f))])
    onlyfiles = list(dict.fromkeys(onlyfiles))
    good_cells = [int(t) for t in onlyfiles]
    number_of_cells = len(pd.read_csv(galaxy_folder_path+'/Photometry_'+run_type+'/Final_for_MAGPHYS.csv', sep = '\t'))
    total_cells = list(np.arange(number_of_cells))
    total_cells = [int(t) for t in total_cells]
    remaining_cells = np.array(list(set(good_cells).symmetric_difference(set(total_cells))))
    if len(remaining_cells) == 0: return print('All cells have been runned.')
    else: print('These cells have been skipped by MAGPHYS', remaining_cells)
    answer = input('Do you want me to run those cells? (y/n)')
    if answer == 'n': return print('Too bad.')
    elif answer == 'y': rerun_magphys(galaxy_name, run_type, galaxy_folder_path, remaining_cells)
    else: return('Give me a proper answer, you asshole')
    return
    
def rerun_magphys(galaxy_name, run_type, galaxy_folder_path, remaining_cells):
    import multiprocessing

    magphys_table = pd.read_csv(galaxy_folder_path+'/Photometry_'+run_type+'/Final_for_MAGPHYS.csv', sep = '\t', dtype={'id': object})
    z_fill_number = int(np.trunc(np.log10(remaining_cells.max())))+2
    remaining_idxs = [idx.zfill(z_fill_number) for idx in remaining_cells.astype(str)]
    Magphys_subtable = magphys_table[magphys_table['id'].isin(remaining_idxs)]
    subprocess.call('mkdir ~/magphys/'+galaxy_name+'/run_'+run_type+'/leftovers', shell = True)
    magphys_folder = '~/magphys/'+galaxy_name+'/run_'+run_type+'/leftovers/'
     
    # Remove previous leftovers files, to avoid confusion
    photometry_folder = '../'+galaxy_name+'/Photometry_'+run_type+'/'
    subprocess.call(['rm '+photometry_folder+'observations/left*'], shell = True)
    subprocess.call(['rm '+magphys_folder+'/left*'], shell = True)
    
    # Split the remaining apertures in maximum 15 cores   
    list_of_subsets = np.array_split(Magphys_subtable, 15)
    for idx, subset in enumerate(list_of_subsets):
        if len(subset) == 0: continue
        subset.to_csv(photometry_folder+'observations/leftovers_'+str(idx)+'.dat', index = False, sep = '\t')
        
    # Copy and paste the splitted observations into the magphys folder.
    subprocess.call(['rm '+magphys_folder+'left*'], shell = True)
    subprocess.call(['cp -rp '+photometry_folder+'observations/left* '+magphys_folder+'.'], shell = True)

    # Re-run magphys.
    run_magphys_leftovers(galaxy_name, run_type)
    return