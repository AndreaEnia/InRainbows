#!/bin/bash

export galaxy_path=$1/run_$2/

# Run the obs_zero, building the galaxy libraries
# then wait TOT seconds for this to be done
sed -i '$ d' .magphys_tcshrc
sed -i '$ d' .magphys_tcshrc
sed -i '$ a setenv USER_FILTERS   '$galaxy_path'filters.dat' .magphys_tcshrc
sed -i '$ a setenv USER_OBS       '$galaxy_path'obs_zero.dat' .magphys_tcshrc
screen -s tcsh -d -m csh run_magphys_zero.csh
sleep 300

# Run the splitted observations_**_**.dat
cd $galaxy_path
for obs in `find -name "obse*"`; do
    echo $obs
    cd '../..'
    sed -i '$ d' .magphys_tcshrc
    sed -i '$ a setenv USER_OBS       '$galaxy_path${obs#*/} .magphys_tcshrc
    screen -s tcsh -d -m csh run_magphys.csh
    sleep 2
    cd $galaxy_path
    done
