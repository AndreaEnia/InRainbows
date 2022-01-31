#!/bin/bash

export galaxy_path=$1/run_$2/leftovers/

# Run the splitted observations_**_**.dat
cd $galaxy_path
for obs in `find -name "left*"`; do
    echo $obs
    cd '../../../'
    sed -i '$ d' .magphys_tcshrc
    sed -i '$ a setenv USER_OBS       '$galaxy_path${obs#*/} .magphys_tcshrc
    screen -s tcsh -d -m csh run_magphys.csh
    sleep 1
    cd $galaxy_path
    done
