#!/bin/sh

# parameters for tidying up
yearFrom=$1
yearTo=$2
# Options for this are TOM5
version=$3

# Get paramters as specificed in setUpData.dat
read -r -a parms < tidy_parms 

spinupStart=${parms[0]}
spinupEnd=${parms[1]}
spinupRestartKeepFrequency=${parms[2]}
spinupOutputKeepFrequency=${parms[3]}
runRestartKeepFrequency=${parms[4]}
runOutputKeepFrequency=${parms[5]}
keepGrid_T=${parms[6]}
keepDiad=${parms[7]}
keepPtrc=${parms[8]}
keepIce=${parms[9]}
keepGrid_V=${parms[10]}

echo $spinupStart
echo $spinupEnd
echo $spinupRestartKeepFrequency
echo $spinupOutputKeepFrequency
echo $runRestartKeepFrequency
echo $runOutputKeepFrequency
echo $keepGrid_T
echo $keepDiad
echo $keepPtrc
echo $keepIce
echo $keepGrid_V

pointsPerYear=5475

model=$(basename "$PWD")
id=${model: -4}
echo "Tidying up year: $1 to $2 for $model"
baseDir="/gpfs/data/greenocean/software/runs/"

# make central area
if [ ! -d $baseDir$model ]; then
    echo "creating directory"
    mkdir $baseDir$model
fi

echo "Coping data centrally"
for (( y=$yearFrom; y<=$yearTo; y++ ))
do

    # perform the totalling up before anything else
    # ln -s ORCA2_1m_${y}0101_${y}1231_grid_T.nc ${id}_${y}_grid.nc
    # ln -s ORCA2_1m_${y}0101_${y}1231_diad_T.nc ${id}_${y}_dia2d.nc
    # ln -s ORCA2_1m_${y}0101_${y}1231_diad_T.nc ${id}_${y}_dia3d.nc
    # ln -s ORCA2_1m_${y}0101_${y}1231_ptrc_T.nc ${id}_${y}_ptrc.nc
    # echo $version > total.arg
    # echo $id >> total.arg
    # echo $y >> total.arg
    # # total up the outputs and append text files.
    # totaltrd_ada
    # # copy the totals over
    # cp total* ${baseDir}${model}/
    # rm -f ${id}_${y}_grid.nc
    # rm -f ${id}_${y}_dia2d.nc
    # rm -f ${id}_${y}_dia3d.nc
    # rm -f ${id}_${y}_ptrc.nc

    # python based totalling
    ./breakdown.py breakdown_parms ${y} ${y}

    #  process output files
    if [[ $y < $spinupEnd ]]; then
    # years from start
        since=$((y-spinupStart))
        remainder=$(( since % spinupOutputKeepFrequency))
    else
        since=$((y-spinupEnd))
        remainder=$(( since % runOutputKeepFrequency))
    fi
    
    echo $points $remainder $since

    if [[ "$remainder"  -eq 0 ]]; then

        echo "copying output $y"
        if [[ $keepGrid_T -eq 1 ]]; then cp ORCA2_1m_${y}0101_${y}1231_grid_T.nc $baseDir$model; fi
        if [[ $keepDiad -eq 1 ]]; then cp ORCA2_1m_${y}0101_${y}1231_diad_T.nc $baseDir$model; fi
        if [[ $keepPtrc -eq 1 ]]; then cp ORCA2_1m_${y}0101_${y}1231_ptrc_T.nc $baseDir$model; fi
        if [[ $keepIce -eq 1 ]]; then cp ORCA2_1m_${y}0101_${y}1231_icemod.nc $baseDir$model; fi
        if [[ $keepGrid_V -eq 1 ]]; then cp ORCA2_1m_${y}0101_${y}1231_grid_V.nc $baseDir$model; fi

        echo "Deleting local data if exists centrally"
        rm -f ORCA2_1m_${y}0101_${y}1231_grid_*.nc
        rm -f ORCA2_1m_${y}0101_${y}1231_diad_T.nc
        rm -f ORCA2_1m_${y}0101_${y}1231_ptrc_T.nc
        rm -f ORCA2_1m_${y}0101_${y}1231_icemod.nc

        if [[ $keepGrid_T -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_1m_${y}0101_${y}1231_grid_T.nc; fi
        if [[ $keepDiad -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_1m_${y}0101_${y}1231_diad_T.nc; fi
        if [[ $keepPtrc -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_1m_${y}0101_${y}1231_ptrc_T.nc; fi
        if [[ $keepIce -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_1m_${y}0101_${y}1231_icemod.nc; fi
        if [[ $keepGrid_V -eq 1 ]]; then ln -s ${baseDir}${model}/ORCA2_1m_${y}0101_${y}1231_grid_V.nc; fi
    else
        rm -f ORCA2_1m_${y}0101_${y}1231_grid_*.nc
        rm -f ORCA2_1m_${y}0101_${y}1231_diad_T.nc
        rm -f ORCA2_1m_${y}0101_${y}1231_ptrc_T.nc
        rm -f ORCA2_1m_${y}0101_${y}1231_icemod.nc
    fi

    # process restart files
    if [[ $y < $spinupEnd ]]; then
    # years from start
        since=$((y-spinupStart))
        remainder=$(( since % spinupRestartKeepFrequency))
    else
        since=$((y-spinupEnd))
        remainder=$(( since % runRestartKeepFrequency))
	# reset since to correctly count points for linking and deleting purposes
	since=$((y-spinupStart))
    fi
    
    points=$((pointsPerYear * since))
    # this is concerned with the previous year's restart files, meaning no 
    # restarts are deleted until the next year is complete
    echo $points $remainder $since

    if [[ "$remainder"  -eq 0 ]]; then
        echo "copying restart $y"
        cp ORCA2_*${points}_restart_*.nc $baseDir$model
        rm -f ORCA2_*${points}_restart_*.nc
        ls -1 ${baseDir}${model}/ORCA2_*${points}_restart_*.nc | awk '{print "ln -s "$1 }' | bash
    else
        echo "Deleting local data"
        rm -f ORCA2_*${points}_restart_*.nc
    fi




done
