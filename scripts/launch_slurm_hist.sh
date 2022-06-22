#!/bin/bash

for year in 2016
do
  for obs in n_bjets
  do

    #data
    sbatch scripts/slurm_HistCreator_data_all.sh ${obs} ${year}

    #inclusive
    #for puchoice in `seq -3 -2`
    for puchoice in -2
    do
      sbatch scripts/slurm_HistCreator_mc_inclusive.sh ${obs} ${year} ${puchoice}
    done
    
    #differential
    sbatch scripts/slurm_HistCreator_mc_differential.sh ${obs} ${year} -2
    for puchoice in `seq 0 23`
    do
      sbatch scripts/slurm_HistCreator_mc_differential.sh ${obs} ${year} ${puchoice}
    done

  done
done
