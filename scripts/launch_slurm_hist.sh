#!/bin/bash

for year in 2016 2017
do
  #for obs in m_dilep pt_emu
  for obs in n_bjets
  do

    #data
    #sbatch scripts/slurm_HistCreator_data_all.sh ${obs} ${year}

    #inclusive
    #for puchoice in `seq -3 -1`
    for puchoice in -2
    do
      sbatch scripts/slurm_HistCreator_mc_inclusive.sh ${obs} ${year} ${puchoice}
      #python bin/color_reco.py ${obs} ${year} inclusive ${puchoice}
    done
    
    #differential
    #sbatch scripts/slurm_HistCreator_mc_differential.sh ${obs} ${year} -3
    sbatch scripts/slurm_HistCreator_mc_differential.sh ${obs} ${year} -2
    #python bin/color_reco.py ${obs} ${year} timed -2
    #sbatch scripts/slurm_HistCreator_mc_differential.sh ${obs} ${year} -1
    for puchoice in `seq 0 23`
    #for puchoice in -3 12 14 19 1 20
    do
      sbatch scripts/slurm_HistCreator_mc_differential.sh ${obs} ${year} ${puchoice}
      #python bin/color_reco.py ${obs} ${year} timed ${puchoice}
    done

  done
done
