#!/bin/bash
#
#SBATCH --ntasks=1

#export EOS_MGM_URL=root://lyoeos.in2p3.fr

#srun -n1 --exclusive ./runLHEtoNanoGEN.sh root://lyoeos.in2p3.fr//home/nchanon/TTbarNanoGEN/Test/PROC_SM_ttbar_emu_13TeV_LHE_test.lhe &
#srun -n1 --exclusive ./runLHEtoNanoGEN.sh /gridgroup/cms/nchanon/CMSSW_10_6_27/src/TTbarLO_test/slurm/PROC_SM_ttbar_emu_13TeV_LHE_test.lhe &

#for i in `seq 0 9`
#do
#  srun -n1 --exclusive ./runLHEtoNanoGEN.sh PROC_SM_ttbar_emu_13TeV_LHE_${i} &
#done


#arg 1: observable
#arg 2: year
#arg 3: pu choice

#srun -n1 --exclusive ./scripts/runHistCreator_inclusive.sh ${1} ${2} ${3} &

srun -n1 --exclusive ./bin/histograms_creator ${1} ${2} dataAll -1 &

wait
