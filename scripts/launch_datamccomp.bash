#!/bin/bash

for i in 2016 2017
#for i in 2016
do
#    python ./bin/data_mc_comparaison_splitted.py m_dilep ${i} "Dilepton mass (GeV)" $1
#    python ./bin/data_mc_comparaison_splitted.py eta_elec ${i} "Electron #eta" $1
#    python ./bin/data_mc_comparaison_splitted.py pt_elec ${i} "Electron p_{T} (GeV)" $1
#    python ./bin/data_mc_comparaison_splitted.py eta_muon ${i} "Muon #eta" $1
#    python ./bin/data_mc_comparaison_splitted.py pt_muon ${i} "Muon p_{T} (GeV)" $1
#    python ./bin/data_mc_comparaison_splitted.py n_bjets ${i} "number of b jets" $1
#    python ./bin/data_mc_comparaison_splitted.py n_jets ${i} "number of jets" $1
#    python ./bin/data_mc_comparaison_splitted.py j1_pt ${i} "First jet p_{T} (GeV)" $1
#    python ./bin/data_mc_comparaison_splitted.py j1_eta ${i} "First jet #eta" $1
#    python ./bin/data_mc_comparaison_splitted.py b1_pt ${i} "First b-jet p_{T} (GeV)" $1
#    python ./bin/data_mc_comparaison_splitted.py b1_eta ${i} "First b-jet #eta" $1
#    python ./bin/data_mc_comparaison_splitted.py j2_pt ${i} "Second jet p_{T} (GeV)" $1
#    python ./bin/data_mc_comparaison_splitted.py j2_eta ${i} "Second jet #eta" $1
#    python ./bin/data_mc_comparaison_splitted.py pt_emu ${i} "Dilepton p_{T} (GeV)" $1

    python ./bin/mc_comparaison.py pt_emu ${i} "Dilepton p_{T} (GeV)" $1
    python ./bin/mc_comparaison.py m_dilep ${i} "Dilepton mass (GeV)" $1
    python ./bin/mc_comparaison.py n_bjets ${i} "number of b jets" $1

done

#python ./bin/data_mc_comparaison.py m_dilep 2016 "Dilepton mass (GeV)" $1
#python ./bin/data_mc_comparaison.py m_dilep 2017 "Dilepton mass (GeV)" $1
#python ./bin/data_mc_comparaison.py eta_elec 2016 "Electron #eta" $1 
#python ./bin/data_mc_comparaison.py eta_elec 2017 "Electron #eta" $1
#python ./bin/data_mc_comparaison.py pt_elec 2016 "Electron p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py pt_elec 2017 "Electron p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py eta_muon 2016 "Muon #eta" $1
#python ./bin/data_mc_comparaison.py eta_muon 2017 "Muon #eta" $1
#python ./bin/data_mc_comparaison.py pt_muon 2016 "Muon p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py pt_muon 2017 "Muon p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py n_bjets 2016 "number of b-jets" $1
#python ./bin/data_mc_comparaison.py n_jets 2017 "number of b-jets" $1
#python ./bin/data_mc_comparaison.py n_jets 2016 "number of jets" $1
#python ./bin/data_mc_comparaison.py n_bjets 2017 "number of jets" $1
#python ./bin/data_mc_comparaison.py j1_pt 2016 "First jet p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py j1_pt 2017 "First jet p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py j1_eta 2016 "First jet #eta" $1
#python ./bin/data_mc_comparaison.py j1_eta 2017 "First jet #eta" $1
#python ./bin/data_mc_comparaison.py b1_pt 2016 "First b-jet p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py b1_pt 2017 "First b-jet p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py b1_eta 2016 "First b-jet #eta" $1
#python ./bin/data_mc_comparaison.py b1_eta 2017 "First b-jet #eta" $1
#python ./bin/data_mc_comparaison.py j2_pt 2016 "Second jet p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py j2_pt 2017 "Second jet p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py j2_eta 2016 "Second jet #eta" $1
#python ./bin/data_mc_comparaison.py j2_eta 2017 "Second jet #eta" $1
#python ./bin/data_mc_comparaison.py pt_emu 2016 "Dilepton p_{T} (GeV)" $1
#python ./bin/data_mc_comparaison.py pt_emu 2017 "Dilepton p_{T} (GeV)" $1

