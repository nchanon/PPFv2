// 2016 samples : 

#include <vector> 
#include <string> 
 
using  namelist = std::vector<std::string>; 
 
#ifndef COMMON_LIST 
#define COMMON_LIST 
namelist ttbarList {
    "signal",
    "ttx",
    "singletop",
    "dibosons",
    "wjets",
    "zjets"
};

namelist systematicList {
    "syst_elec_reco",
    "syst_elec_id",
    "syst_muon_id",
    "syst_muon_iso",
    "syst_em_trig",
    "syst_pu"
};

namelist systematicTimeList {
    "lumi_inclusive",
    "lumi_stability",
    "lumi_linearity"
};

namelist systematicRate {
    "1.2",
    "1.3",
    "1.3",
    "1.3",
    "1.2"
};

#endif 

namelist triggerList_2016 {
    "trg_muon_electron_mu23ele12_fired",
    "trg_muon_electron_mu23ele12DZ_fired",
    "trg_muon_electron_mu8ele23_fired",
    "trg_muon_electron_mu8ele23DZ_fired",
    "trg_muon_mu24_fired",
    "trg_muon_mutk24_fired",
    "trg_electron_ele27_fired"
};

namelist sampleList_MC_2016 {
    "MC_dibosons_WW",
    "MC_dibosons_WW2",
    "MC_dibosons_WZ",
    "MC_dibosons_WZ2",
    "MC_dibosons_ZZ",
    "MC_dibosons_ZZ2",
    "MC_signal_dilep",
    "MC_signal_hadronic",
    "MC_signal_semilep",
    "MC_singletop_ST_s",
    "MC_singletop_ST_tW_antitop",
    "MC_singletop_ST_tW_top",
    "MC_singletop_ST_t_antitop",
    "MC_singletop_ST_t_top",
    "MC_ttx_TTW",
    "MC_ttx_TTW2",
    "MC_ttx_TTW3",
    "MC_ttx_TTZ",
    "MC_ttx_TTZ2",
    "MC_ttx_TTZ3",
    "MC_ttx_TTZ4",
    "MC_wjets_WJets",
    "MC_wjets_WJets2",
    "MC_zjets_DY_50"
};

std::vector<double> mc_rescale_2016 {
    0.840900712543,
    0.840900712543,
    0.514671934349,
    0.514671934349,
    0.298363410657,
    0.298363410657,
    0.000657611972209,
    0.000631102864445,
    0.000406965344278,
    0.0100761139063,
    0.00739198265276,
    0.00740052052694,
    0.163526347845,
    0.153384581332,
    0.00311306057978,
    0.00311306057978,
    0.00618955069753,
    0.00234533182629,
    0.00234533182629,
    0.00234533182629,
    0.00491230632022,
    5.47500908161e-05,
    5.47500908161e-05,
    0.000117804520985
};

namelist sampleList_DATA_2016 {
    "MuonEG_Run2016B_17Jul2018",
    "MuonEG_Run2016C_17Jul2018",
    "MuonEG_Run2016D_17Jul2018",
    "MuonEG_Run2016E_17Jul2018",
    "MuonEG_Run2016F_17Jul2018",
    "MuonEG_Run2016G_17Jul2018",
    "MuonEG_Run2016H_17Jul2018",
    "SingleElectron_Run2016B_17Jul2018",
    "SingleElectron_Run2016C_17Jul2018",
    "SingleElectron_Run2016D_17Jul2018",
    "SingleElectron_Run2016E_17Jul2018",
    "SingleElectron_Run2016F_17Jul2018",
    "SingleElectron_Run2016G_17Jul2018",
    "SingleElectron_Run2016H_17Jul2018",
    "SingleMuon_Run2016B_17Jul2018",
    "SingleMuon_Run2016C_17Jul2018",
    "SingleMuon_Run2016D_17Jul2018",
    "SingleMuon_Run2016E_17Jul2018",
    "SingleMuon_Run2016F_17Jul2018",
    "SingleMuon_Run2016G_17Jul2018",
    "SingleMuon_Run2016H_17Jul2018"
};

std::vector<double> succedJobs_2016 {
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0
};

namelist data_2016 {
    "Run2016"
};

