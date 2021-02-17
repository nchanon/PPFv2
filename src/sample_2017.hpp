// 2017 samples : 

#include <vector> 
#include <string> 
 
using  namelist = std::vector<std::string>; 
 
#ifndef COMMON_LIST 
#define COMMON_LIST 
namelist triggerList {
    "trg_muon_electron_mu8ele23DZ_fired",
    "trg_muon_electron_mu23ele12_fired",
    "trg_muon_mu27_fired",
    "trg_electron_ele35_fired"
};

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
    "1.3",
    "1.2",
    "1.5",
    "1.5",
    "1.5"
};

#endif 

namelist sampleList_MC_2017 {
    "MC_dibosons_WW",
    "MC_dibosons_WZ",
    "MC_dibosons_ZZ",
    "MC_signal_dilep",
    "MC_signal_hadronic",
    "MC_signal_semilep",
    "MC_singletop_ST_s",
    "MC_singletop_ST_s2",
    "MC_singletop_ST_tW_antitop",
    "MC_singletop_ST_tW_antitop2",
    "MC_singletop_ST_tW_top",
    "MC_singletop_ST_tW_top2",
    "MC_singletop_ST_t_antitop",
    "MC_singletop_ST_t_top",
    "MC_ttx_TTW",
    "MC_ttx_TTW2",
    "MC_ttx_TTW3",
    "MC_ttx_TTZ",
    "MC_ttx_TTZ2",
    "MC_ttx_TTZ3",
    "MC_wjets_WJets",
    "MC_wjets_WJets2",
    "MC_zjets_DY_1050",
    "MC_zjets_DY_50",
    "MC_zjets_DY_502"
};

std::vector<double> mc_rescale_2017 {
    2.9838524718,
    1.76815199913,
    0.356295313799,
    0.000741891445048,
    0.000386197572497,
    0.000710616081968,
    0.00595755051361,
    0.00595755051361,
    0.00268170342106,
    0.00268170342106,
    0.00268388061318,
    0.00268388061318,
    0.914563604658,
    0.954264128071,
    0.00325732310809,
    0.00325732310809,
    0.00647638104016,
    0.00209806378462,
    0.00209806378462,
    0.00439440247811,
    0.000217302970276,
    0.000217302970276,
    102.645679589,
    7.38985888724e-05,
    7.38985888724e-05
};

namelist sampleList_DATA_2017 {
    "MuonEG_Run2017B_31Mar2018",
    "MuonEG_Run2017C_31Mar2018",
    "MuonEG_Run2017D_31Mar2018",
    "MuonEG_Run2017E_31Mar2018",
    "MuonEG_Run2017F_31Mar2018",
    "SingleElectron_Run2017B_31Mar2018",
    "SingleElectron_Run2017C_31Mar2018",
    "SingleElectron_Run2017D_31Mar2018",
    "SingleElectron_Run2017E_31Mar2018",
    "SingleElectron_Run2017F_31Mar2018",
    "SingleMuon_Run2017B_31Mar2018",
    "SingleMuon_Run2017C_31Mar2018",
    "SingleMuon_Run2017D_31Mar2018",
    "SingleMuon_Run2017E_31Mar2018",
    "SingleMuon_Run2017F_31Mar2018"
};

std::vector<double> succedJobs_2017 {
    1.0,
    1.0,
    1.0,
    1.04384133612,
    1.0,
    1.14155251142,
    1.0,
    1.0,
    1.20048019208,
    1.30718954248,
    1.0,
    1.0,
    1.0,
    1.0,
    1.0
};

namelist data_2017 {
    "Run2017"
};

