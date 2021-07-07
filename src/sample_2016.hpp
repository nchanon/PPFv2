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

namelist jecList {
    "TotalUp",
    "TotalDown"
};

namelist systematicList {
    "syst_elec_reco",
    "syst_elec_id",
    "syst_muon_id",
    "syst_muon_iso",
    "syst_pu",
    "syst_b",
    "syst_pt_top"
};

namelist systematicAltList {
    "CP5Up",
    "CP5Down",
    "hdampUp",
    "hdampDown",
    "mtop169",
    "mtop175",
    "erdOn",
    "QCD",
    "GluonMove"
};

namelist systematicTimeList {
    "lumi_flat",
    "lumi_stability",
    "lumi_linearity",
    "emu_trig"
};

namelist systematicRate {
    "1.05",
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
    "MC_zjets_DY_1050",
    "MC_zjets_DY_10502",
    "MC_zjets_DY_10503",
    "MC_zjets_DY_50"
};

namelist sampleList_ALT_2016 {
    "alt_MC_CP5Down",
    "alt_MC_CP5Down_pmx",
    "alt_MC_CP5Up",
    "alt_MC_CP5Up_pmx",
    "alt_MC_GluonMove",
    "alt_MC_QCDbased",
    "alt_MC_erdOn",
    "alt_MC_erdOn_ext",
    "alt_MC_hdampDown",
    "alt_MC_hdampDown_pmx",
    "alt_MC_hdampUp",
    "alt_MC_hdampUp_pmx",
    "alt_MC_mtop169p5",
    "alt_MC_mtop175p5"
};

std::vector<double> mc_rescale_2016 {
    0.533846696909,
    0.533846696909,
    0.423248767814,
    0.423248767814,
    0.298363410657,
    0.298363410657,
    0.00104418753183,
    0.000631102864445,
    0.000416266269225,
    0.0100761139063,
    0.00739198265276,
    0.00740052052694,
    0.163526347845,
    0.155780372634,
    0.00311306057978,
    0.00311306057978,
    0.00618955069753,
    0.00234533182629,
    0.00234533182629,
    0.00234533182629,
    0.00491230632022,
    5.47500908161e-05,
    5.47500908161e-05,
    5.32604427151e-05,
    0.000193651424377,
    0.000193651424377,
    0.000428329396525
};

std::vector<double> alt_mc_rescale_2016 {
    0.0030867982601,
    0.0030867982601,
    0.00298863598207,
    0.00298863598207,
    0.00302745295971,
    0.0029870083958,
    0.00304513167774,
    0.00304513167774,
    0.00297464152152,
    0.00297464152152,
    0.00302421492156,
    0.00302421492156,
    0.00171083151691,
    0.00171083151691
};

std::vector<std::vector<double> > jec_mc_rescale_2016 {
    {
        0.400489478536,
        0.400489478536,
        0.799693438601,
        0.799693438601,
        1.60798476687,
        1.60798476687,
        0.00104418752714,
        0.000148689229888,
        0.000101025699085,
        0.0869438046804,
        0.0183610303736,
        0.0183822376971,
        0.179885498635,
        0.101984750525,
        1.35688948025,
        1.35688948025,
        1.35688948025,
        0.82581114106,
        0.82581114106,
        0.82581114106,
        0.82581114106,
        7.92403694431e-08,
        7.92403694431e-08,
        7.61838262353e-07,
        7.61838262353e-07,
        7.61838262353e-07,
        1.68507783618e-06
    },
    {
        0.400489478536,
        0.400489478536,
        0.799693438601,
        0.799693438601,
        1.60798476687,
        1.60798476687,
        0.00104418752714,
        0.000148689229888,
        0.000101025699085,
        0.0869438046804,
        0.0183610303736,
        0.0183822376971,
        0.179885498635,
        0.101984750525,
        1.35688948025,
        1.35688948025,
        1.35688948025,
        0.82581114106,
        0.82581114106,
        0.82581114106,
        0.82581114106,
        7.92403694431e-08,
        7.92403694431e-08,
        7.61838262353e-07,
        7.61838262353e-07,
        7.61838262353e-07,
        1.68507783618e-06
    }
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

