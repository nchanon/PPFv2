#include <vector>
#include <string>

using  namelist = std::vector<std::string>;

namelist triggerList{
    "trg_muon_electron_mu8ele23DZ_fired",
    "trg_muon_electron_mu23ele12_fired",
    "trg_muon_mu27_fired",
    "trg_electron_ele35_fired"
};

namelist ttbarList{
    "signal",
    "ttx",
    "singletop",
    "dibosons",
    "wjets",
    "zjets"
};

namelist data2017{
    "Run2017"
};

namelist sampleList_MC{
    "MC_signal_dilep",
    "MC_signal_hadronic",
    "MC_signal_semilep",
    //
    "MC_ttx_TTW",
    "MC_ttx_TTW2",
    "MC_ttx_TTW3",
    "MC_ttx_TTZ",
    "MC_ttx_TTZ2",
    "MC_ttx_TTZ3",
    //
    "MC_singletop_ST_s",
    "MC_singletop_ST_s2",
    "MC_singletop_ST_t_antitop",
    "MC_singletop_ST_t_top",
    "MC_singletop_ST_tW_antitop",
    "MC_singletop_ST_tW_antitop2",
    "MC_singletop_ST_tW_top",
    "MC_singletop_ST_tW_top2",
    //
    "MC_dibosons_WW",
    "MC_dibosons_WZ",
    "MC_dibosons_ZZ",
    //
    "MC_wjets_WJets",
    "MC_wjets_WJets2",
    "MC_zjets_DY_1050",
    "MC_zjets_DY_50",
    "MC_zjets_DY_502"
};

std::vector<double> effectiveN0{
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1,
    1
};

namelist sampleList_DATA{
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