// 2017 samples : 

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
    "vjets"
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
    "1.3"
};

#endif 

namelist triggerList_2017 {
    "trg_muon_electron_mu8ele23DZ_fired",
    "trg_muon_electron_mu23ele12_fired",
    "trg_muon_mu27_fired",
    "trg_electron_ele35_fired"
};

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
    "MC_zjets_DY_50"
};

namelist sampleList_ALT_2017 {
    "alt_MC_CP5Down",
    "alt_MC_CP5Down_pmx",
    "alt_MC_CP5Up",
    "alt_MC_CP5Up_pmx",
    "alt_MC_GluonMove",
    "alt_MC_QCDbased",
    "alt_MC_QCDbased_ext",
    "alt_MC_erdOn",
    "alt_MC_hdampDown",
    "alt_MC_hdampDown_pmx",
    "alt_MC_hdampUp",
    "alt_MC_hdampUp_pmx",
    "alt_MC_mtop169p5",
    "alt_MC_mtop169p5_pmx",
    "alt_MC_mtop175p5",
    "alt_MC_mtop175p5_pmx"
};

std::vector<double> mc_rescale_2017 {
    1.02264407904,
    0.785531524662,
    0.356295313799,
    0.000772359005614,
    0.00369278988411,
    0.000880993019701,
    0.00626403787079,
    0.00626403787079,
    0.0027541587068,
    0.0027541587068,
    0.00273496196486,
    0.00273496196486,
    0.940173472641,
    0.965740972171,
    0.0026728219182,
    0.0026728219182,
    0.0301071397037,
    0.00227191128267,
    0.00227191128267,
    0.0574278760564,
    32.9146594319,
    32.9146594319,
    25.487904781,
    0.000642188670664
};

std::vector<double> alt_mc_rescale_2017 {
    0.00339552026081,
    0.00339552026081,
    0.00331173759109,
    0.00331173759109,
    0.00685274348143,
    0.00332763984557,
    0.00332763984557,
    0.0100765946902,
    0.00332308013899,
    0.00332308013899,
    0.00392992624072,
    0.00392992624072,
    0.00174653736358,
    0.00174653736358,
    0.00174653736358,
    0.00174653736358
};

std::vector<std::vector<double> > jec_mc_rescale_2017 {
    {
        0.76718315635,
        1.48419665657,
        1.92020005346,
        0.000772343615678,
        0.000870029459444,
        0.000213812509647,
        0.0540505288261,
        0.0540505288261,
        0.00684108635597,
        0.00684108635597,
        0.0067934033488,
        0.0067934033488,
        1.03422828283,
        0.632241728871,
        1.16500268801,
        1.16500268801,
        6.60016585059,
        0.799959147657,
        0.799959147657,
        9.65423912177,
        0.0476377250632,
        0.0476377250632,
        0.100271379186,
        9.18587747342e-06
    },
    {
        0.76718315635,
        1.48419665657,
        1.92020005346,
        0.000772343615678,
        0.000870029459444,
        0.000213812509647,
        0.0540505288261,
        0.0540505288261,
        0.00684108635597,
        0.00684108635597,
        0.0067934033488,
        0.0067934033488,
        1.03422828283,
        0.632241728871,
        1.16500268801,
        1.16500268801,
        6.60016585059,
        0.799959147657,
        0.799959147657,
        9.65423912177,
        0.0476377250632,
        0.0476377250632,
        0.100271379186,
        9.18587747342e-06
    }
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

namelist data_2017 {
    "Run2017"
};

