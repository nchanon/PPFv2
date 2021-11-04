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
    "syst_pt_top",
    "syst_prefiring"
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
    "trg_muon_electron_mu8ele23_fired",
    "trg_muon_mu24_fired",
    "trg_muon_mutk24_fired",
    "trg_electron_ele27_fired",
    "trg_muon_electron_mu8ele23DZ_fired",
    "trg_muon_electron_mu23ele12DZ_fired"
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
    "alt_MC_mtop175p5",
    "alt_MChad_CP5Down",
    "alt_MChad_CP5Up",
    "alt_MChad_GluonMove",
    "alt_MChad_QCDbased",
    "alt_MChad_erdOn",
    "alt_MChad_hdampDown",
    "alt_MChad_hdampUp",
    "alt_MChad_mtop169p5",
    "alt_MChad_mtop175p5",
    "alt_MCsemilep_CP5Down",
    "alt_MCsemilep_CP5Up",
    "alt_MCsemilep_GluonMove",
    "alt_MCsemilep_QCDbased",
    "alt_MCsemilep_erdOn",
    "alt_MCsemilep_hdampDown",
    "alt_MCsemilep_hdampUp",
    "alt_MCsemilep_mtop169p5",
    "alt_MCsemilep_mtop175p5"
};

std::vector<double> mc_rescale_2016 {
    0.533846696909,
    0.533846696909,
    0.423248767814,
    0.423248767814,
    0.298363410657,
    0.298363410657,
    0.0035400931205,
    0.000664652750811,
    0.00083920145827,
    0.0105939647742,
    0.00739198265276,
    0.00740052052694,
    0.163526347845,
    0.153384581332,
    0.00410526703258,
    0.00410526703258,
    0.0256093460723,
    0.00261284391532,
    0.00261284391532,
    0.00261284391532,
    0.047979475341,
    0.000607119104493,
    0.000277344617309,
    0.000277344617309,
    0.000277344617309,
    0.000127799473009
};

std::vector<double> alt_mc_rescale_2016 {
    0.00308679827397,
    0.00308679827397,
    0.00298863599549,
    0.00298863599549,
    0.00295894094832,
    0.00298700840922,
    0.00307288618321,
    0.00307288618321,
    0.00297464153488,
    0.00297464153488,
    0.00297853098459,
    0.00297853098459,
    0.00282261782195,
    0.00425593585093,
    0.00154871621744,
    0.0015477093498,
    0.00161011742597,
    0.00173656818093,
    0.0015817217071,
    0.00156411904217,
    0.00161398287655,
    0.00417129316054,
    0.0049023491269,
    0.00157487982356,
    0.00156633462456,
    0.0016544799152,
    0.00149931617346,
    0.0015114282423,
    0.00153021921529,
    0.00152804663569,
    0.00150275260967,
    0.00216875844666
};

std::vector<std::vector<double> > jec_mc_rescale_2016 {
    {
        0.533846696909,
        0.533846696909,
        0.423248767814,
        0.423248767814,
        0.298363410657,
        0.298363410657,
        0.0035400931205,
        0.000664652750811,
        0.00083920145827,
        0.0105939647742,
        0.00739198265276,
        0.00740052052694,
        0.163526347845,
        0.153384581332,
        0.00410526703258,
        0.00410526703258,
        0.0256093460723,
        0.00261284391532,
        0.00261284391532,
        0.00261284391532,
        0.047979475341,
        0.000607119104493,
        0.000277344617309,
        0.000277344617309,
        0.000277344617309,
        0.000127799473009
    },
    {
        0.533846696909,
        0.533846696909,
        0.423248767814,
        0.423248767814,
        0.298363410657,
        0.298363410657,
        0.0035400931205,
        0.000664652750811,
        0.00083920145827,
        0.0105939647742,
        0.00739198265276,
        0.00740052052694,
        0.163526347845,
        0.153384581332,
        0.00410526703258,
        0.00410526703258,
        0.0256093460723,
        0.00261284391532,
        0.00261284391532,
        0.00261284391532,
        0.047979475341,
        0.000607119104493,
        0.000277344617309,
        0.000277344617309,
        0.000277344617309,
        0.000127799473009
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

