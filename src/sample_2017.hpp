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

namelist systematicList {
    "syst_elec_reco",
    "syst_elec_id",
    "syst_muon_id",
    "syst_muon_iso",
    "syst_pu",
    "syst_b_correlated",
    "syst_b_uncorrelated",
    "syst_l_correlated",
    "syst_l_uncorrelated",
    "syst_pt_top",
    "syst_prefiring",
    "syst_em_trig",
    "syst_ps_isr",
    "syst_ps_fsr",
    "syst_qcdscale",
    "syst_pdfas"
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
    "1.04",
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

namelist jecList_2017 {
    "Total_up",
    "Total_down",
    "Absolute_up",
    "Absolute_down",
    "Absolute_2017_up",
    "Absolute_2017_down",
    "FlavorQCD_up",
    "FlavorQCD_down",
    "BBEC1_up",
    "BBEC1_down",
    "BBEC1_2017_up",
    "BBEC1_2017_down",
    "RelativeBal_up",
    "RelativeBal_down",
    "RelativeSample_2017_up",
    "RelativeSample_2017_down"
};

namelist sampleList_MC_2017 {
    "MC_dibosons_WW",
    "MC_dibosons_WZ",
    "MC_dibosons_ZZ",
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
    "MC_wjets_WJets",
    "MC_wjets_WJets2",
    "MC_zjets_DY_1050",
    "MC_zjets_DY_50",
    "MC_zjets_DY_502"
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
    "alt_MC_mtop175p5_pmx",
    "alt_MChad_CP5Down",
    "alt_MChad_CP5Up",
    "alt_MChad_GluonMove",
    "alt_MChad_QCDbased",
    "alt_MChad_erdOn",
    "alt_MChad_hdampDown",
    "alt_MChad_hdampUp",
    "alt_MChad_mtop169p5",
    "alt_MChad_mtop175p5_pmx",
    "alt_MCsemilep_CP5Down",
    "alt_MCsemilep_CP5Up",
    "alt_MCsemilep_GluonMove",
    "alt_MCsemilep_QCDbased",
    "alt_MCsemilep_erdOn",
    "alt_MCsemilep_hdampDown",
    "alt_MCsemilep_hdampUp",
    "alt_MCsemilep_mtop169p5",
    "alt_MCsemilep_mtop169p5_pmx",
    "alt_MCsemilep_mtop175p5",
    "alt_MCsemilep_mtop175p5_pmx"
};

std::vector<double> number_of_events_2017 {
    7765828,
    3057400,
    1925931,
    69155808,
    130262440,
    64412112,
    9103042,
    7745276,
    7945242,
    3675910,
    4971954,
    9903448,
    9903448,
    811306,
    19024650,
    19024650,
    750000,
    77700506,
    77700506,
    39536839,
    186329507,
    186329507
};

std::vector<double> mc_rescale_2017 {
    0.634777257,
    0.640187381435,
    0.356295313799,
    0.000741891445048,
    0.000384013200586,
    0.000786496035051,
    0.0125975934813,
    0.0054987272033,
    0.00537023826546,
    0.914563604658,
    1.13615504086,
    0.00249455236929,
    0.00249455236929,
    0.0301071397036,
    0.00227191128268,
    0.00227191128268,
    0.0574278760567,
    32.9146594319,
    32.9146594319,
    23.7951788546,
    7.77598770742e-05,
    7.77598770742e-05
};

std::vector<double> alt_mc_rescale_2017 {
    0.00878306761891,
    0.00878306761891,
    0.00785540671835,
    0.00785540671835,
    0.0093941424701,
    0.00331023637607,
    0.00331023637607,
    0.0100765947355,
    0.00332308015388,
    0.00332308015388,
    0.00389827953276,
    0.00389827953276,
    0.00322904290599,
    0.00322904290599,
    0.00380413022256,
    0.00380413022256,
    0.00183550378008,
    0.00184525974403,
    0.00191859348332,
    0.00192550691559,
    0.00193724245117,
    0.00184460449991,
    0.00183491936006,
    0.00129926910718,
    0.00129926910718,
    0.00310307149414,
    0.00276646766535,
    0.00193013079237,
    0.00186758520179,
    0.0476021635796,
    0.00257384888712,
    0.00242160171032,
    0.00175284892958,
    0.00175284892958,
    0.00232322070763,
    0.00232322070763
};

std::vector<std::vector<double> > jec_mc_rescale_2017 {
    {
        0.634777257,
        0.640187381435,
        0.356295313799,
        0.000741891445048,
        0.000384013200586,
        0.000786496035051,
        0.0125975934813,
        0.0054987272033,
        0.00537023826546,
        0.914563604658,
        1.13615504086,
        0.00249455236929,
        0.00249455236929,
        0.0301071397036,
        0.00227191128268,
        0.00227191128268,
        0.0574278760567,
        32.9146594319,
        32.9146594319,
        23.7951788546,
        7.77598770742e-05,
        7.77598770742e-05
    },
    {
        0.634777257,
        0.640187381435,
        0.356295313799,
        0.000741891445048,
        0.000384013200586,
        0.000786496035051,
        0.0125975934813,
        0.0054987272033,
        0.00537023826546,
        0.914563604658,
        1.13615504086,
        0.00249455236929,
        0.00249455236929,
        0.0301071397036,
        0.00227191128268,
        0.00227191128268,
        0.0574278760567,
        32.9146594319,
        32.9146594319,
        23.7951788546,
        7.77598770742e-05,
        7.77598770742e-05
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

