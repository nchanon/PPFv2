#include "generator.hpp"
#include "debug.h"

#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <string>
#include <TCanvas.h>
#include <TLeaf.h>
#include <TString.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <ROOT/RDataFrame.hxx>

//using namespace ROOT;

using namespace std;

const double OMEGA_GMST = 7.2722e-5;
const double OMEGA_UTC  = 7.2921e-5;
const double T0_2016    = 1451606400.;
const double T0_2017    = 1483228800.;
const double PHASE      = 3.2830;


/////////////////////////////////
// Constructor
/////////////////////////////////

Generator::Generator(std::string      const& observable_p,
                     std::vector<int> const& binning_p,
                     std::string      const& year_p
                    ) 
    : observable(observable_p), nBin(binning_p[0]), minBin(binning_p[1]), maxBin(binning_p[2]),  year(year_p)
{

    //doLoop = true;
    doLoop = false;

/*
    std::ifstream infile("./inputs/timed/LumiData_"+year+".csv");
    if(!infile.is_open()){
        std::cout << "Problem with lumi csv file" << std::endl;
    }
    else{
        std::string line;
        do
        {
            std::getline(infile, line);
            if(line[0] == '#') continue;
            std::istringstream parse( line );
            while(parse)
            {
                int count = 0;
                std::string s;
                while(std::getline( parse, s, ',' )){
                    if(count == 2)
                        timestamp.push_back(utcConverter(s));
                    if(count == 6){
                        instLumi.push_back(std::stod(s));
                    }
                    count++;
                }
            }
        }
        while(!infile.eof());
    }
    for(size_t i = 0; i < 10; ++i){
        
        std::cout << i << " : " << timestamp[i] << " -> " << instLumi[i] << std::endl;
    }
*/    

}

/////////////////////////////////
// Private methods
/////////////////////////////////

time_t Generator::utcConverter(std::string const& time)
{
    //std::istringstream parse(time);
    std::tm t{};
    std::istringstream ss(time);
    ss >> std::get_time(&t, "%D %X");
    if (ss.fail()) {
        throw std::runtime_error{"failed to parse time string"};
    }
    t.tm_year += 2000-1900; // c normalisation
    t.tm_hour += 1; // summer hour normalisation
    std::time_t time_stamp = mktime(&t);
    return time_stamp;
}


double Generator::siderealHour(double time_p)
{
    if(year == "2017")
        return (OMEGA_UTC * (time_p - T0_2016) + PHASE)/OMEGA_GMST;
    else if(year == "2016")
        return (OMEGA_UTC * (time_p - T0_2016) + PHASE)/OMEGA_GMST;
    else
        return 0;
}

double Generator::luminositySumOfWeight(TTree *tree_p)
{
    double sumoflumi = 0;
    for(int i = 0; i < tree_p->GetEntriesFast(); ++i){
        tree_p->GetEntry(i);
	sumoflumi += tree_p->GetLeaf("lumi")->GetValue();
    }
    double lumiavg = sumoflumi / tree_p->GetEntriesFast();

    return lumiavg;
}

double Generator::luminosityCorrection(TTree *tree_p, double lumiavg)
{
    double w = lumiavg/tree_p->GetLeaf("lumi")->GetValue();

    return w;
}

void Generator::drawHisto1D(TTree* tree, std::string obs, std::string string_eventSelection, std::string string_weight, std::string string_triggered, TH1F* hist){

    std::string string_cut = "(" + string_eventSelection + ")*" + string_weight + "*" + string_triggered;
    std::cout << "Cut: " << string_cut<<std::endl;
    std::string string_obs_redirected = obs + " >> " + hist->GetName();

    tree->Draw(string_obs_redirected.c_str(), string_cut.c_str());

    return;
}

void Generator::drawHisto2D(TTree* tree, std::string obs1, std::string obs2, std::string string_eventSelection, std::string string_weight, std::string string_triggered, TH2F* hist){

    std::string string_cut = "(" + string_eventSelection + ")*" + string_weight + "*" + string_triggered;
    //std::cout << "Cut: " << string_cut<<std::endl;
    std::string string_obs_redirected = obs1 + ":" + obs2 + " >> " + hist->GetName();
    tree->Draw(string_obs_redirected.c_str(), string_cut.c_str());

    return;
}

double Generator::generateWeight(TTree *tree_p,  bool isTimed)
{
    double w = 1;
    w *= tree_p->GetLeaf("weight_punew")->GetValue();
    w *= tree_p->GetLeaf("weight_generator")->GetValue();
    w *= tree_p->GetLeaf("weight_top")->GetValue(); //should be applied only to signal, not ttx
    w *= tree_p->GetLeaf("weight_prefiring")->GetValue();
    //lepton
    w *= tree_p->GetLeaf("weight_sfe_id")->GetValue();
    w *= tree_p->GetLeaf("weight_sfe_reco")->GetValue();
    w *= tree_p->GetLeaf("weight_sfm_id")->GetValue();
    w *= tree_p->GetLeaf("weight_sfm_iso")->GetValue();
    //hadron
    w *= tree_p->GetLeaf("weight_sfb")->GetValue();
    //w *= tree_p->GetLeaf("weight_sfl")->GetValue();
    //w *= tree_p->GetLeaf("weight_sfc")->GetValue();

    //trigger
    if(!isTimed)
        w *= tree_p->GetLeaf("weight_sf_em_trig")->GetValue();
    return w;
}

std::string Generator::generateWeightString(bool isTimed, int timebin)
{

  std::string string_weight = "";

  if (timebin==-1) string_weight = "weight_pu";
  else if (timebin==-2) string_weight = "weight_punew";
  else if (timebin==-3) string_weight = "weight_puinc";
  else if (timebin>=0) string_weight = "weight_putime" + std::to_string(timebin);

  string_weight += "*weight_generator";
  string_weight += "*weight_top"; //should be only for ttbar
  string_weight += "*weight_prefiring";
  string_weight += "*weight_sfe_id";
  string_weight += "*weight_sfe_reco";
  string_weight += "*weight_sfm_id";
  string_weight += "*weight_sfm_iso";
  string_weight += "*weight_sfb";

  if (!isTimed)
      string_weight += "*weight_sf_em_trig";

  return string_weight;
}

double Generator::generateSystematics(TTree            * tree_p,
                                      std::string const& systematicName,
                                      bool               isUp
                                     )
{

    double weight = -999.;

    if(systematicName == "syst_pu"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_punew_up")->GetValue(0)/tree_p->GetLeaf("weight_pu")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_punew_down")->GetValue(0)/tree_p->GetLeaf("weight_pu")->GetValue(0);
    }
    else if(systematicName == "syst_pt_top"){
        if(isUp)
            weight = 1;
        else
            weight = 1./tree_p->GetLeaf("weight_top")->GetValue(0);
    }
    //else if(systematicName == "syst_b"){
    //    if(isUp)
    //        weight = tree_p->GetLeaf("weight_sfb_up")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
    //    else
    //        weight = tree_p->GetLeaf("weight_sfb_down")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
    //}
    else if(systematicName == "syst_b_correlated"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_sfb_up_correlated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_sfb_down_correlated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
    }
    else if(systematicName == "syst_b_uncorrelated"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_sfb_up_uncorrelated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_sfb_down_uncorrelated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
    }
    //else if(systematicName == "syst_c"){
    //    if(isUp)
    //        weight = tree_p->GetLeaf("weight_sfc_up")->GetValue(0)/tree_p->GetLeaf("weight_sfc")->GetValue(0);
    //    else
    //        weight = tree_p->GetLeaf("weight_sfc_down")->GetValue(0)/tree_p->GetLeaf("weight_sfc")->GetValue(0);
    //}
    //else if(systematicName == "syst_l"){
    //    if(isUp)
    //        weight = tree_p->GetLeaf("weight_sfl_up")->GetValue(0)/tree_p->GetLeaf("weight_sfl")->GetValue(0);
    //    else
    //        weight = tree_p->GetLeaf("weight_sfl_down")->GetValue(0)/tree_p->GetLeaf("weight_sfl")->GetValue(0);
    //}
    else if(systematicName == "syst_l_correlated"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_sfl_up_correlated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_sfl_down_correlated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
    }
    else if(systematicName == "syst_l_uncorrelated"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_sfl_up_uncorrelated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_sfl_down_uncorrelated")->GetValue(0)/tree_p->GetLeaf("weight_sfb")->GetValue(0);
    }
    else if(systematicName == "syst_prefiring"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_prefiring_up")->GetValue(0)/tree_p->GetLeaf("weight_prefiring")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_prefiring_down")->GetValue(0)/tree_p->GetLeaf("weight_prefiring")->GetValue(0);
    }
    else if(systematicName == "syst_ps_isr"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_ps_variation4")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_ps_variation6")->GetValue(0);
	if (weight==-99 || weight==0) weight = 1;
    }
    else if(systematicName == "syst_ps_fsr"){
        if(isUp)
            weight = tree_p->GetLeaf("weight_ps_variation5")->GetValue(0);
        else
            weight = tree_p->GetLeaf("weight_ps_variation7")->GetValue(0);
	if (weight==-99 || weight==0) weight = 1;
    }
    else{
        double foo = tree_p->GetLeaf(systematicName.c_str())->GetValue(0);
        if(isUp)
            weight = (1+foo);
        else
            weight = (1-foo);
    }

    return weight;
}

std::string Generator::generateSystematicsString(std::string const& systematicName,
                                      	    bool               isUp,
					    int 	       timebin
                                     	   )
{
    std::string string_weight_syst = "";

    if(systematicName == "syst_pu"){
	if (timebin==-1){
            if(isUp) string_weight_syst = "weight_pu_up/weight_pu";
	    else string_weight_syst = "weight_pu_down/weight_pu";
	}
        if (timebin==-2){
            if(isUp) string_weight_syst = "weight_punew_up/weight_punew";
            else string_weight_syst = "weight_punew_down/weight_punew";
        }
        if (timebin==-3){
            if(isUp) string_weight_syst = "weight_puinc_up/weight_puinc";
            else string_weight_syst = "weight_puinc_down/weight_puinc";
        }
	else if (timebin>=0) {
            if(isUp) string_weight_syst = "weight_putime"+std::to_string(timebin)+"_up/weight_putime" + std::to_string(timebin);
            else string_weight_syst = "weight_putime"+std::to_string(timebin)+"_down/weight_putime"+std::to_string(timebin);
	}
    }
    else if(systematicName == "syst_pt_top"){
        if(isUp) string_weight_syst = "1";
	else string_weight_syst = "1./weight_top";
    }
    else if(systematicName == "syst_b_correlated"){
        if(isUp) string_weight_syst = "weight_sfb_up_correlated/weight_sfb";
	else string_weight_syst = "weight_sfb_down_correlated/weight_sfb";
    }
    else if(systematicName == "syst_b_uncorrelated"){
        if(isUp) string_weight_syst = "weight_sfb_up_uncorrelated/weight_sfb";
	else string_weight_syst = "weight_sfb_down_uncorrelated/weight_sfb";
    }
    else if(systematicName == "syst_l_correlated"){
        if(isUp) string_weight_syst = "weight_sfl_up_correlated/weight_sfb";
	else string_weight_syst = "weight_sfl_down_correlated/weight_sfb";
    }
    else if(systematicName == "syst_l_uncorrelated"){
        if(isUp) string_weight_syst = "weight_sfl_up_uncorrelated/weight_sfb";
	else string_weight_syst = "weight_sfl_down_uncorrelated/weight_sfb";
    }
    else if(systematicName == "syst_prefiring"){
        if(isUp) string_weight_syst = "weight_prefiring_up/weight_prefiring";
	else string_weight_syst = "weight_prefiring_down/weight_prefiring";
    }
    else if(systematicName == "syst_ps_isr"){
        if(isUp) string_weight_syst = "weight_ps_variation4";
	else string_weight_syst = "weight_ps_variation6";
	string_weight_syst = "((" + string_weight_syst + "==-99 || " + string_weight_syst + "==0)?1:"+string_weight_syst + ")";
    }
    else if(systematicName == "syst_ps_fsr"){
        if(isUp) string_weight_syst = "weight_ps_variation5";
        else string_weight_syst = "weight_ps_variation7";
        string_weight_syst = "((" + string_weight_syst + "==-99 || " + string_weight_syst + "==0)?1:"+string_weight_syst + ")";
    }
    else {
        if(isUp) string_weight_syst = "(1+" + systematicName + ")";
	else string_weight_syst = "(1-" + systematicName + ")";
    }

    return string_weight_syst;
}


void Generator::generateLHEweightSystematics(TTree            * tree_p,
                                             std::string const& systematicName,
					     double	      * systList)
{

    std::string leafname;
    if(systematicName == "syst_qcdscale"){
	for (int k=0; k<6; k++) {
	    leafname = "weight_qcdscale_variation" + std::to_string(k);
	    systList[k] = tree_p->GetLeaf(leafname.c_str())->GetValue(0);
	    if (systList[k]==-99 || systList[k]==0) systList[k] = 1;
	}
    }
    else if (systematicName == "syst_pdfas")
    {
        for (int k=0; k<102; k++) {
            leafname = "weight_pdfas_variation_" + std::to_string(k);
            systList[k] = tree_p->GetLeaf(leafname.c_str())->GetValue(0);
	    if (systList[k]==-99 || systList[k]==0) systList[k] = 1;
        }
    }

    return;
}

void Generator::generateLHEweightSystematicsStrings(std::string const& systematicName,
                                                    std::string        * systList)
{
    std::string string_weight_syst = ""; 

    if(systematicName == "syst_qcdscale"){
        for (int k=0; k<6; k++) {
	    string_weight_syst = "weight_qcdscale_variation" + std::to_string(k);
            systList[k] = "((" + string_weight_syst + "==-99 || " + string_weight_syst + "==0)?1:"+string_weight_syst + ")";
	}
    }
    else if (systematicName == "syst_pdfas")
    {
        for (int k=0; k<102; k++) {
            string_weight_syst =  "weight_pdfas_variation_" + std::to_string(k);
	    systList[k] = "((" + string_weight_syst + "==-99 || " + string_weight_syst + "==0)?1:"+string_weight_syst + ")";
	}
    }

    return;
}

void Generator::generateTimeSystematics(std::vector<double>      & weightsUp,
                                        std::vector<double>      & weightsDown
                                        )
{
    std::string timeName = "./inputs/timed/LumiUncertainties_"+year+".root";
    TFile* timeFile = new TFile(timeName.c_str());
    std::vector<TH1F*> timeUp;
    std::vector<TH1F*> timeDown;
    for(auto const& key : *timeFile->GetListOfKeys()){
        if(TString(key->GetName()).Contains("Up"))
            timeUp.push_back((TH1F*)timeFile->Get(key->GetName()));
        if(TString(key->GetName()).Contains("Down"))
            timeDown.push_back((TH1F*)timeFile->Get(key->GetName()));
    }
    double nbin = timeUp[0]->GetNbinsX();
    for(size_t i = 0; i < timeUp.size(); ++i){
        timeUp[i]->Rebin(timeUp[i]->GetNbinsX());
        timeDown[i]->Rebin(timeDown[i]->GetNbinsX());
        if(TString(timeUp[i]->GetName()).Contains("Inclusive")){
            weightsUp[0] = timeUp[i]->Integral()/nbin; 
            weightsDown[0] = timeDown[i]->Integral()/nbin;
        }
        if(TString(timeUp[i]->GetName()).Contains("Stability")){
            weightsUp[1] = timeUp[i]->Integral()/nbin; 
            weightsDown[1] = timeDown[i]->Integral()/nbin;
        }
        if(TString(timeUp[i]->GetName()).Contains("Linearity")){
            weightsUp[2] = timeUp[i]->Integral()/nbin; 
            weightsDown[2] = timeDown[i]->Integral()/nbin;
        }
    }
}

bool Generator::isTriggerPassed(TTree         * tree_p,
                                namelist const& triggerList_p,
                                bool            is2016H
                               )
{
    int trig = 0; 
    if(!is2016H){
        for(size_t i = 0; i < triggerList_p.size(); ++i){
            if(triggerList_p[i] == "trg_muon_electron_mu8ele23DZ_fired" or
               triggerList_p[i] == "trg_muon_electron_mu23ele12DZ_fired")
               continue;
            if(tree_p->GetLeaf(triggerList_p[i].c_str())->GetValue(0) == 1){
                ++trig;
            }
        }
    }
    else{
        for(size_t i = 0; i < triggerList_p.size(); ++i){
            //if(triggerList_p[i] == "trg_muon_electron_mu8ele23_fired" or
            //   triggerList_p[i] == "trg_muon_electron_mu23ele12_fired")
            //   continue;
            if(tree_p->GetLeaf(triggerList_p[i].c_str())->GetValue(0) == 1){
                ++trig;
            }
        }    
    }

    return(trig >= 1);
}

std::string Generator::isTriggerPassedString(namelist const&    triggerList_p,
                                  bool               is2016H
                                 )
{

    std::string string_triggered = "(1";

    if(!is2016H){
        for(size_t i = 0; i < triggerList_p.size(); ++i){
            if(triggerList_p[i] == "trg_muon_electron_mu8ele23DZ_fired" or
               triggerList_p[i] == "trg_muon_electron_mu23ele12DZ_fired")
               continue;
	    //string_triggered +=  "*(" + triggerList_p[i] + "==1)"; 
	    string_triggered +=  "||(" + triggerList_p[i] + "!=0)"; 
	}
    }
    else{
        for(size_t i = 0; i < triggerList_p.size(); ++i){
	    //string_triggered +=  "*(" + triggerList_p[i] + "==1)";
	    string_triggered +=  "||(" + triggerList_p[i] + "[0]==1)";
	}
    }
    
    string_triggered += ")";

    return string_triggered;
}

bool Generator::eventSelection(TTree      	* tree_p) 
{
 
  double n_jets = tree_p->GetLeaf("n_jets")->GetValue(0);
  double n_bjets = tree_p->GetLeaf("n_bjets")->GetValue(0);

  if (n_jets>=2 && n_bjets>=1) return 1;
  else return 0;

}

bool Generator::eventSelection(TTree           * tree_p,
                                std::string     jecName)
{
  std::string nJetCut = "n_jets_"+jecName; 
  std::string nBJetCut = "n_bjets_"+jecName;

  double n_jets = tree_p->GetLeaf(nJetCut.c_str())->GetValue(0);
  double n_bjets = tree_p->GetLeaf(nBJetCut.c_str())->GetValue(0);

  if (n_jets>=2 && n_bjets>=1) return 1;
  else return 0;

}

std::string Generator::eventSelectionString()
{
    std::string string_selection = "n_jets>=2 && n_bjets>=1";
    return string_selection;

}

std::string Generator::eventSelectionString(std::string     jecName)
{
    std::string nJetCut = "n_jets_"+jecName;
    std::string nBJetCut = "n_bjets_"+jecName;

    std::string string_selection = nJetCut+">=2 && " + nBJetCut + ">=1";
    return string_selection;
}

float Generator::getGenObservableValue(TTree         * tree_p)
{
    std::string gen_observable = "Events.gen_"+observable;
    float val = -999;
    if (tree_p->GetLeaf("Events.gen_matched")->GetValue(0)==0) val = -1;
    else if (observable=="m_dilep" || observable=="pt_emu" || observable=="n_bjets") val = tree_p->GetLeaf(gen_observable.c_str())->GetValue(0);
    return val;
}

float Generator::getObservableValue(TTree         * tree_p)
{

    float val = -999;

    if (observable=="m_lblj"){  //mTttbar: demande de connnaitre la masse des jets (ajouter dans les heppy outputs)
        TLorentzVector elec, muon, j1, j2;
        elec.SetPtEtaPhiM(tree_p->GetLeaf("pt_elec")->GetValue(0), tree_p->GetLeaf("eta_elec")->GetValue(0), tree_p->GetLeaf("phi_elec")->GetValue(0),0);
        muon.SetPtEtaPhiM(tree_p->GetLeaf("pt_muon")->GetValue(0), tree_p->GetLeaf("eta_muon")->GetValue(0), tree_p->GetLeaf("phi_muon")->GetValue(0),0);
        j1.SetPtEtaPhiM(tree_p->GetLeaf("j1_pt")->GetValue(0),tree_p->GetLeaf("j1_eta")->GetValue(0),tree_p->GetLeaf("j1_phi")->GetValue(0),tree_p->GetLeaf("j1_m")->GetValue(0));
        j2.SetPtEtaPhiM(tree_p->GetLeaf("j2_pt")->GetValue(0),tree_p->GetLeaf("j2_eta")->GetValue(0),tree_p->GetLeaf("j2_phi")->GetValue(0),tree_p->GetLeaf("j2_m")->GetValue(0));
	val = (elec+muon+j1+j2).M();
    }
    else if (observable=="pt_ttbar"){
        TLorentzVector elec, muon, j1, j2;
        elec.SetPtEtaPhiM(tree_p->GetLeaf("pt_elec")->GetValue(0), tree_p->GetLeaf("eta_elec")->GetValue(0), tree_p->GetLeaf("phi_elec")->GetValue(0),0);
        muon.SetPtEtaPhiM(tree_p->GetLeaf("pt_muon")->GetValue(0), tree_p->GetLeaf("eta_muon")->GetValue(0), tree_p->GetLeaf("phi_muon")->GetValue(0),0);
        j1.SetPtEtaPhiM(tree_p->GetLeaf("j1_pt")->GetValue(0),tree_p->GetLeaf("j1_eta")->GetValue(0),tree_p->GetLeaf("j1_phi")->GetValue(0),tree_p->GetLeaf("j1_m")->GetValue(0));
        j2.SetPtEtaPhiM(tree_p->GetLeaf("j2_pt")->GetValue(0),tree_p->GetLeaf("j2_eta")->GetValue(0),tree_p->GetLeaf("j2_phi")->GetValue(0),tree_p->GetLeaf("j2_m")->GetValue(0));
        val = (elec+muon+j1+j2).Pt();
    }
    //else if (observable=="pt_emu"){
    //    TLorentzVector elec, muon;
    //    elec.SetPtEtaPhiM(tree_p->GetLeaf("pt_elec")->GetValue(0), tree_p->GetLeaf("eta_elec")->GetValue(0), tree_p->GetLeaf("phi_elec")->GetValue(0),0);
    //    muon.SetPtEtaPhiM(tree_p->GetLeaf("pt_muon")->GetValue(0), tree_p->GetLeaf("eta_muon")->GetValue(0), tree_p->GetLeaf("phi_muon")->GetValue(0),0);
    // 	  val = (elec+muon).Pt();
    //}
    else if (observable=="m_lb"){
	TLorentzVector elec, muon, b1, b2, j1, j2;
        elec.SetPtEtaPhiM(tree_p->GetLeaf("pt_elec")->GetValue(0), tree_p->GetLeaf("eta_elec")->GetValue(0), tree_p->GetLeaf("phi_elec")->GetValue(0),0);
        muon.SetPtEtaPhiM(tree_p->GetLeaf("pt_muon")->GetValue(0), tree_p->GetLeaf("eta_muon")->GetValue(0), tree_p->GetLeaf("phi_muon")->GetValue(0),0);
        b1.SetPtEtaPhiM(tree_p->GetLeaf("b1_pt")->GetValue(0),tree_p->GetLeaf("b1_eta")->GetValue(0),tree_p->GetLeaf("b1_phi")->GetValue(0),tree_p->GetLeaf("b1_m")->GetValue(0));
	double b2_pt = tree_p->GetLeaf("b2_pt")->GetValue(0);
	//b2.SetPtEtaPhiM(tree_p->GetLeaf("b2_pt")->GetValue(0),tree_p->GetLeaf("b2_eta")->GetValue(0),tree_p->GetLeaf("b2_phi")->GetValue(0),tree_p->GetLeaf("b2_m")->GetValue(0));
        //j1.SetPtEtaPhiM(tree_p->GetLeaf("j1_pt")->GetValue(0),tree_p->GetLeaf("j1_eta")->GetValue(0),tree_p->GetLeaf("j1_phi")->GetValue(0),tree_p->GetLeaf("j1_m")->GetValue(0));
        //j2.SetPtEtaPhiM(tree_p->GetLeaf("j2_pt")->GetValue(0),tree_p->GetLeaf("j2_eta")->GetValue(0),tree_p->GetLeaf("j2_phi")->GetValue(0),tree_p->GetLeaf("j2_m")->GetValue(0));
	//std::cout << "New event"<<std::endl;
	//std::cout << "b1 pt="<<b1.Pt()<<std::endl;
        //std::cout << "b2 pt="<<b2.Pt()<<std::endl;
        //std::cout << "j1 pt="<<j1.Pt()<<std::endl;
        //std::cout << "j2 pt="<<j2.Pt()<<std::endl;
	TLorentzVector b_top1, lep_top1, b_top2, lep_top2;
	if (b2_pt==99){
            //deltaR(b,lep)  le plus proche : top1  => autre  jet de plus haut pt et l'autre lep : top2
	    b_top1 = b1;
	    if (b1.DeltaR(elec)<b1.DeltaR(muon)) { lep_top1 = elec; lep_top2 = muon; }
	    else { lep_top1 = muon; lep_top2 = elec; }
	    double j1_pt = tree_p->GetLeaf("j1_pt")->GetValue(0);
	    double j2_pt = tree_p->GetLeaf("j2_pt")->GetValue(0);
            //j1.SetPtEtaPhiM(tree_p->GetLeaf("j1_pt")->GetValue(0),tree_p->GetLeaf("j1_eta")->GetValue(0),tree_p->GetLeaf("j1_phi")->GetValue(0),tree_p->GetLeaf("j1_m")->GetValue(0));
            //j2.SetPtEtaPhiM(tree_p->GetLeaf("j2_pt")->GetValue(0),tree_p->GetLeaf("j2_eta")->GetValue(0),tree_p->GetLeaf("j2_phi")->GetValue(0),tree_p->GetLeaf("j2_m")->GetValue(0));
	    if (j1_pt!=b_top1.Pt()) {
                j1.SetPtEtaPhiM(j1_pt,tree_p->GetLeaf("j1_eta")->GetValue(0),tree_p->GetLeaf("j1_phi")->GetValue(0),tree_p->GetLeaf("j1_m")->GetValue(0));
		b_top2 = j1;
	    }
	    else if (j2_pt!=b_top1.Pt()) {
                j2.SetPtEtaPhiM(j2_pt,tree_p->GetLeaf("j2_eta")->GetValue(0),tree_p->GetLeaf("j2_phi")->GetValue(0),tree_p->GetLeaf("j2_m")->GetValue(0));
		b_top2 = j2;
	    }
	}
	else  {
	    b2.SetPtEtaPhiM(b2_pt,tree_p->GetLeaf("b2_eta")->GetValue(0),tree_p->GetLeaf("b2_phi")->GetValue(0),tree_p->GetLeaf("b2_m")->GetValue(0));
            //deltaR(b,lep)  le plus proche : top1  => l'autre b et l'autre lep : top2
	    double dR = b1.DeltaR(elec);
	    b_top1 = b1; lep_top1 = elec; b_top2 = b2; lep_top2 = muon;
	    if (b1.DeltaR(muon) < dR){
		dR = b1.DeltaR(muon);
		b_top1 = b1; b_top2 = b2; 
		lep_top1 = muon; lep_top2 = elec;
	    }
	    if (b2.DeltaR(muon)<dR){
		dR = b2.DeltaR(muon);
		b_top1 = b2; b_top2 = b1;
                lep_top1 = muon; lep_top2 = elec;
	    }
	    if (b2.DeltaR(elec)<dR){
		dR = b2.DeltaR(elec);
                b_top1 = b2; b_top2 = b1;
		lep_top1 = elec; lep_top2 = muon;		
	    }
	}
	val = (b_top2+lep_top2).M();
	//std::cout << "m_lb="<<val<<std::endl;
	    //if (b1.DeltaR(elec)<b1.DeltaR(muon)) val  = (b1+elec).M();
	    //else val = (b1+muon).M();

    }
    else val = tree_p->GetLeaf(observable.c_str())->GetValue(0);

    return val;
}				    

void Generator::write(std::string       const& filename,
                      std::vector<TH1F>      & listObject,
                      std::string       const& option_p
                     )
{
    TFile *output = new TFile(filename.c_str(), option_p.c_str());
    for(size_t i = 0; i < listObject.size(); ++i){
        listObject[i].SetBinContent(nBin, listObject[i].GetBinContent(nBin)+listObject[i].GetBinContent(nBin+1));
        listObject[i].Write();
    }
    output->Close();
}

void Generator::write(std::string       const& filename,
                      std::vector<TH2F>      & listObject,
                      std::string       const& option_p
                     )
{
 
    TFile *output = new TFile(filename.c_str(), option_p.c_str());
    for(size_t i = 0; i < listObject.size(); ++i){
        int nBinX = listObject[i].GetNbinsX();
        int nBinY = listObject[i].GetNbinsY();
	for (int ix=0; ix< nBinX ; ix++){
            listObject[i].SetBinContent(1+ix, nBinY, listObject[i].GetBinContent(1+ix, nBinY)+listObject[i].GetBinContent(1+ix, nBinY+1));
	}
        for (int iy=0; iy< nBinY ; iy++){
            listObject[i].SetBinContent(nBinX, 1+iy, listObject[i].GetBinContent(nBinX, 1+iy)+listObject[i].GetBinContent(nBinX+1, 1+iy));
        }
	listObject[i].SetBinContent(nBinX, nBinY, listObject[i].GetBinContent(nBinX, nBinY)+listObject[i].GetBinContent(nBinX+1, nBinY+1));
        listObject[i].Write();
    }
    output->Close();
}

void Generator::write(std::string       const& filename,
                      std::vector<std::vector<TH1F>> & listObject,
                      std::string       const& option_p
                     )
{
    TFile *output = new TFile(filename.c_str(), option_p.c_str());
    for(size_t l = 0; l < listObject.size(); ++l){
        for(size_t i = 0; i < listObject[l].size(); ++i){        
            listObject[l][i].SetBinContent(nBin, listObject[l][i].GetBinContent(nBin)+listObject[l][i].GetBinContent(nBin+1));
            listObject[l][i].Write();
        }
    }

    output->Close();
}


void Generator::groupingMC(std::vector<TH1F>      & list,
                           namelist          const& groupList_p,
                           bool                     clean
                          )
{

    double nbin = list[0].GetNbinsX();
    double min = list[0].GetXaxis()->GetXmin();
    double max = list[0].GetXaxis()->GetXmax();

    for(std::string group : groupList_p){
        TH1F h((group).c_str(), (group).c_str(), nbin, min, max);
        std::string grp = group.substr(1,group.size()-1);
        for(size_t i = 0; i < list.size(); ++i){
            if(TString(list[i].GetName()).Contains(grp))
                h.Add(&list[i]);
        }
        list.push_back(h);
    }
    if(clean)
        list.erase(list.begin(), list.end()-groupList_p.size());
}

void Generator::groupingMC(std::vector<TH2F>      & list,
                           namelist          const& groupList_p,
                           std::string       const& name,
                           bool                     clean
                          )
{
    
    double nbinX = list[0].GetNbinsX();
    double minX = list[0].GetXaxis()->GetXmin();
    double maxX = list[0].GetXaxis()->GetXmax();
    double nbinY = list[0].GetNbinsY();
    double minY = list[0].GetYaxis()->GetXmin();
    double maxY = list[0].GetYaxis()->GetXmax();

    for(std::string group : groupList_p){
        TH2F h((group+"_"+name).c_str(), (group+"_"+name).c_str(), nbinX, minX, maxX, nbinY, minY, maxY);
        std::string grp = group.substr(1,group.size()-1);
        for(size_t i = 0; i < list.size(); ++i){
            if(TString(list[i].GetName()).Contains(grp))
                h.Add(&list[i]);
        }
	/*
	for (int ix=0; ix<nbinX; ix++){
	    float area = h.Integral(1+ix, 1+ix, 1, nbinY);
	    std::cout << "ix+1="<<ix+1<<" area="<<area<<std::endl;
	    for (int iy=0; iy<nbinY; iy++){
		float bincontent = h.GetBinContent(ix+1,iy+1);
		std::cout << "ix+1="<<ix+1<<" iy+1="<< iy+1<<" old="<<bincontent<<" new="<<bincontent/area <<std::endl;
		h.SetBinContent(ix+1, iy+1, bincontent/area);	
	    }
	}
	*/
        list.push_back(h);
    }
    if(clean)
        list.erase(list.begin(), list.end()-groupList_p.size());
}

void Generator::groupingMC(std::vector<TH1F>      & list,
                           namelist          const& groupList_p,                  
                           std::string       const& name,
                           bool                     clean
                          )
{

    double nbin = list[0].GetNbinsX();
    double min = list[0].GetXaxis()->GetXmin();
    double max = list[0].GetXaxis()->GetXmax();

    for(std::string group : groupList_p){
        TH1F h((group+'_'+name).c_str(), (group+'_'+name).c_str(), nbin, min, max);
        std::string grp = group.substr(1,group.size()-1);
        for(size_t i = 0; i < list.size(); ++i){
            if(TString(list[i].GetName()).Contains(grp))
                h.Add(&list[i]);
        }
        list.push_back(h);
    }
    if(clean)
        list.erase(list.begin(), list.end()-groupList_p.size());
}


void Generator::groupingData(std::vector<TH1F>      & list,
                             namelist          const& groupList_p,                       bool                     clean
                            )
{
    double nbin = list[0].GetNbinsX();
    double min = list[0].GetXaxis()->GetXmin();
    double max = list[0].GetXaxis()->GetXmax();

    for(std::string group : groupList_p){
        TH1F h("data_obs", "data_obs", nbin, min, max);
        for(size_t i = 0; i < list.size(); ++i){
            if(TString(list[i].GetName()).Contains(group))
                h.Add(&list[i]);
        }
        list.push_back(h);
    }
    if(clean)
        list.erase(list.begin(), list.end()-1);
}

void Generator::groupingDataTimed(std::vector<TH1F>      & list,
                                  namelist          const& groupList_p,
                                  int                      bin,            
                                  int                      nBin,           bool                     clean
                                 )
{
    double nbin = list[0].GetNbinsX();
    double min = list[0].GetXaxis()->GetXmin();
    double max = list[0].GetXaxis()->GetXmax();

    for(std::string group : groupList_p){
        TH1F h(("data_obs_bin"+std::to_string(bin)).c_str(), ("data_obs_bin"+std::to_string(bin)).c_str(), nbin, min, max);
        for(size_t i = 0; i < list.size(); ++i){
            if(TString(list[i].GetName()).Contains(group) and 
               TString(list[i].GetName()).Contains("bin"+std::to_string(bin)+'.'))
                h.Add(&list[i]);
        }
        list.push_back(h);
    }
    if(bin == nBin-1){// To clean after all group histogram creation
        if(clean)
            list.erase(list.begin(), list.end()-nBin);
    }
}

void Generator::groupingSystematics(std::vector<TH1F>      & list,
                                    namelist          const& groupList_p,
                                    namelist          const& systematicsList_p,
                                    bool                     isUp,
                                    bool                     clean
                                   )
{
    double nbin = list[0].GetNbinsX();
    double min = list[0].GetXaxis()->GetXmin();
    double max = list[0].GetXaxis()->GetXmax();

    std::string updown;
    if(isUp)
        updown = "Up";
    else
        updown = "Down";

    for(std::string group : groupList_p){
        for(std::string syst : systematicsList_p){
            TH1F h((group+"_"+syst+updown).c_str(), (group+"_"+syst+updown).c_str(), nbin, min, max);  
            std::string grp = group.substr(1,group.size()-1);
            for(size_t i = 0; i < list.size(); ++i){
                if(TString(list[i].GetName()).Contains(grp) and TString(list[i].GetName()).Contains(syst))

                    h.Add(&list[i]);        
            }
            list.push_back(h); 
        }
    }
    if(clean)
        list.erase(list.begin(), list.end()-(groupList_p.size()*systematicsList_p.size()));
}

void Generator::groupingSystematics(std::vector<TH2F>      & list,
                                    namelist          const& groupList_p,
                           	    std::string       const& name,
                                    namelist          const& systematicsList_p,
                                    bool                     isUp,
                                    bool                     clean
                                   )
{
    double nbinX = list[0].GetNbinsX();
    double minX = list[0].GetXaxis()->GetXmin();
    double maxX = list[0].GetXaxis()->GetXmax();
    double nbinY = list[0].GetNbinsY();
    double minY = list[0].GetYaxis()->GetXmin();
    double maxY = list[0].GetYaxis()->GetXmax();

    std::string updown;
    if(isUp)
        updown = "Up";
    else
        updown = "Down";

    for(std::string group : groupList_p){
        for(std::string syst : systematicsList_p){
            TH2F h((group+"_"+name+"_"+syst+updown).c_str(), (group+"_"+name+"_"+syst+updown).c_str(), nbinX, minX, maxX, nbinY, minY, maxY);
            //TH1F h((group+"_"+syst+updown).c_str(), (group+"_"+syst+updown).c_str(), nbin, min, max);
            std::string grp = group.substr(1,group.size()-1);
            for(size_t i = 0; i < list.size(); ++i){
                if(TString(list[i].GetName()).Contains(grp) and TString(list[i].GetName()).Contains(syst))

                    h.Add(&list[i]);
            }
            list.push_back(h);
        }
    }
    if(clean)
        list.erase(list.begin(), list.end()-(groupList_p.size()*systematicsList_p.size()));
}

void Generator::groupingLHEweightSystematics(std::vector<TH1F>      & listLHE,
					     std::vector<TH1F>      & list,
                                    	     namelist          const& groupList_p,
                                    	     std::string           syst,
                                   	     bool                     clean
                                   	    )
{
    double nbin = listLHE[0].GetNbinsX();
    double binmin = listLHE[0].GetXaxis()->GetXmin();
    double binmax = listLHE[0].GetXaxis()->GetXmax();

    int nvar = 0;
    if (syst=="syst_qcdscale") nvar=6;
    if (syst=="syst_pdfas") nvar=102;

    for(std::string group : groupList_p){
        TH1F h_up((group+"_"+syst+"Up").c_str(), (group+"_"+syst+"Up").c_str(), nbin, binmin, binmax);
        TH1F h_down((group+"_"+syst+"Down").c_str(), (group+"_"+syst+"Down").c_str(), nbin, binmin, binmax);
	TH1F h_nom;
        for(size_t i = 0; i < list.size(); ++i){
            if(list[i].GetName()==group) h_nom = list[i];
        }
        std::string grp = group.substr(1,group.size()-1);

	TH1F* h = new TH1F[nvar];
	for (int k=0; k<nvar; k++) h[k] = TH1F((group+"_"+syst+std::to_string(k)).c_str(), (group+"_"+syst+std::to_string(k)).c_str(), nbin, binmin, binmax);
        for(size_t i = 0; i < listLHE.size(); ++i){
	    for (int k=0; k<nvar; k++){
		std::string suffix = std::string(listLHE[i].GetName());
		std::size_t found = suffix.find("pdfas");
		suffix = suffix.substr(found+5, suffix.size()-1);
		//std::cout << "Group "<<group<<" pdf"<<k <<" suffix="<<suffix<<std::endl; 
                if(syst=="syst_qcdscale" && TString(listLHE[i].GetName()).Contains(grp) and TString(listLHE[i].GetName()).Contains(syst+std::to_string(k))) 
		    h[k].Add(&listLHE[i]);

                if(syst=="syst_pdfas" && TString(listLHE[i].GetName()).Contains(grp) and suffix==std::to_string(k)){
		    //std::cout << "Group "<<group<<" pdf"<<k <<" add "<<listLHE[i].GetName()<<std::endl;
                    h[k].Add(&listLHE[i]);
		}
	    }
        }
	for (int j=0; j<=nbin; j++){ //also for the overflow bin
	    double diff = 0;
	    double diff_as = 0;
	    //double diff_min = 0;
	    //double diff_max = 0;
            double min = h_nom.GetBinContent(1+j);
            double max = h_nom.GetBinContent(1+j);
	    //if (syst=="syst_pdfas") std::cout << "bin"<<j<<" nom val="<<h_nom.GetBinContent(1+j) << std::endl;
	    for (int k=0; k<nvar; k++){
	        if (syst=="syst_qcdscale"){ //enveloppe
		    if (h[k].GetBinContent(1+j)<min) min=h[k].GetBinContent(1+j);
                    if (h[k].GetBinContent(1+j)>max) max=h[k].GetBinContent(1+j);
		}
                else if (syst=="syst_pdfas"){ //sum of difference squared
		    if  (k<100){
		      //std::cout << "bin"<<j<<" pdf"<<k<<" val="<<h[k].GetBinContent(1+j) << std::endl;
		      diff += (h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j))*(h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j));
		      //if (h[k].GetBinContent(1+j)<h_nom.GetBinContent(1+j)) diff_min += (h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j))*(h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j));
		      //if (h[k].GetBinContent(1+j)>h_nom.GetBinContent(1+j)) diff_max += (h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j))*(h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j));
		      //std::cout << "bin"<<j<<" pdf"<<k<<" val="<<h[k].GetBinContent(1+j) <<" sqrt(diff)="<<sqrt(diff)<<std::endl;
		    }
		    else if (k==100) diff_as += h[k].GetBinContent(1+j);
		    else if (k==101) diff_as -= h[k].GetBinContent(1+j);
                }
	    }
	    if (syst=="syst_qcdscale"){
	        h_up.SetBinContent(1+j,max);
	        h_down.SetBinContent(1+j,min);
	    }
	    else if (syst=="syst_pdfas"){
		diff_as /= 2.;
		diff_as = diff_as*diff_as;
		h_up.SetBinContent(1+j,h_nom.GetBinContent(1+j)+sqrt(diff+diff_as));
		h_down.SetBinContent(1+j,h_nom.GetBinContent(1+j)-sqrt(diff+diff_as));
		//h_up.SetBinContent(1+j,h_nom.GetBinContent(1+j)+sqrt(diff_max));
		//h_down.SetBinContent(1+j,h_nom.GetBinContent(1+j)-sqrt(diff_min));
	    }    
        }
	//up_integ = h_up.Integral()
	//down_integ = h_down.Integral()
	//h_up.Scale(h_nom.Integral()/up_integ)
	//h_down.Scale(h_nom.Integral()/down_integ)
        listLHE.push_back(h_up);
        listLHE.push_back(h_down);
    }

    if(clean)
        listLHE.erase(listLHE.begin(), listLHE.end()-(groupList_p.size()*2));

}

void Generator::groupingLHEweightSystematics(std::vector<TH2F>      & listLHE,
					     std::vector<TH2F>      & list,
                                    	     namelist          const& groupList_p,
	                                     std::string       const& name,
                                    	     std::string           syst,
                                   	     bool                     clean
                                   	    )
{
    //double nbin = listLHE[0].GetNbinsX();
    //double binmin = listLHE[0].GetXaxis()->GetXmin();
    //double binmax = listLHE[0].GetXaxis()->GetXmax();
    double nbinX = list[0].GetNbinsX();
    double minX = list[0].GetXaxis()->GetXmin();
    double maxX = list[0].GetXaxis()->GetXmax();
    double nbinY = list[0].GetNbinsY();
    double minY = list[0].GetYaxis()->GetXmin();
    double maxY = list[0].GetYaxis()->GetXmax();

    int nvar = 0;
    if (syst=="syst_qcdscale") nvar=6;
    if (syst=="syst_pdfas") nvar=102;

    for(std::string group : groupList_p){
        TH2F h_up((group+"_"+name+"_"+syst+"Up").c_str(), (group+"_"+name+"_"+syst+"Up").c_str(), nbinX, minX, maxX, nbinY, minY, maxY);
        TH2F h_down((group+"_"+name+"_"+syst+"Down").c_str(), (group+"_"+name+"_"+syst+"Down").c_str(), nbinX, minX, maxX, nbinY, minY, maxY);
	TH2F h_nom;
        for(size_t i = 0; i < list.size(); ++i){
            if(list[i].GetName()==group+"_"+name) h_nom = list[i];
        }
        std::string grp = group.substr(1,group.size()-1);

	TH2F* h = new TH2F[nvar];
	for (int k=0; k<nvar; k++) h[k] = TH2F((group+"_"+name+"_"+syst+std::to_string(k)).c_str(), (group+"_"+name+"_"+syst+std::to_string(k)).c_str(), nbinX, minX, maxX, nbinY, minY, maxY);
        for(size_t i = 0; i < listLHE.size(); ++i){
	    for (int k=0; k<nvar; k++){
		std::string suffix = std::string(listLHE[i].GetName());
		std::size_t found = suffix.find("pdfas");
		suffix = suffix.substr(found+5, suffix.size()-1);
		//std::cout << "Group "<<group<<" pdf"<<k <<" suffix="<<suffix<<std::endl; 
                if(syst=="syst_qcdscale" && TString(listLHE[i].GetName()).Contains(grp) and TString(listLHE[i].GetName()).Contains(syst+std::to_string(k))) 
		    h[k].Add(&listLHE[i]);

                if(syst=="syst_pdfas" && TString(listLHE[i].GetName()).Contains(grp) and suffix==std::to_string(k)){
		    //std::cout << "Group "<<group<<" pdf"<<k <<" add "<<listLHE[i].GetName()<<std::endl;
                    h[k].Add(&listLHE[i]);
		}
	    }
        }
	for (int j=0; j<=nbinX; j++){ //also for the overflow bin
	  for (int i=0; i<=nbinY; i++){ 
	    double diff = 0;
	    double diff_as = 0;
	    //double diff_min = 0;
	    //double diff_max = 0;
            double min = h_nom.GetBinContent(1+j, 1+i);
            double max = h_nom.GetBinContent(1+j, 1+i);
	    //if (syst=="syst_pdfas") std::cout << "bin"<<j<<" nom val="<<h_nom.GetBinContent(1+j) << std::endl;
	    for (int k=0; k<nvar; k++){
	        if (syst=="syst_qcdscale"){ //enveloppe
		    if (h[k].GetBinContent(1+j,1+i)<min) min=h[k].GetBinContent(1+j, 1+i);
                    if (h[k].GetBinContent(1+j,1+i)>max) max=h[k].GetBinContent(1+j, 1+i);
		}
                else if (syst=="syst_pdfas"){ //sum of difference squared
		    if  (k<100){
		      //std::cout << "bin"<<j<<" pdf"<<k<<" val="<<h[k].GetBinContent(1+j) << std::endl;
		      diff += (h[k].GetBinContent(1+j, 1+i)-h_nom.GetBinContent(1+j, 1+i))*(h[k].GetBinContent(1+j,1+i)-h_nom.GetBinContent(1+j,1+i));
		      //if (h[k].GetBinContent(1+j)<h_nom.GetBinContent(1+j)) diff_min += (h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j))*(h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j));
		      //if (h[k].GetBinContent(1+j)>h_nom.GetBinContent(1+j)) diff_max += (h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j))*(h[k].GetBinContent(1+j)-h_nom.GetBinContent(1+j));
		      //std::cout << "bin"<<j<<" pdf"<<k<<" val="<<h[k].GetBinContent(1+j) <<" sqrt(diff)="<<sqrt(diff)<<std::endl;
		    }
		    else if (k==100) diff_as += h[k].GetBinContent(1+j,1+i);
		    else if (k==101) diff_as -= h[k].GetBinContent(1+j,1+i);
                }
	    }
	    if (syst=="syst_qcdscale"){
	        h_up.SetBinContent(1+j,1+i,max);
	        h_down.SetBinContent(1+j,1+i,min);
	    }
	    else if (syst=="syst_pdfas"){
		diff_as /= 2.;
		diff_as = diff_as*diff_as;
		h_up.SetBinContent(1+j,1+i,h_nom.GetBinContent(1+j,1+i)+sqrt(diff+diff_as));
		h_down.SetBinContent(1+j,1+i,h_nom.GetBinContent(1+j,1+i)-sqrt(diff+diff_as));
		//h_up.SetBinContent(1+j,h_nom.GetBinContent(1+j)+sqrt(diff_max));
		//h_down.SetBinContent(1+j,h_nom.GetBinContent(1+j)-sqrt(diff_min));
	    }
	  }    
        }
	//up_integ = h_up.Integral()
	//down_integ = h_down.Integral()
	//h_up.Scale(h_nom.Integral()/up_integ)
	//h_down.Scale(h_nom.Integral()/down_integ)
        listLHE.push_back(h_up);
        listLHE.push_back(h_down);
    }

    if(clean)
        listLHE.erase(listLHE.begin(), listLHE.end()-(groupList_p.size()*2));

}


/////////////////////////////////
// Public methods
/////////////////////////////////


void Generator::generateJecMC(namelist            const& sampleList_p,
                               namelist            const& jecList_p,
                               namelist            const& groupList_p,
                               namelist            const& triggerList_p,
                               //std::vector<std::vector<double>> const& correction_p,
                               std::vector<double> const& correction_p,
                               bool                clean_p,
                               bool                isTimed_p,
                               int                 timebin
                    )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<std::vector<TH1F> > list(jecList_p.size());
    std::vector<std::vector<TH2F> > listResponseMatrix(jecList_p.size());

    std::string cleaned="";
    if(!clean_p) cleaned = "_unclean";

    std::string sinc="";
    if (!isTimed_p) sinc = "_inclusive";

    std::string stimebin="";
    if (timebin==-1) stimebin = "_puold";
    if (timebin==-2) stimebin = "_punew";
    if (timebin==-3) stimebin = "_puinc";
    if (timebin>=0) stimebin = "_put"+std::to_string(timebin);

    std::string filename_p = "./results/"+year+"/flattree/"+observable+"_jec"+sinc+cleaned+stimebin+".root";
    
    for(size_t jl = 0; jl < jecList_p.size(); ++jl)
    {
        for(size_t n = 0; n < sampleList_p.size(); ++n){
            
            //td::cout << " -> " << jecList_p[jl] << " : " << sampleList_p[n] << "  " << correction_p[jl][n] << std::endl;
            std::cout << " -> " << jecList_p[jl] << " : " << sampleList_p[n] << "  " << correction_p[n] << std::endl;
            
            //std::string filename = "./inputs/"+year+"/JEC/"+jecList_p[jl]+'/'+sampleList_p[n]+"/NtupleProducer/tree.root";
	    std::string filename = "./inputs/"+year+"/MC/"+sampleList_p[n]+"/NtupleProducer/tree.root";

            bool doResponseMatrix = false;
            if ((TString(sampleList_p[n]).Contains("signal") || TString(sampleList_p[n]).Contains("singletop")) && isTimed_p) doResponseMatrix = true;

            TFile* file = new TFile(filename.c_str());
            TTree *tree;
            file->GetObject("events", tree);
            TCanvas *canvas = new TCanvas((jecList_p[jl]+sampleList_p[n]).c_str());
            TH1F* hist      = new TH1F((jecList_p[jl]+sampleList_p[n]).c_str(), (jecList_p[jl]+sampleList_p[n]).c_str(), nBin, minBin, maxBin);

	    TFile* fileGEN;
	    TTree* treeGEN;
	    TH2F* hist_responseMatrix;
	    std::string filename_gen;
	    if (doResponseMatrix){
		std::string genSample = sampleList_p[n];
		genSample = genSample.substr(3, genSample.size()-1);
		filename_gen = "./inputs/"+year+"/GEN/"+genSample+"_NanoGEN_"+year+"_selection_particle.root";
		std::cout <<"Gen events from "<<filename_gen<<std::endl;
		fileGEN = new TFile(filename_gen.c_str(),"READ");
		treeGEN = (TTree*)fileGEN->Get("Events");
		tree->AddFriend(treeGEN);
	    }
            if (doResponseMatrix){
                if (observable=="n_bjets") hist_responseMatrix = new TH2F((jecList_p[jl]+sampleList_p[n] + "_responseMatrix").c_str(), observable.c_str(), nBin, minBin, maxBin, nBin+2, -1, maxBin);
                else hist_responseMatrix = new TH2F((jecList_p[jl]+sampleList_p[n] + "_responseMatrix").c_str(), observable.c_str(), nBin+1, minBin-1, maxBin, nBin, minBin, maxBin);
	    }

	    if (doLoop){
		//int nEvents = 10;
		int nEvents = tree->GetEntriesFast();
		for(int i = 0; i < nEvents; ++i){
		    tree->GetEntry(i);
		    if (eventSelection(tree, jecList_p[jl])!=1) continue;
		    double weight = generateWeight(tree,isTimed_p);
		    if(isTriggerPassed(tree, triggerList_p,true)){
			    hist->Fill(getObservableValue(tree), weight);
			    //hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight);
		    }
		    if(i % 100000 == 0)
			std::cout << "100 000 events passed" << std::endl;
		}
	    }
            else if (!doLoop){
                std::string string_eventSelection = eventSelectionString(jecList_p[jl]);
                std::string string_weight = generateWeightString(isTimed_p,timebin);
                std::string string_triggered = isTriggerPassedString(triggerList_p,true);
                drawHisto1D(tree, observable, string_eventSelection, string_weight, string_triggered, hist);
                if (doResponseMatrix)
                    drawHisto2D(tree, observable, "Events.gen_"+observable, string_eventSelection, string_weight+"*(Events.gen_matched!=0)", string_triggered, hist_responseMatrix);
	    }

            //hist->Scale(correction_p[jl][n]);
            hist->Scale(correction_p[n]);
            list[jl].push_back(*hist);
            if (doResponseMatrix){
                hist_responseMatrix->Scale(correction_p[n]);
                listResponseMatrix[jl].push_back(*hist_responseMatrix);
	    }

            delete hist;
            delete canvas;
            delete tree;
            delete file;
        }
        groupingMC(list[jl], groupList_p, jecList_p[jl], clean_p);

        std::vector<std::string> groupList_responseMatrix;
	if(isTimed_p){
        groupList_responseMatrix.push_back("signal");
        groupList_responseMatrix.push_back("singletop");
        groupingMC(listResponseMatrix[jl], groupList_responseMatrix, "responseMatrix_"+jecList_p[jl], clean_p);
	}
    }
    write(filename_p, list, "RECREATE");
    if(isTimed_p){
    for(size_t jl = 0; jl < jecList_p.size(); ++jl) write(filename_p, listResponseMatrix[jl], "UPDATE");
    }
}


void Generator::generateAltMC(namelist            const& sampleList_p,
                           namelist            const& groupList_p,
                    namelist            const& triggerList_p,
                    std::vector<double> const& correction_p,
                    bool                       clean_p,
                    bool                       isTimed_p,
		    int   	               timebin
                    )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;

    std::string cleaned="";
    if(!clean_p) cleaned = "_unclean";

    std::string sinc="";
    if (!isTimed_p) sinc = "_inclusive";

    std::string stimebin="";
    if (timebin==-1) stimebin = "_puold";
    if (timebin==-2) stimebin = "_punew";
    if (timebin==-3) stimebin = "_puinc";
    if (timebin>=0) stimebin = "_put"+std::to_string(timebin);

    std::string filename_p = "./results/"+year+"/flattree/"+observable+"_alt"+sinc+cleaned+stimebin+".root";
    

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/ALT/"+sampleList_p[n]+"/NtupleProducer/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);

        std::cout << "ALT -> " << sampleList_p[n] << "  " << correction_p[n] << std::endl;
	
	if (doLoop){
	    for(int i = 0; i < tree->GetEntriesFast(); ++i){
		tree->GetEntry(i);
		if (eventSelection(tree)!=1) continue;
		double weight = generateWeight(tree,isTimed_p);
		//if(sampleList_p[n].find("alt") == std::string::npos){
		    if(isTriggerPassed(tree, triggerList_p,true)){
			hist->Fill(getObservableValue(tree), weight);
			//hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight);
		    }
		//}
		//else{
		//    hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight);
		//}
		
		if(i % 100000 == 0)
		    std::cout << "100 000 events passed" << std::endl;
	    }
	}
	else if (!doLoop){
            std::string string_eventSelection = eventSelectionString();
            std::string string_weight = generateWeightString(isTimed_p, timebin);
            std::string string_triggered = isTriggerPassedString(triggerList_p,true);
            drawHisto1D(tree, observable, string_eventSelection, string_weight, string_triggered, hist);
	}
        hist->Scale(correction_p[n]);
        list.push_back(*hist);

        delete hist;
        delete canvas;
        delete tree;
        delete file;
    }
    groupingMC(list, groupList_p, clean_p);

    for (unsigned int i=0; i<list.size(); i++){
      std::cout << list.at(i).GetName() << std::endl;
      std::string histname = "signal_" + std::string(list.at(i).GetName());
      int pos = histname.find("169");
      if (pos!=-1) histname.replace(pos, 3, "Down");
      pos = histname.find("175");
      if (pos!=-1) histname.replace(pos, 3, "Up");
      list.at(i).SetName(histname.c_str());
    }

    write(filename_p, list, "RECREATE");
}


void Generator::generateMC(namelist            const& sampleList_p,
                           namelist            const& triggerList_p,
                           namelist            const& groupList_p,
                           namelist            const& systematicsList_p,
                           namelist            const& systematicsTimeList_p,
			   std::vector<double>    const& numberofevents_p,
                           std::vector<double> const& correction_p,
                           std::string         const& option_p,
                           bool                       clean_p,
                           bool                       isTimed_p,
	                   int                        timebin
                          )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::vector<TH1F> listUp;    
    std::vector<TH1F> listDown;
    //std::vector<TH1F> listTimeUp;    
    //std::vector<TH1F> listTimeDown;
    std::vector<TH1F> listQCDscale;
    std::vector<TH1F> listPDFas;
    std::vector<TH2F> listResponseMatrix;
    std::vector<TH2F> listResponseMatrixUp;
    std::vector<TH2F> listResponseMatrixDown;
    std::vector<TH2F> listResponseMatrixQCDscale;
    std::vector<TH2F> listResponseMatrixPDFas;

    std::string cleaned="";
    if(!clean_p) cleaned = "_unclean";

    std::string sinc=""; 
    if (!isTimed_p) sinc = "_inclusive";

    std::string stimebin="";
    if (timebin==-1) stimebin = "_puold";
    if (timebin==-2) stimebin = "_punew";
    if (timebin==-3) stimebin = "_puinc";
    if (timebin>=0) stimebin = "_put"+std::to_string(timebin);

    std::string filename_p = "./results/"+year+"/flattree/"+observable+sinc+cleaned+stimebin+".root";
    //std::string filenameTime_p = "./results/"+year+"/flattree/"+observable+"_timed"+cleaned+".root";

    //std::vector<double> timeSystWeightUp(systematicsTimeList_p.size(), 0);
    //std::vector<double> timeSystWeightDown(systematicsTimeList_p.size(), 0);
    //generateTimeSystematics(timeSystWeightUp, timeSystWeightDown);

    for(size_t i = 0; i < systematicsList_p.size(); ++i)
            std::cout  <<  systematicsList_p[i] << std::endl;

    bool doResponseMatrix = false;
    for(size_t n = 0; n < sampleList_p.size(); ++n){
      doResponseMatrix = false;
    
      if ((TString(sampleList_p[n]).Contains("signal") || TString(sampleList_p[n]).Contains("singletop")) && isTimed_p) doResponseMatrix = true;
	
        std::string filename = "./inputs/"+year+"/MC/"+sampleList_p[n]+"/NtupleProducer/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);

        TFile* fileGEN;
        TTree* treeGEN;
        std::string filename_gen;
	if (doResponseMatrix){
            std::string genSample = sampleList_p[n];
	    genSample = genSample.substr(3, genSample.size()-1);
            filename_gen = "./inputs/"+year+"/GEN/"+genSample+"_NanoGEN_"+year+"_selection_particle.root";
	    std::cout <<"Gen events from "<<filename_gen<<std::endl;
	    fileGEN = new TFile(filename_gen.c_str(),"READ");
	    treeGEN = (TTree*)fileGEN->Get("Events");
	    tree->AddFriend(treeGEN);
        }

        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);
        TH1F* hist_events      = new TH1F((sampleList_p[n] + "_events").c_str(), observable.c_str(), nBin, minBin, maxBin);
	TH2F* hist_responseMatrix;
	if (doResponseMatrix){
	    if (observable=="n_bjets") hist_responseMatrix = new TH2F((sampleList_p[n] + "_responseMatrix").c_str(), observable.c_str(), nBin, minBin, maxBin, nBin+2, -1, maxBin);
	    else hist_responseMatrix = new TH2F((sampleList_p[n] + "_responseMatrix").c_str(), observable.c_str(), nBin+1, minBin-1, maxBin, nBin, minBin, maxBin);
	}

        std::vector<TH1F*> histUp(systematicsList_p.size()-2);
        std::vector<TH1F*> histDown(systematicsList_p.size()-2);
        //std::vector<TH1F*> histUpTime(systematicsList_p.size());
        //std::vector<TH1F*> histDownTime(systematicsList_p.size());
        std::vector<TH1F*> histVarQCDscale(6);
        std::vector<TH1F*> histVarPDFas(102);
	
	std::vector<TH2F*> hist_responseMatrixUp(systematicsList_p.size()-2);
	std::vector<TH2F*> hist_responseMatrixDown(systematicsList_p.size()-2);
	std::vector<TH2F*> hist_responseMatrixVarQCDscale(6);
	std::vector<TH2F*> hist_responseMatrixVarPDFas(102);
        for(size_t i = 0; i < systematicsList_p.size(); ++i){
	    //std::cout  <<  systematicsList_p[i] << std::endl;
            if (systematicsList_p[i]!="syst_qcdscale" && systematicsList_p[i]!="syst_pdfas"){
               histUp[i]   = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+"Up").c_str(), (observable+"Up").c_str(), nBin, minBin, maxBin);
               histDown[i] = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+"Down").c_str(), (observable+"Down").c_str(), nBin, minBin, maxBin);
	       if (doResponseMatrix){
		   hist_responseMatrixUp[i] = new TH2F((sampleList_p[n] + "_responseMatrix_"+systematicsList_p[i]+"Up").c_str(), (observable+"Up").c_str(), nBin, minBin, maxBin, nBin+2, -1, maxBin);
                   hist_responseMatrixDown[i] = new TH2F((sampleList_p[n] + "_responseMatrix_"+systematicsList_p[i]+"Down").c_str(), (observable+"Down").c_str(), nBin, minBin, maxBin, nBin+2, -1, maxBin);
	       } 
	    }
            else if (systematicsList_p[i]=="syst_qcdscale"){
	       for(size_t k = 0; k < 6; ++k) {
		   histVarQCDscale[k] = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+std::to_string(k)).c_str(), (observable+std::to_string(k)).c_str(), nBin, minBin, maxBin);
		   if (doResponseMatrix)
			hist_responseMatrixVarQCDscale[k] = new TH2F((sampleList_p[n]+"_responseMatrix_"+systematicsList_p[i]+std::to_string(k)).c_str(), (observable+std::to_string(k)).c_str(), nBin, minBin, maxBin, nBin+2, -1, maxBin); 
	       }
	    }
            else if (systematicsList_p[i]=="syst_pdfas"){
	       for (size_t k = 0; k < 102; ++k) {
		   histVarPDFas[k] = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+std::to_string(k)).c_str(), (observable+std::to_string(k)).c_str(), nBin, minBin, maxBin);
		   if (doResponseMatrix)
			hist_responseMatrixVarPDFas[k] = new TH2F((sampleList_p[n]+"_responseMatrix_"+systematicsList_p[i]+std::to_string(k)).c_str(), (observable+std::to_string(k)).c_str(), nBin, minBin, maxBin, nBin+2, -1, maxBin);

	       }
	    }
        }
	
        //for(size_t i = 0; i < systematicsTimeList_p.size(); ++i){
        //    histUpTime[i]   = new TH1F((sampleList_p[n]+"_"+systematicsTimeList_p[i]+"Up").c_str(), (observable+"Up").c_str(), nBin, minBin, maxBin);
        //    histDownTime[i] = new TH1F((sampleList_p[n]+"_"+systematicsTimeList_p[i]+"Down").c_str(), (observable+"Down").c_str(), nBin, minBin, maxBin);
        //}

        std::cout << " -> " << sampleList_p[n] << "  " << correction_p[n] << std::endl;

	double systUp = 0;
	double systDown = 0;
	double* systVarQCDscale = new double[6];
	double* systVarPDFas = new double[102];
        std::string* strings_systVarQCDscale = new std::string[6];
        std::string* strings_systVarPDFas = new std::string[102];

	//int nEvents = tree->GetEntriesFast();
	int nEvents = 1000;

        if (doLoop){
	    for(int i = 0; i < nEvents; ++i){
		tree->GetEntry(i);
		if (eventSelection(tree)!=1) continue;
		double weight = generateWeight(tree,isTimed_p);
		if(isTriggerPassed(tree, triggerList_p,true)){
		    hist_events->Fill(getObservableValue(tree));
		    hist->Fill(getObservableValue(tree), weight);
		    if (doResponseMatrix) 
			hist_responseMatrix->Fill(getObservableValue(tree), getGenObservableValue(tree), weight);
		    for(size_t j = 0; j < systematicsList_p.size(); ++j){
			if (systematicsList_p[j]!="syst_qcdscale" && systematicsList_p[j]!="syst_pdfas"){
			    systUp = generateSystematics(tree, systematicsList_p[j], true);
			    systDown = generateSystematics(tree, systematicsList_p[j], false);
			    histUp[j]->Fill(getObservableValue(tree), weight*systUp);
			    histDown[j]->Fill(getObservableValue(tree), weight*systDown);
			    if (doResponseMatrix){
				hist_responseMatrixUp[j]->Fill(getObservableValue(tree), getGenObservableValue(tree), weight*systUp);
				hist_responseMatrixDown[j]->Fill(getObservableValue(tree), getGenObservableValue(tree), weight*systDown);
			    }
			}   
			else if (systematicsList_p[j]=="syst_qcdscale"){
			    generateLHEweightSystematics(tree, systematicsList_p[j], systVarQCDscale);
			    for (int k=0; k<6; k++) {
			       histVarQCDscale[k]->Fill(getObservableValue(tree), weight*systVarQCDscale[k]);
			       if (doResponseMatrix) 
				  hist_responseMatrixVarQCDscale[k]->Fill(getObservableValue(tree), getGenObservableValue(tree), weight*systVarQCDscale[k]);
			    }
			}                     
			else if (systematicsList_p[j]=="syst_pdfas"){
			    generateLHEweightSystematics(tree, systematicsList_p[j], systVarPDFas);
			    for (int k=0; k<102; k++) {
			       histVarPDFas[k]->Fill(getObservableValue(tree), weight*systVarPDFas[k]);
			       if (doResponseMatrix) 
				   hist_responseMatrixVarPDFas[k]->Fill(getObservableValue(tree), getGenObservableValue(tree), weight*systVarPDFas[k]);
			    }
			}
		    }
		    //for(size_t j = 0; j < systematicsTimeList_p.size(); ++j){
		    //    double systUp = timeSystWeightUp[j];
		    //    double systDown = timeSystWeightDown[j];
		    //    histUpTime[j]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight*systUp);
		    //    histDownTime[j]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight*systDown);                        
		    //}
		}
		if(i % 100000 == 0)
		    std::cout << "Event "<<i<<" / "<<nEvents << std::endl;
	    }
	}
	else if (!doLoop){
            std::string string_eventSelection = eventSelectionString();
            std::string string_weight = generateWeightString(isTimed_p, timebin);
            std::string string_triggered = isTriggerPassedString(triggerList_p,true);
            //Nominal hist
            drawHisto1D(tree, observable, string_eventSelection, string_weight, string_triggered, hist);
	    //Hist without weight
	    drawHisto1D(tree, observable, string_eventSelection, "1", string_triggered, hist_events);
	    //Nominal response matrix
	    if (doResponseMatrix)
		drawHisto2D(tree, observable, "Events.gen_"+observable, string_eventSelection, string_weight+"*(Events.gen_matched!=0)", string_triggered, hist_responseMatrix);
	    //Systematics
	    for(size_t j = 0; j < systematicsList_p.size(); ++j){
                if (systematicsList_p[j]!="syst_qcdscale" && systematicsList_p[j]!="syst_pdfas"){
		    std::string string_syst_up = generateSystematicsString(systematicsList_p[j], true, timebin);
                    std::string string_syst_down = generateSystematicsString(systematicsList_p[j], false, timebin);
	            drawHisto1D(tree, observable, string_eventSelection, string_weight+"*"+string_syst_up, string_triggered, histUp[j]);
                    drawHisto1D(tree, observable, string_eventSelection, string_weight+"*"+string_syst_down, string_triggered, histDown[j]);
		    if (doResponseMatrix){
			drawHisto2D(tree, observable, "Events.gen_"+observable, string_eventSelection, string_weight+"*(Events.gen_matched!=0)*"+string_syst_up, string_triggered, hist_responseMatrixUp[j]);
                        drawHisto2D(tree, observable, "Events.gen_"+observable, string_eventSelection, string_weight+"*(Events.gen_matched!=0)*"+string_syst_down, string_triggered, hist_responseMatrixDown[j]);
		    }
		}
		else if (systematicsList_p[j]=="syst_qcdscale"){
		    generateLHEweightSystematicsStrings(systematicsList_p[j], strings_systVarQCDscale);
		    for (int k=0; k<6; k++) {
			drawHisto1D(tree, observable, string_eventSelection, string_weight+"*"+strings_systVarQCDscale[k], string_triggered, histVarQCDscale[k]);
			if (doResponseMatrix)
			    drawHisto2D(tree, observable, "Events.gen_"+observable, string_eventSelection, string_weight+"*(Events.gen_matched!=0)*"+strings_systVarQCDscale[k], string_triggered, hist_responseMatrixVarQCDscale[k]);
		    }
		}
		else if (systematicsList_p[j]=="syst_pdfas"){
                    generateLHEweightSystematicsStrings(systematicsList_p[j], strings_systVarPDFas);
                    for (int k=0; k<102; k++) {
                        drawHisto1D(tree, observable, string_eventSelection, string_weight+"*"+strings_systVarPDFas[k], string_triggered, histVarPDFas[k]);
                        if (doResponseMatrix)
                            drawHisto2D(tree, observable, "Events.gen_"+observable, string_eventSelection, string_weight+"*(Events.gen_matched!=0)*"+strings_systVarPDFas[k], string_triggered, hist_responseMatrixVarPDFas[k]);
                    }
		}
	    }
	}
        hist->Scale(correction_p[n]);
	if (doResponseMatrix) 
	    hist_responseMatrix->Scale(correction_p[n]);


        bool doRecomputeUncert = false;

	double bin_uncertainty_weight = 0;
	double bin_uncertainty_binomial = 0;
	for (int iobs=1; iobs<=nBin; iobs++){
	    if (doRecomputeUncert) {
		bin_uncertainty_weight = (hist->GetBinError(iobs))*(hist->GetBinError(iobs));
		if (hist_events->GetBinContent(iobs)>0) bin_uncertainty_binomial = hist->GetBinContent(iobs)/((double)hist_events->GetBinContent(iobs)) * sqrt(((double)hist_events->GetBinContent(iobs)) * (1 - ((double)hist_events->GetBinContent(iobs))/((double)numberofevents_p[n])));
		hist->SetBinError(iobs, sqrt(bin_uncertainty_weight*bin_uncertainty_weight+bin_uncertainty_binomial*bin_uncertainty_binomial));
	    }
	}

        list.push_back(*hist);
	if (doResponseMatrix)
	    listResponseMatrix.push_back(*hist_responseMatrix);
        for(size_t i = 0; i < systematicsList_p.size(); ++i){
	    if (systematicsList_p[i]!="syst_qcdscale" && systematicsList_p[i]!="syst_pdfas"){
                histUp[i]->Scale(correction_p[n]);
                histDown[i]->Scale(correction_p[n]);
                listUp.push_back(*histUp[i]);
                listDown.push_back(*histDown[i]);
		if (doResponseMatrix){
		    hist_responseMatrixUp[i]->Scale(correction_p[n]);
		    hist_responseMatrixDown[i]->Scale(correction_p[n]);
		    listResponseMatrixUp.push_back(*hist_responseMatrixUp[i]);
                    listResponseMatrixDown.push_back(*hist_responseMatrixDown[i]);
		}
	    }
	    else if (systematicsList_p[i]=="syst_qcdscale"){
		for (int k=0; k<6; k++) {
		    histVarQCDscale[k]->Scale(correction_p[n]);
		    listQCDscale.push_back(*histVarQCDscale[k]);
		    if (doResponseMatrix){
			hist_responseMatrixVarQCDscale[k]->Scale(correction_p[n]);
			listResponseMatrixQCDscale.push_back(*hist_responseMatrixVarQCDscale[k]);
		    }
		}
	    }
	    else if (systematicsList_p[i]=="syst_pdfas"){
                for (int k=0; k<102; k++) { 
		    histVarPDFas[k]->Scale(correction_p[n]);
		    listPDFas.push_back(*histVarPDFas[k]);
   		    if (doResponseMatrix){
			hist_responseMatrixVarPDFas[k]->Scale(correction_p[n]);
			listResponseMatrixPDFas.push_back(*hist_responseMatrixVarPDFas[k]);
		      }
		}
	    }
	}
        //for(size_t i = 0; i < systematicsTimeList_p.size(); ++i){
            //histUpTime[i]->Scale(correction_p[n]);
            //histDownTime[i]->Scale(correction_p[n]);
            //listTimeUp.push_back(*histUpTime[i]);
            //listTimeDown.push_back(*histDownTime[i]);
        //}


        delete hist;
	if (doResponseMatrix) delete hist_responseMatrix;
        for(size_t i = 0; i < systematicsList_p.size()-2; ++i){
            delete histUp[i];
            delete histDown[i];
	    if (doResponseMatrix){
	      delete hist_responseMatrixUp[i];
	      delete hist_responseMatrixDown[i];
	    }
	  //delete histUpTime[i];
	  //delete histDownTime[i];
        }
	for (int k=0; k<6; k++) {
	    delete histVarQCDscale[k];
	    if (doResponseMatrix)
	        delete hist_responseMatrixVarQCDscale[k];
	}
	for (int k=0; k<102; k++) {
	    delete histVarPDFas[k];
            if (doResponseMatrix)
		delete hist_responseMatrixVarPDFas[k];
	}
        delete canvas;
        delete tree;
        delete file;
    }

 
    groupingMC(list, groupList_p, clean_p);
    groupingSystematics(listUp, groupList_p, systematicsList_p, true, clean_p);    // isUp = true
    groupingSystematics(listDown, groupList_p, systematicsList_p, false, clean_p); // isUp = false
    
    std::vector<std::string> groupList_responseMatrix;
    if(isTimed_p){
      groupList_responseMatrix.push_back("signal");
      groupList_responseMatrix.push_back("singletop");
      groupingMC(listResponseMatrix, groupList_responseMatrix, "responseMatrix", clean_p);
      groupingSystematics(listResponseMatrixUp, groupList_responseMatrix, "responseMatrix", systematicsList_p, true, clean_p); // isUp = true
      groupingSystematics(listResponseMatrixDown, groupList_responseMatrix, "responseMatrix", systematicsList_p, false, clean_p); // isUp = false
      }
    
    unsigned int it=0;
    while (it < listUp.size()){
	if (TString(listUp[it].GetName()).Contains("syst_qcdscale") || TString(listUp[it].GetName()).Contains("syst_pdfas")) {
	    listUp.erase(listUp.begin()+it);
	}
	else it++;
    }
    it=0;
    while (it < listDown.size()){
        if (TString(listDown[it].GetName()).Contains("syst_qcdscale") || TString(listDown[it].GetName()).Contains("syst_pdfas")) {
	    listDown.erase(listDown.begin()+it);
	}
	else it++;
    }
    
    if(isTimed_p){
      it=0;
      while (it < listResponseMatrixUp.size()){
        if (TString(listResponseMatrixUp[it].GetName()).Contains("syst_qcdscale") || TString(listResponseMatrixUp[it].GetName()).Contains("syst_pdfas")) {
	  listResponseMatrixUp.erase(listResponseMatrixUp.begin()+it);
        }
        else it++;
      }
    

      it=0;
      while (it < listResponseMatrixDown.size()){
        if (TString(listResponseMatrixDown[it].GetName()).Contains("syst_qcdscale") || TString(listResponseMatrixDown[it].GetName()).Contains("syst_pdfas")) {
	  listResponseMatrixDown.erase(listResponseMatrixDown.begin()+it);
        }
        else it++;
      }
      }
    
    groupingLHEweightSystematics(listQCDscale, list, groupList_p, "syst_qcdscale", clean_p);
    groupingLHEweightSystematics(listPDFas, list, groupList_p, "syst_pdfas", clean_p);
    //groupingSystematics(listTimeUp, groupList_p, systematicsTimeList_p, true, clean_p);    // isUp = true
    //groupingSystematics(listTimeDown, groupList_p, systematicsTimeList_p, false, clean_p); // isUp = false
    if(isTimed_p){
      groupingLHEweightSystematics(listResponseMatrixQCDscale, listResponseMatrix, groupList_responseMatrix, "responseMatrix", "syst_qcdscale", clean_p);
      groupingLHEweightSystematics(listResponseMatrixPDFas, listResponseMatrix, groupList_responseMatrix, "responseMatrix", "syst_pdfas", clean_p);
       }
    write(filename_p, list, option_p);
    write(filename_p, listUp, "UPDATE");
    write(filename_p, listDown, "UPDATE");
    write(filename_p, listQCDscale, "UPDATE");
    write(filename_p, listPDFas, "UPDATE");
    //write(filenameTime_p, listTimeUp, "RECREATE");
    //write(filenameTime_p, listTimeDown, "UPDATE"); 
    if(isTimed_p){
      write(filename_p, listResponseMatrix, "UPDATE");
      write(filename_p, listResponseMatrixUp, "UPDATE");
      write(filename_p, listResponseMatrixDown, "UPDATE");
      write(filename_p, listResponseMatrixQCDscale, "UPDATE");
      write(filename_p, listResponseMatrixPDFas, "UPDATE");
      }
}


void Generator::generateMCforComp(namelist            const& sampleList_p,
                           namelist            const& triggerList_p,
                           namelist            const& groupList_p,
                           std::vector<double> const& correction_p,
                           std::string         const& option_p,
                           bool                       clean_p,
			   int			      timebin
                          )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    //std::vector<TH1D> list2;

    std::string cleaned;
    if(!clean_p) cleaned = "_unclean";

    std::string stimebin="";
    if (timebin==-1) stimebin = "_puold";
    if (timebin==-2) stimebin = "_punew";
    if (timebin==-3) stimebin = "_puinc";

    std::string filename_p = "./results/"+year+"/flattree/"+observable+cleaned+"_forComp"+stimebin+".root";

    //ROOT::EnableImplicitMT();

    for(size_t n = 0; n < sampleList_p.size(); ++n){
    
        std::string filename = "./inputs/"+year+"/MC/"+sampleList_p[n]+"/NtupleProducer/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        //TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
	TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);
        //TH1D* hist      = new TH1D(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);

        std::cout << " -> " << sampleList_p[n] << "  " << correction_p[n] << std::endl;

	if (doLoop){
            TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);
            for(int i = 0; i < tree->GetEntriesFast(); ++i){
                tree->GetEntry(i);
                if (eventSelection(tree)!=1) continue;
                double weight = generateWeight(tree, false);
                if(isTriggerPassed(tree, triggerList_p,true)){
		  //hist->Fill(getObservableValue(tree), weight);
		  hist->Fill(getObservableValue(tree));
                    //hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight);
                }
                if(i % 100000 == 0)
                    std::cout << "100 000 events passed" << std::endl;
            }
	    //hist->Scale(correction_p[n]);
            //list.push_back(*hist);
	    //delete hist;
	}
	else if (!doLoop){
	    std::string string_eventSelection = eventSelectionString();
	    std::string string_weight = generateWeightString(false,timebin);
	    std::string string_triggered = isTriggerPassedString(triggerList_p,true);
	    drawHisto1D(tree, observable, string_eventSelection, string_weight, string_triggered, hist);
	    //std::string string_cut = "(" + string_eventSelection + ")*" + string_weight + "*" + string_triggered;
	    //std::cout << "Cut: " << string_cut<<std::endl;
	    //std::string string_obs_redirected = observable + " >> " + sampleList_p[n];
	    //tree->Draw(string_obs_redirected.c_str(), string_cut.c_str());
	    /*
	    ROOT::RDataFrame df(*tree);
	    std::cout << "A"<<std::endl;
	    df.Filter(string_eventSelection.c_str()).Filter(string_triggered.c_str()).Count();
	    std::cout << "B"<<std::endl;
	    //auto df2 = df.Define("weight_all", string_weight.c_str());
	    //auto hist = df.Histo1D(observable.c_str(), string_weight.c_str());
            std::cout << "C"<<std::endl;
	    auto hist = df.Histo1D(observable.c_str(), "weight_pu");
	    //auto hist = h->GetValue();
	    hist->DrawClone();
	    //auto hist = df.Fill<float>(TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin), {observable.c_str()});
	    std::cout << "DD"<<std::endl;
	    //hist->Scale(correction_p[n]);
            //std::cout << "D"<<std::endl;
            list2.push_back(*hist);
	    std::cout << "E"<<std::endl;
	    list2[n].Scale(correction_p[n]);
	    //delete hist;
	    */
	}
        hist->Scale(correction_p[n]);
        list.push_back(*hist);

        delete hist;
        //delete canvas;
        delete tree;
        delete file;
    }
    std::cout << "E"<<std::endl;
    groupingMC(list, groupList_p, clean_p);
    std::cout << "F"<<std::endl;
    write(filename_p, list, option_p);
    std::cout << "Wrote "<<filename_p<<std::endl;
}

void Generator::generateData(namelist            const& sampleList_p,
                             namelist            const& triggerList_p,
                             namelist            const& groupList_p,
                             std::vector<double> const& correction_p,
                             std::string         const& rootOption_p,
                             bool                       correctedLumi,
                             bool                       clean_p
                            )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::vector<TH1F> listUp;    
    std::vector<TH1F> listDown;

    double luminosityWeight = 1;

    std::string cleaned;
    if(!clean_p) cleaned = "_unclean";
    std::string filename_p = "./results/"+year+"/flattree/"+observable;
    if(correctedLumi){
            filename_p += "_lumicorrected";
    }
    filename_p += cleaned+"_data.root";

    std::map<unsigned int, bool> isFilled;  
    for(size_t n = 0; n < sampleList_p.size(); ++n){
        bool is2016H = false;
        if(sampleList_p[n].find("Run2016H") != std::string::npos){
            is2016H = true;
            std::cout << "Run H !" << std::endl;
        }
        std::string filename = "./inputs/"+year+"/DATA/"+sampleList_p[n]+"/NtupleProducer/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);

	double lumiSumOfOneOverWeight = 0;
	std::vector<int> useForPlot;

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            if (eventSelection(tree)!=1) continue;
            if(isTriggerPassed(tree, triggerList_p, is2016H)){
                if (sampleList_p[n].find("MuonEG") != std::string::npos){
                    isFilled.emplace(tree->GetLeaf("event")->GetValue(0),1);
		    lumiSumOfOneOverWeight += 1./tree->GetLeaf("lumi")->GetValue(0);
		    useForPlot.push_back(i);
                }
                if (sampleList_p[n].find("SingleElectron") != std::string::npos){
                    if(!isFilled[tree->GetLeaf("event")->GetValue(0)]){
                        isFilled.emplace(tree->GetLeaf("event")->GetValue(0),1);
                        lumiSumOfOneOverWeight += 1./tree->GetLeaf("lumi")->GetValue(0);
			useForPlot.push_back(i);
                    }
                }
                if (sampleList_p[n].find("SingleMuon") != std::string::npos){
                    if(!isFilled[tree->GetLeaf("event")->GetValue(0)]){
                        isFilled.emplace(tree->GetLeaf("event")->GetValue(0),1);
                        lumiSumOfOneOverWeight += 1./tree->GetLeaf("lumi")->GetValue(0);
			useForPlot.push_back(i);
                    }
                }
            }
            if(i % 100000 == 0)
                std::cout << "100 000 events passed" << std::endl;
        }
	double lumiavg_oneoverweight = lumiSumOfOneOverWeight / ((double)useForPlot.size());

	int nPassingEventsCheck = 0;
	double avg_luminosityWeight = 0;
        for(unsigned int i = 0; i < useForPlot.size(); ++i){
            tree->GetEntry(useForPlot[i]);
            if (correctedLumi) luminosityWeight = (1./tree->GetLeaf("lumi")->GetValue(0)) / lumiavg_oneoverweight;
	    else luminosityWeight = 1;
	    avg_luminosityWeight += luminosityWeight;
            hist->Fill(getObservableValue(tree),luminosityWeight);
	    //hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0),luminosityWeight);
	    nPassingEventsCheck++;
	}
	std::cout << "nPassingEventsCheck="<<nPassingEventsCheck<<" avg_luminosityWeight="<<avg_luminosityWeight/((double)nPassingEventsCheck)<< std::endl;

        hist->Scale(correction_p[n]);
        list.push_back(*hist);

        delete hist;
        delete canvas;
        delete tree;
        delete file;
    }
    groupingData(list, groupList_p, clean_p);
    //for (auto& x: isFilled)
    //    std::cout << " [" << x.first << ':' << x.second << ']' << '\n';
    write(filename_p, list, rootOption_p);
}

void Generator::generateDataTimed(namelist            const& sampleList_p,
                                  namelist            const& triggerList_p,
                                  namelist            const& groupList_p,
                                  std::vector<double> const& correction_p,
                                  int                        nBin_p,
                                  bool                       correctedLumi,
                                  bool                       clean_p
                                 )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;

    double luminosityWeight = 1;

    std::string cleaned;
    std::string lumicorr;
    if(correctedLumi) lumicorr = "_lumicorrected";
    if(!clean_p) cleaned = "_unclean";
    std::string filename_p = "./results/"+year+"/flattree/"+observable+ lumicorr+"_data_timed"+std::to_string(nBin_p)+cleaned+".root";
    
    std::map<unsigned int, bool> isFilled;  
    for(size_t n = 0; n < sampleList_p.size(); ++n){
        bool is2016H = false;
        if(sampleList_p[n].find("Run2016H") != std::string::npos){
            is2016H = true;
        }

        std::string filename = "./inputs/"+year+"/DATA/"+sampleList_p[n]+"/NtupleProducer/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        std::vector<TH1F*> hist(nBin_p);
        for(size_t i = 0; i < hist.size(); ++i){
            hist[i] = new TH1F((sampleList_p[n]+"_bin"+std::to_string(i)+'.').c_str(), (observable+"_bin"+std::to_string(i)+'.').c_str(), nBin, minBin, maxBin);
        }

        double lumiSumOfOneOverWeight = 0;
        std::vector<int> useForPlot;

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            if (eventSelection(tree)!=1) continue;
            //int whichBin = int(siderealHour(tree->GetLeaf("unix_time")->GetValue(0)))%nBin_p;
            if(isTriggerPassed(tree, triggerList_p, is2016H)){
                if (sampleList_p[n].find("MuonEG") != std::string::npos){
                    //hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                    isFilled.emplace(tree->GetLeaf("event")->GetValue(0),1);
		    lumiSumOfOneOverWeight += 1./tree->GetLeaf("lumi")->GetValue(0);
		    useForPlot.push_back(i);
                }
                if (sampleList_p[n].find("SingleElectron") != std::string::npos){
                    if(!isFilled[tree->GetLeaf("event")->GetValue(0)]){
                        //hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                        isFilled.emplace(tree->GetLeaf("event")->GetValue(0),1);
			lumiSumOfOneOverWeight += 1./tree->GetLeaf("lumi")->GetValue(0);
			useForPlot.push_back(i);
                    }
                }
                if (sampleList_p[n].find("SingleMuon") != std::string::npos){
                    if(!isFilled[tree->GetLeaf("event")->GetValue(0)]){
                        //hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                        isFilled.emplace(tree->GetLeaf("event")->GetValue(0),1);
			lumiSumOfOneOverWeight += 1./tree->GetLeaf("lumi")->GetValue(0);
			useForPlot.push_back(i);
                    }
                }
            }
            if(i % 100000 == 0)
                std::cout << "100 000 events passed" << std::endl;
        }
        double lumiavg_oneoverweight = lumiSumOfOneOverWeight / ((double)useForPlot.size());

        int nPassingEventsCheck = 0;
        double avg_luminosityWeight = 0;
        for(unsigned int i = 0; i < useForPlot.size(); ++i){
            int whichBin = int(siderealHour(tree->GetLeaf("unix_time")->GetValue(0)))%nBin_p;
            tree->GetEntry(useForPlot[i]);
            if (correctedLumi) luminosityWeight = (1./tree->GetLeaf("lumi")->GetValue(0)) / lumiavg_oneoverweight;
            else luminosityWeight = 1;
            avg_luminosityWeight += luminosityWeight;
            hist[whichBin]->Fill(getObservableValue(tree),luminosityWeight);
            //hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0),luminosityWeight);
            nPassingEventsCheck++;
        }
        std::cout << "nPassingEventsCheck="<<nPassingEventsCheck<<" avg_luminosityWeight="<<avg_luminosityWeight/((double)nPassingEventsCheck)<< std::endl;

        for(size_t i = 0; i < hist.size(); ++i){
            hist[i]->Scale(correction_p[n]);
            list.push_back(*hist[i]);
            delete hist[i];
        }
        delete canvas;
        delete tree;
        delete file;
    }
    for(int i = 0; i < nBin_p; ++i){
        groupingDataTimed(list, groupList_p, i, nBin_p, clean_p);
    }
    write(filename_p, list, "RECREATE");
}
