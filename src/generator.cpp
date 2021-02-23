#include "generator.hpp"
#include "debug.h"

#include <iostream>
#include <cmath>

#include <TCanvas.h>
#include <TLeaf.h>
#include <TString.h>


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
{}

/////////////////////////////////
// Private methods
/////////////////////////////////

double Generator::siderealHour(double time_p)
{
    if(year == "2017")
        return (OMEGA_UTC * (time_p - T0_2017))/OMEGA_GMST;
    else if(year == "2016")
        return (OMEGA_UTC * (time_p - T0_2016))/OMEGA_GMST;
    else
        return 0;
}


double Generator::generateWeight(TTree *tree_p)
{
    double w = 1;
    w *= tree_p->GetLeaf("weight_pu")->GetValue();
    w *= tree_p->GetLeaf("weight_generator")->GetValue();
    w *= tree_p->GetLeaf("weight_top")->GetValue();
    //lepton
    w *= tree_p->GetLeaf("weight_sfe_id")->GetValue();
    w *= tree_p->GetLeaf("weight_sfe_reco")->GetValue();
    w *= tree_p->GetLeaf("weight_sfm_id")->GetValue();
    w *= tree_p->GetLeaf("weight_sfm_iso")->GetValue();
    //hadron
    for(int i = 0; i< tree_p->GetLeaf("n_bjets")->GetValue(); ++i){
        w *= tree_p->GetLeaf("weight_sfb")->GetValue();
    }
    //trigger
    w *= tree_p->GetLeaf("weight_sf_em_trig")->GetValue();
    return w;
}

double Generator::generateSystematics(TTree            * tree_p,
                                      std::string const& systematicName,
                                      bool               isUp
                                     )
{
    if(systematicName == "syst_pu"){
        if(isUp)
            return tree_p->GetLeaf("weight_pu_up")->GetValue(0)/tree_p->GetLeaf("weight_pu")->GetValue(0);
        else
            return tree_p->GetLeaf("weight_pu_down")->GetValue(0)/tree_p->GetLeaf("weight_pu")->GetValue(0);
    }
    else{
        double foo = tree_p->GetLeaf(systematicName.c_str())->GetValue(0);
        if(isUp)
            return (1+foo);
        else
            return (1-foo);
    }
}

bool Generator::isTriggerPassed(TTree         * tree_p,
                                namelist const& triggerList_p
                               )
{
    int trig = 0; 
    for(size_t i = 0; i < triggerList_p.size(); ++i){
        if(tree_p->GetLeaf(triggerList_p[i].c_str())->GetValue(0) == 1){
            ++trig;
        }
    }
    return(trig >= 1);
}

void Generator::write(std::string       const& filename,
                      std::vector<TH1F>      & listObject,
                      std::string       const& option_p
                     )
{
    TFile *output = new TFile(filename.c_str(), option_p.c_str());
    for(size_t i = 0; i < listObject.size(); ++i){
        listObject[i].Write();
    }
    output->Close();
}

void Generator::groupingMC(std::vector<TH1F>      & list,
                           namelist          const& groupList_p,                        bool                     clean
                          )
{
    double nbin = list[0].GetNbinsX();
    double min = list[0].GetXaxis()->GetXmin();
    double max = list[0].GetXaxis()->GetXmax();

    for(std::string group : groupList_p){
        TH1F h((group).c_str(), (group).c_str(), nbin, min, max);
        for(size_t i = 0; i < list.size(); ++i){
            if(TString(list[i].GetName()).Contains(group))
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
            for(size_t i = 0; i < list.size(); ++i){
                if(TString(list[i].GetName()).Contains(group) and TString(list[i].GetName()).Contains(syst))
                    h.Add(&list[i]);        
            }
            list.push_back(h); 
        }
    }
    if(clean)
        list.erase(list.begin(), list.end()-(groupList_p.size()*systematicsList_p.size()));
}


/////////////////////////////////
// Public methods
/////////////////////////////////

void Generator::generateMC(namelist            const& sampleList_p,
                           namelist            const& triggerList_p,
                           namelist            const& groupList_p,
                           namelist            const& systematicsList_p,
                           std::vector<double> const& correction_p,
                           std::string         const& option_p,
                           bool                       clean_p
                          )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::vector<TH1F> listUp;    
    std::vector<TH1F> listDown;

    std::string cleaned;
    if(!clean_p) cleaned = "_unclean";
    std::string filename_p = "./results/"+year+"/flattree/"+observable+cleaned+".root";

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/MC/"+sampleList_p[n]+"/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);
        std::vector<TH1F*> histUp(systematicsList_p.size());
        std::vector<TH1F*> histDown(systematicsList_p.size());
        for(size_t i = 0; i < systematicsList_p.size(); ++i){
            histUp[i]   = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+"Up").c_str(), (observable+"Up").c_str(), nBin, minBin, maxBin);
            histDown[i] = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+"Down").c_str(), (observable+"Down").c_str(), nBin, minBin, maxBin);
        }

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            double weight = generateWeight(tree);
            if(isTriggerPassed(tree, triggerList_p)){
                hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight);
                for(size_t j = 0; j < systematicsList_p.size(); ++j){
                    double systUp = generateSystematics(tree, systematicsList_p[j], true);
                    double systDOwn = generateSystematics(tree, systematicsList_p[j], false);
                    histUp[j]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight*systUp);
                    histDown[j]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0), weight*systDOwn);                        
                }
            }
            if(i % 100000 == 0)
                std::cout << "100 000 events passed" << std::endl;
        }
        hist->Scale(correction_p[n]);
        list.push_back(*hist);
        for(size_t i = 0; i < systematicsList_p.size(); ++i){
            histUp[i]->Scale(correction_p[n]);
            histDown[i]->Scale(correction_p[n]);
            listUp.push_back(*histUp[i]);
            listDown.push_back(*histDown[i]);
        }

        delete hist;
        for(size_t i = 0; i < systematicsList_p.size(); ++i){
            delete histUp[i];
            delete histDown[i];
        }
        delete canvas;
        delete tree;
        delete file;
    }
    groupingMC(list, groupList_p, clean_p);
    groupingSystematics(listUp, groupList_p, systematicsList_p, true, clean_p);    // isUp = true
    groupingSystematics(listDown, groupList_p, systematicsList_p, false, clean_p); // isUp = false
    write(filename_p, list, option_p);
    write(filename_p, listUp, "UPDATE");
    write(filename_p, listDown, "UPDATE");
}

void Generator::generateData(namelist            const& sampleList_p,
                             namelist            const& triggerList_p,
                             namelist            const& groupList_p,
                             std::vector<double> const& correction_p,
                             std::string         const& rootOption_p,
                             bool                       clean_p
                            )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::vector<TH1F> listUp;    
    std::vector<TH1F> listDown;

    std::string cleaned;
    if(!clean_p) cleaned = "_unclean";
    std::string filename_p = "./results/"+year+"/flattree/"+observable+cleaned+".root";

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/DATA/"+sampleList_p[n]+"/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable.c_str(), nBin, minBin, maxBin);

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            if(isTriggerPassed(tree, triggerList_p)){
                if (sampleList_p[n].find("MuonEG") != std::string::npos){
                    hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                    continue;
                    if (sampleList_p[n].find("SingleMuon") != std::string::npos){
                        hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                        continue;
                        if (sampleList_p[n].find("SingleElectron") != std::string::npos){
                            hist->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                            continue;
                        }
                    }
                }
            }
            if(i % 100000 == 0)
                std::cout << "100 000 events passed" << std::endl;
        }
        hist->Scale(correction_p[n]);
        list.push_back(*hist);

        delete hist;
        delete canvas;
        delete tree;
        delete file;
    }
    groupingData(list, groupList_p, clean_p);

    write(filename_p, list, rootOption_p);
}

void Generator::generateDataTimmed(namelist            const& sampleList_p,
                                   namelist            const& triggerList_p,
                                   namelist            const& groupList_p,
                                   std::vector<double> const& correction_p,
                                   int                        nBin_p,
                                   bool                       clean_p
                                  )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;

    std::string cleaned;
    if(!clean_p) cleaned = "_unclean";
    std::string filename_p = "./results/"+year+"/flattree/"+observable+"_data_timed"+std::to_string(nBin_p)+cleaned+".root";

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/DATA/"+sampleList_p[n]+"/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        std::vector<TH1F*> hist(nBin_p);
        for(size_t i = 0; i < hist.size(); ++i){
            hist[i] = new TH1F((sampleList_p[n]+"_bin"+std::to_string(i)+'.').c_str(), (observable+"_bin"+std::to_string(i)+'.').c_str(), nBin, minBin, maxBin);
        }

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            int whichBin = int(siderealHour(tree->GetLeaf("unix_time")->GetValue(0)))%nBin_p;
            if(isTriggerPassed(tree, triggerList_p)){
                if (sampleList_p[n].find("MuonEG") != std::string::npos){
                    hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                    continue;
                    if (sampleList_p[n].find("SingleMuon") != std::string::npos){
                        hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                        continue;
                        if (sampleList_p[n].find("SingleElectron") != std::string::npos){
                            hist[whichBin]->Fill(tree->GetLeaf(observable.c_str())->GetValue(0));
                            continue;
                        }
                    }
                }
            }
            if(i % 100000 == 0)
                std::cout << "100 000 events passed" << std::endl;
        }

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