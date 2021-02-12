#include "generator.hpp"
#include "debug.h"

#include <iostream>
#include <cmath>

#include <TCanvas.h>
#include <TLeaf.h>
#include <TString.h>


const double OMEGA_UT1  = 7.2921e-5;
const double OMEGA_GMST = 7.2722e-5;
const double T0         = 1483228800.;
const double PHASE      = 3.2830;


/////////////////////////////////
// Constructor
/////////////////////////////////

Generator::Generator(std::string const& year_p) : year(year_p)
{}

/////////////////////////////////
// Private methods
/////////////////////////////////

double Generator::siderealHour(double time_p)
{
    return (OMEGA_UT1 * (time_p - T0))/OMEGA_GMST;
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
                           namelist          const& groupList_p
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
}

void Generator::groupingData(std::vector<TH1F>      & list,
                             namelist          const& groupList_p
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
}

void Generator::groupingDataTimed(std::vector<TH1F>      & list,
                                  namelist          const& groupList_p,
                                  int                      bin
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
}

void Generator::groupingSystematics(std::vector<TH1F>      & list,
                                    namelist          const& groupList_p,
                                    namelist          const& systematicsList_p,
                                    bool                     isUp
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
}


/////////////////////////////////
// Public methods
/////////////////////////////////

void Generator::generateMC(std::string         const& observable_p,
                           namelist            const& sampleList_p,
                           namelist            const& triggerList_p,
                           namelist            const& groupList_p,
                           namelist            const& systematicsList_p,
                           std::vector<double> const& correction_p,
                           std::string         const& option_p
                          )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::vector<TH1F> listUp;    
    std::vector<TH1F> listDown;
    std::string filename_p = "./results/"+year+"/flattree/"+observable_p+".root";

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/MC/"+sampleList_p[n]+"/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable_p.c_str(), 25, 0, 500);
        std::vector<TH1F*> histUp(systematicsList_p.size());
        std::vector<TH1F*> histDown(systematicsList_p.size());
        for(size_t i = 0; i < systematicsList_p.size(); ++i){
            histUp[i]   = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+"Up").c_str(), (observable_p+"Up").c_str(), 25, 0, 500);
            histDown[i] = new TH1F((sampleList_p[n]+"_"+systematicsList_p[i]+"Down").c_str(), (observable_p+"Down").c_str(), 25, 0, 500);
        }

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            double weight = generateWeight(tree);
            if(isTriggerPassed(tree, triggerList_p)){
                hist->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0), weight);
                for(size_t j = 0; j < systematicsList_p.size(); ++j){
                    double systUp = generateSystematics(tree, systematicsList_p[j], true);
                    double systDOwn = generateSystematics(tree, systematicsList_p[j], false);
                    histUp[j]->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0), weight*systUp);
                    histDown[j]->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0), weight*systDOwn);                        
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
    groupingMC(list, groupList_p);
    groupingSystematics(listUp, groupList_p, systematicsList_p, true);    // isUp = true
    groupingSystematics(listDown, groupList_p, systematicsList_p, false); // isUp = false
    write(filename_p, list, option_p);
    write(filename_p, listUp, "UPDATE");
    write(filename_p, listDown, "UPDATE");
}

void Generator::generateData(std::string         const& observable_p,
                             namelist            const& sampleList_p,
                             namelist            const& triggerList_p,
                             namelist            const& groupList_p,
                             std::vector<double> const& correction_p,
                             std::string         const& rootOption_p
                            )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::vector<TH1F> listUp;    
    std::vector<TH1F> listDown;
    std::string filename_p = "./results/"+year+"/flattree/"+observable_p+".root";

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/DATA/"+sampleList_p[n]+"/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        TH1F* hist      = new TH1F(sampleList_p[n].c_str(), observable_p.c_str(), 25, 0, 500);

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            if(isTriggerPassed(tree, triggerList_p)){
                if (sampleList_p[n].find("MuonEG") != std::string::npos){
                    hist->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0));
                    continue;
                    if (sampleList_p[n].find("SingleMuon") != std::string::npos){
                        hist->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0));
                        continue;
                        if (sampleList_p[n].find("SingleElectron") != std::string::npos){
                            hist->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0));
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
    groupingData(list, groupList_p);
    write(filename_p, list, rootOption_p);
}

void Generator::generateDataTimmed(std::string         const& observable_p,
                                   namelist            const& sampleList_p,
                                   namelist            const& triggerList_p,
                                   namelist            const& groupList_p,
                                   std::vector<double> const& correction_p,
                                   int                        nBin_p
                                  )
{
    TH1F::SetDefaultSumw2(1);
    std::vector<TH1F> list;
    std::string filename_p = "./results/"+year+"/flattree/"+observable_p+"_data_timed"+std::to_string(nBin_p)+".root";

    for(size_t n = 0; n < sampleList_p.size(); ++n){
        std::string filename = "./inputs/"+year+"/DATA/"+sampleList_p[n]+"/tree.root";
        TFile* file = new TFile(filename.c_str());
        TTree *tree;
        file->GetObject("events", tree);
        TCanvas *canvas = new TCanvas(sampleList_p[n].c_str());
        std::vector<TH1F*> hist(nBin_p);
        for(size_t i = 0; i < hist.size(); ++i){
            hist[i] = new TH1F((sampleList_p[n]+"_bin"+std::to_string(i)+'.').c_str(), (observable_p+"_bin"+std::to_string(i)+'.').c_str(), 25, 0, 500);
        }

        std::cout << " -> " << sampleList_p[n] << std::endl;
        for(int i = 0; i < tree->GetEntriesFast(); ++i){
            tree->GetEntry(i);
            int whichBin = int(siderealHour(tree->GetLeaf("unix_time")->GetValue(0)))%nBin_p;
            if(isTriggerPassed(tree, triggerList_p)){
                if (sampleList_p[n].find("MuonEG") != std::string::npos){
                    hist[whichBin]->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0));
                    continue;
                    if (sampleList_p[n].find("SingleMuon") != std::string::npos){
                        hist[whichBin]->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0));
                        continue;
                        if (sampleList_p[n].find("SingleElectron") != std::string::npos){
                            hist[whichBin]->Fill(tree->GetLeaf(observable_p.c_str())->GetValue(0));
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
        groupingDataTimed(list, groupList_p, i);
    }
    write(filename_p, list, "RECREATE");
}