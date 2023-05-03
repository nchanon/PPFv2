#include "../src/sme.hpp"
#include "../src/debug.h"

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>


int main(int argc, char** argv){

    std::string obs;

    if(argc != 2){
        Log("Wrong number of arguments !!!");
        Log(" missing arguments : observable, year and option");
        return 0;
    }
    else{
         obs = argv[1];
         //year = argv[2];
         //launch = argv[3];
    }

    int nKinBin = -1;
    if (obs=="n_bjets") nKinBin = 6;
    if (obs=="pt_emu" || obs=="m_dilep") nKinBin = 8;

    std::cout << "Observable: "<<obs << std::endl;



    //int binage = 86400;
    int binage = 24;
    std::string output;
    std::vector<std::string> year{"2016", "2017"};

    int t0 = 0; //Choose 1st Jan 2016 for both 2016 and 2017?

    for(size_t i = 0; i < year.size(); ++i)
    {    
	if (year[i]=="2016") t0 = 1451606400;
	//if (year[i]=="2017") t0 = 1483228800;
	if (year[i]=="2017") t0 = 1451606400;

        if(binage != 24)
            output = "./results/"+year[i]+"/flattree/"+obs+"_sme"+std::to_string(binage)+".root";
        else
            output = "./results/"+year[i]+"/flattree/"+obs+"_sme.root";


        TFile * file = new TFile(output.c_str(), "RECREATE");

	//ttbar production and decay
        std::vector<SME> sme{
            SME(Wilson::L, obs, false),
            SME(Wilson::R, obs, false),
            SME(Wilson::C, obs, false),
            SME(Wilson::D, obs, false)
        };
        for(size_t j = 0; j < 4; ++j){
            sme[j].generateModulation(t0, binage, false);
        }

	for (int m=0; m<nKinBin; m++){
	  std::vector<SME> smePerMassBin{
              SME(Wilson::L, m, obs, false, "particle"),
              SME(Wilson::R, m, obs, false, "particle"),
              SME(Wilson::C, m, obs, false, "particle"),
              SME(Wilson::D, m, obs, false, "particle")
          };
	  for(size_t j = 0; j < 4; ++j){
            smePerMassBin[j].generateModulationPerMassBin(t0, binage, m, false, "particle");
          }
	}

        for (int m=1; m<=4; m++){
          std::vector<SME> smePerMassBin{
              SME(Wilson::L, m, obs, false, "reco_"+year[i]),
              SME(Wilson::R, m, obs, false, "reco_"+year[i]),
              SME(Wilson::C, m, obs, false, "reco_"+year[i]),
              SME(Wilson::D, m, obs, false, "reco_"+year[i])
          };
          for(size_t j = 0; j < 4; ++j){
            smePerMassBin[j].generateModulationPerMassBin(t0, binage, m, false, "reco_"+year[i]);
          }
        }

        for (int m=0; m<nKinBin; m++){
          std::vector<SME> smePerMassBin{
              SME(Wilson::L, m, obs, false, "signal_miniAOD_LO_"+year[i]),
              SME(Wilson::R, m, obs, false, "signal_miniAOD_LO_"+year[i]),
              SME(Wilson::C, m, obs, false, "signal_miniAOD_LO_"+year[i]),
              SME(Wilson::D, m, obs, false, "signal_miniAOD_LO_"+year[i])
          };
          for(size_t j = 0; j < 4; ++j){
            smePerMassBin[j].generateModulationPerMassBin(t0, binage, m, false, "signal_miniAOD_LO_"+year[i]);
          }
        }

        for (int m=0; m<nKinBin; m++){
          std::vector<SME> smePerMassBin{
              SME(Wilson::L, m, obs, false, "signal_miniAOD_LO_Comb"),
              SME(Wilson::R, m, obs, false, "signal_miniAOD_LO_Comb"),
              SME(Wilson::C, m, obs, false, "signal_miniAOD_LO_Comb"),
              SME(Wilson::D, m, obs, false, "signal_miniAOD_LO_Comb")
          };
          for(size_t j = 0; j < 4; ++j){
            smePerMassBin[j].generateModulationPerMassBin(t0, binage, m, false, "signal_miniAOD_LO_Comb");
          }
        }

	//single top decay
        std::vector<SME> sme_singletop{
            SME(Wilson::L, obs, true),
            SME(Wilson::R, obs, true),
            SME(Wilson::C, obs, true),
            SME(Wilson::D, obs, true)
        };
        for(size_t j = 0; j < 4; ++j){
            sme_singletop[j].generateModulation(t0, binage, true);
        }

        for (int m=0; m<nKinBin; m++){
          std::vector<SME> smePerMassBin_singletop{
              SME(Wilson::L, m, obs, true, "particle"),
              SME(Wilson::R, m, obs, true, "particle"),
              SME(Wilson::C, m, obs, true, "particle"),
              SME(Wilson::D, m, obs, true, "particle")
          };
          for(size_t j = 0; j < 4; ++j){
            smePerMassBin_singletop[j].generateModulationPerMassBin(t0, binage, m, true, "particle");
          }
        }
	

        file->Close();
        delete file;
    }
    std::cout << "Fini !" << std::endl;
    return 0;
}
