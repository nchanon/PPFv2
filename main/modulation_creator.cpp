#include "../src/sme.hpp"
#include "../src/debug.h"

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>


int main(){

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
            output = "./results/"+year[i]+"/flattree/sme"+std::to_string(binage)+".root";
        else
            output = "./results/"+year[i]+"/flattree/sme.root";


        TFile * file = new TFile(output.c_str(), "RECREATE");

        std::vector<SME> sme{
            SME(Wilson::L),
            SME(Wilson::R),
            SME(Wilson::C),
            SME(Wilson::D)
        };
        for(size_t j = 0; j < 4; ++j){
            sme[j].generateModulation(t0, binage);
        }

	for (int m=0; m<8; m++){
	  std::vector<SME> smePerMassBin{
              SME(Wilson::L, m),
              SME(Wilson::R, m),
              SME(Wilson::C, m),
              SME(Wilson::D, m)
          };
	  for(size_t j = 0; j < 4; ++j){
            smePerMassBin[j].generateModulationPerMassBin(t0, binage, m);
          }
	}

        file->Close();
        delete file;
    }
    std::cout << "Fini !" << std::endl;
    return 0;
}
