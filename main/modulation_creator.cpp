#include "../src/sme.hpp"
#include "../src/debug.h"

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>


int main(){

    std::vector<std::string> year{"2016", "2017"};

    for(size_t i = 0; i < year.size(); ++i)
    {    
        TFile * file = new TFile(("./results/"+year[i]+"/flattree/sme.root").c_str(), "RECREATE");
        std::vector<SME> sme{
            SME(Wilson::L),
            SME(Wilson::R),
            SME(Wilson::C),
            SME(Wilson::D)
        };
        for(size_t i = 0; i < 4; ++i){
            sme[i].generateModulation(2);
        }

        file->Close();
        delete file;
    }
    
    return 0;
}
