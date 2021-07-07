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

    for(size_t i = 0; i < year.size(); ++i)
    {    
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
        for(size_t i = 0; i < 4; ++i){
            sme[i].generateModulation(binage);
        }

        file->Close();
        delete file;
    }
    std::cout << "Fini !" << std::endl;
    return 0;
}
