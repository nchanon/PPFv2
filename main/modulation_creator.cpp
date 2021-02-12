#include "../src/sme.hpp"
#include "../src/sample_2017.hpp"
#include "../src/debug.h"

#include <iostream>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1F.h>


int main(int argc, char** argv){

    std::string year;

    if(argc != 2){
        Log("Wrong number of arguments !!!");
        return 0;
    }
    else{
         year       = argv[1];
    }

    
    TFile * file = new TFile(("./results/"+year+"/flattree/sme.root").c_str(), "RECREATE");
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

    return 0;
}
