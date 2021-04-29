#include "../src/amunu_analyser.hpp"


unsigned int nF      = 2084528;
unsigned int nPqqbar =  237058;
unsigned int nP2g    = 1847470;

float ratioQqbar = float(nPqqbar)/nF;
float ratio2g    = float(nP2g)/nF;

int main()
{

    std::vector<std::array<std::array<double, 4>, 4> > vmF = readVecMatrix("VMAF13TeVCMS", 5000000);
    std::vector<std::array<std::array<double, 4>, 4> > vmP2g = readVecMatrix("VMAP2g13TeVCMS", 5000000);
    std::vector<std::array<std::array<double, 4>, 4> > vmPqqbar = readVecMatrix("VMAPqqbar13TeVCMS", 5000000);

    std::vector<std::array<std::array<double, 4>, 4> > vmTotal(vmF.size());
    for(size_t i = 0; i < vmTotal.size(); ++i){
        for(size_t j = 0; j < vmTotal[i].size(); ++j){
            for(size_t k = 0; k < vmTotal[i][j].size(); ++k){
                vmTotal[i][j][k] = 0.5*ratio2g*vmP2g[i][j][k] 
                                  +0.5*ratioQqbar*vmPqqbar[i][j][k]
                                  +vmF[i][j][k];
            }
        }
    }
    


    int nbin = 100, minBin = -50, maxBin = 50;
    CreateTH1(vmP2g,"P2g",nbin, minBin, maxBin);
    CreateTH1(vmPqqbar,"Pqqbar",nbin, minBin, maxBin);
    CreateTH1(vmF,"F",nbin, minBin, maxBin);
    CreateTH1(vmTotal,"Total",nbin, minBin, maxBin);

    CreateComparaison(vmTotal, "fXX");


    return 0;
}