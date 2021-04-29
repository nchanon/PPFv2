#include "amunu_analyser.hpp"
#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <cmath>

constexpr double LATITUDE = 46.309/180*M_PI;
constexpr double AZIMUT   = 101.2790/180*M_PI;
constexpr double OMEGA_GMST = 7.2722e-5;

std::vector<std::array<std::array<double, 4>, 4> > readVecMatrix(std::string const& s, int n)
{
    std::vector<std::array<std::array<double, 4>, 4> > foo(n);
    std::ifstream f("./inputs/pheno/VecMatrix/"+s+".txt");
    for(int i=0; i<n; i++){
        for(int j=0; j<4; j++)
            for(int k=0; k<4; k++){
                f>>foo[i][j][k];
            }
        if(i%1000000 == 0)
            std::cout << "1 000 000 events passed" << std::endl;
    }
    return foo;
}

void CreateTH1(std::vector<std::array<std::array<double, 4>, 4> > const& vm,
               std::string const& name, int nbin, int binMin, int binMax)
{

    std::vector<std::vector<TH1F*> > hist(4);
    for(size_t i = 0; i < hist.size(); ++i){
        hist[i] = std::vector<TH1F*>(4);
    }
    
    for(size_t i = 0; i < hist.size(); ++i){
       for(size_t j = 0; j < hist[i].size(); ++j){
           std::string histname = name +'_' + std::to_string(i)+ std::to_string(j);
           hist[i][j] = new TH1F(histname.c_str(), histname.c_str(), nbin, binMin, binMax);
           for(size_t k = 0; k < vm.size(); ++k){
               if(!isZero(vm[k]))
                hist[i][j]->Fill(vm[k][i][j]);
           }
       }
    }

    std::string rootfileName = "./results/pheno/base/amunu"+name+".root";
    TFile rootfile(rootfileName.c_str(), "RECREATE");
    for(size_t i = 0; i < hist.size(); ++i){
        for(size_t j = 0; j < hist[i].size(); ++j){
            hist[i][j]->Write();
        }
    }
    for(size_t i = 0; i < hist.size(); ++i){
       for(size_t j = 0; j < hist[i].size(); ++j){
           delete hist[i][j];
       }
    }
}

bool isZero(std::array<std::array<double, 4>, 4> const& matrix)
{
    for(size_t i = 0; i < matrix.size(); ++i){
        for(size_t j = 0; j < matrix[i].size(); ++j){
           if(matrix[i][j] != 0)
                return false; 
        }
    }
    return true;
}

void CreateComparaison(std::vector<std::array<std::array<double, 4>, 4> > const& matrix,
                       std::string const& wilson)
{
    double part1 = sin(AZIMUT)*sin(AZIMUT)*sin(LATITUDE)*sin(LATITUDE)+cos(LATITUDE)*cos(LATITUDE);
    double part2 = sin(LATITUDE)*sin(LATITUDE)*cos(AZIMUT)*cos(AZIMUT);

    TH1F *histo = new TH1F(wilson.c_str(), wilson.c_str(), 86400, 0, 24);


    for(size_t i = 0; i < matrix.size(); ++i){
        float Axx = matrix[i][0][0];
        float Azz = matrix[i][2][2];
        if(Axx == 0 and Azz == 0)
            continue;
        float a1 = part1*Axx + part2*Azz;
        float a2 = cos(AZIMUT)*cos(AZIMUT)*Axx + sin(AZIMUT)*sin(AZIMUT)*Azz;
        float a3 = sin(LATITUDE)*cos(AZIMUT)*sin(AZIMUT)*(Azz-Axx);
        float a4 = cos(LATITUDE)*sin(LATITUDE)*cos(AZIMUT)*cos(AZIMUT)*(Azz-Axx);
        float a5 = cos(LATITUDE)*cos(AZIMUT)*sin(AZIMUT)*(Azz-Axx);
        if(wilson == "fXX")
        {
            double part1 = (a1-a2/2)*sin(2*OMEGA_GMST*i); 
            double part2 = a3*cos(2*OMEGA_GMST*i);
            histo->SetBinContent(i+1 ,2*(part1+part2));
        }
        else if(wilson == "fXY")
        {
            double part1 = (a1-a2/2)*sin(2*OMEGA_GMST*i); 
            double part2 = a3*cos(2*OMEGA_GMST*i);
            histo->Fill(2*(part1-part2));
        }  
        else if(wilson == "fXZ")
        {
            histo->Fill(2*(a4*cos(OMEGA_GMST*i) + a5*sin(OMEGA_GMST*i)));
        }  
        else if(wilson == "fYZ")
        {
            histo->Fill(2*(a4*sin(OMEGA_GMST*i) - a5*cos(OMEGA_GMST*i)));
        }  
    }   

    std::string filename = "./results/pheno/base/amunu_"+wilson+".root";
    TFile file(filename.c_str(), "RECREATE");
    histo->Write();

    delete histo;
}
    


