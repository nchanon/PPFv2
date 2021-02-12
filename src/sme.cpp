#include "sme.hpp"
#include "debug.h"

#include <fstream>
#include <iostream>

#include <TH1F.h>

/////////////////////////////////
// Constructor and Operators
/////////////////////////////////

SME::SME()
{}

SME::SME(Wilson wilson_p)
    : wilson(wilson_p)
{
    if(wilson != Wilson::L and wilson != Wilson::R and 
       wilson != Wilson::C and wilson != Wilson::D){
        Log("Error with Wilson coefficent choice");
        exit(0);
    }

    std::vector<double> matrix = generateMatrix(wilson);
    Axx = matrix[0];
    Azz = matrix[10];
}

SME::SME(SME const& other)
{
    wilson = other.wilson;
    Axx    = other.Axx;
    Azz    = other.Azz;
}

SME &SME::operator=(SME const& other)
{
    wilson = other.wilson;
    Axx    = other.Axx;
    Azz    = other.Azz;
    return *this;
}


/////////////////////////////////
// Private methods
/////////////////////////////////

double SME::a1() const
{
    double part1 = sin(AZIMUT)*sin(AZIMUT)*sin(LATITUDE)*sin(LATITUDE)+cos(LATITUDE)*cos(LATITUDE);
    double part2 = sin(LATITUDE)*sin(LATITUDE)*cos(AZIMUT)*cos(AZIMUT);
    return  part1*Axx + part2*Azz;
}

double SME::a2() const
{
    return cos(AZIMUT)*cos(AZIMUT)*Axx + sin(AZIMUT)*sin(AZIMUT)*Azz;
}

double SME::a3() const
{
    return sin(LATITUDE)*cos(AZIMUT)*sin(AZIMUT)*(Azz-Axx);
}

double SME::a4() const
{
    return cos(LATITUDE)*sin(LATITUDE)*cos(AZIMUT)*cos(AZIMUT)*(Azz-Axx);
}

double SME::a5() const
{
    return cos(LATITUDE)*cos(AZIMUT)*sin(AZIMUT)*(Azz-Axx);
}

std::vector<double> SME::generateMatrix(Wilson wilson_p) const
{
    double nPqq              = readNumberOfEvents("./inputs/pheno/13TeVCMSPqqbar.txt");
    std::vector<double> mPqq = readElementMatrix("./inputs/pheno/13TeVCMSPqqbar.txt");
    double nP2g              = readNumberOfEvents("./inputs/pheno/13TeVCMSP2g.txt");
    std::vector<double> mP2g = readElementMatrix("./inputs/pheno/13TeVCMSP2g.txt");
    double nF                = readNumberOfEvents("./inputs/pheno/13TeVCMSF.txt");
    std::vector<double> mF   = readElementMatrix("./inputs/pheno/13TeVCMSF.txt");



    std::vector<double> matrix(mF.size());

    for(size_t i = 0; i < matrix.size(); ++i){
        if(wilson_p == Wilson::L)
            matrix[i] = 0.5*(nPqq/nF)*mPqq[i] + 0.5*(nP2g/nF)*mP2g[i] + mF[i];
        else if(wilson_p == Wilson::R)
            matrix[i] = 0.5*(nPqq/nF)*mPqq[i] + 0.5*(nP2g/nF)*mP2g[i];
        else if(wilson_p == Wilson::C)
            matrix[i] = (nPqq/nF)*mPqq[i] + (nP2g/nF)*mP2g[i] + mF[i];
        else if(wilson_p == Wilson::D)
            matrix[i] = mF[i];
    }
    return matrix;
}

double SME::readNumberOfEvents(std::string const& path_p) const
{
    std::string tmp;
    double foo = 0;

    std::ifstream file(path_p);
    file >> tmp >> tmp >> tmp;
    file >> foo;
    file.close();
    return foo;
}


std::vector<double> SME::readElementMatrix(std::string const& path_p) const
{
    std::vector<double> matrix(16);
    std::string tmp;

    std::ifstream file(path_p);
    file >> tmp >> tmp >> tmp >> tmp;
    for(size_t i = 0; i < matrix.size(); ++i){
        file >> matrix[i];
    }
    file.close();
    return matrix;
}

double SME::fXX(double t) const
{
    double part1 = ((a1()-a2())/2)*cos(2*OMEGA_GMST*t); 
    double part2 = a3()*sin(2*OMEGA_GMST*t);
    return 2*(part1+part2);
}

double SME::fXY(double t) const
{
    double part1 = ((a1()-a2())/2)*sin(2*OMEGA_GMST*t); 
    double part2 = a3()*cos(2*OMEGA_GMST*t);
    return 2*(part1-part2);
}

double SME::fXZ(double t) const{
    return 2*(a4()*cos(OMEGA_GMST*t) + a5()*sin(OMEGA_GMST*t));
}

double SME::fYZ(double t) const{
    return 2*(a4()*sin(OMEGA_GMST*t) - a5()*cos(OMEGA_GMST*t));
}

/////////////////////////////////
// Public methods
/////////////////////////////////

void SME::generateModulation(int t0, int nBin)
{
    std::string wilsonName;
    if(wilson == Wilson::L)
        wilsonName = "L";
    else if(wilson == Wilson::R)
        wilsonName = "R";
    else if(wilson == Wilson::C)
        wilsonName = "C";
    else if(wilson == Wilson::D)
        wilsonName = "D";

    TH1F *hXX = new TH1F(("fXX_"+wilsonName).c_str(), ("fXX_"+wilsonName).c_str(), nBin, 0, nBin);
    TH1F *hXY = new TH1F(("fXY_"+wilsonName).c_str(), ("fXY_"+wilsonName).c_str(), nBin, 0, nBin);
    TH1F *hXZ = new TH1F(("fXZ_"+wilsonName).c_str(), ("fXZ_"+wilsonName).c_str(), nBin, 0, nBin);
    TH1F *hYZ = new TH1F(("fYZ_"+wilsonName).c_str(), ("fYZ_"+wilsonName).c_str(), nBin, 0, nBin);
    for(int i = 0; i < nBin; ++i){
        hXX->SetBinContent(i+1, fXX((t0+i)%24*3600));
        hXY->SetBinContent(i+1, fXY((t0+i)%24*3600));
        hXZ->SetBinContent(i+1, fXZ((t0+i)%24*3600));
        hYZ->SetBinContent(i+1, fYZ((t0+i)%24*3600));
    }
    hXX->Write();
    hXY->Write();
    hXZ->Write();
    hYZ->Write();
}

