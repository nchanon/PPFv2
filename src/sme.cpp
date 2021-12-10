#include "sme.hpp"
#include "debug.h"

#include <fstream>
#include <iostream>

#include <TH1F.h>
#include <TF1.h>

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


SME::SME(Wilson wilson_p, int bin)
{

    wilson = wilson_p;

    if(wilson != Wilson::L and wilson != Wilson::R and
       wilson != Wilson::C and wilson != Wilson::D){
        Log("Error with Wilson coefficent choice");
        exit(0);
    }

    std::vector<double> matrix = generateMatrixPerMassBin(wilson, bin);
    Axx = matrix[0];
    Azz = matrix[10];
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

    std::string cutName = "13TeVCMSnew";
    std::string suffix = "_inc";

    std::string File_qqbar = "./inputs/pheno/" + cutName + "Pqqbar" + suffix + ".txt";
    std::string File_gg = "./inputs/pheno/" + cutName + "P2g" + suffix + ".txt";
    std::string File_F = "./inputs/pheno/" + cutName + "F" + suffix + ".txt";

    double nPqq              = readNumberOfEvents(File_qqbar);
    std::vector<double> mPqq = readElementMatrix(File_qqbar);
    double nP2g              = readNumberOfEvents(File_gg);
    std::vector<double> mP2g = readElementMatrix(File_gg);
    double nF                = readNumberOfEvents(File_F);
    std::vector<double> mF   = readElementMatrix(File_F);

    //double nPqq              = readNumberOfEvents("./inputs/pheno/13TeVCMSPqqbar.txt");
    //std::vector<double> mPqq = readElementMatrix("./inputs/pheno/13TeVCMSPqqbar.txt");
    //double nP2g              = readNumberOfEvents("./inputs/pheno/13TeVCMSP2g.txt");
    //std::vector<double> mP2g = readElementMatrix("./inputs/pheno/13TeVCMSP2g.txt");
    //double nF                = readNumberOfEvents("./inputs/pheno/13TeVCMSF.txt");
    //std::vector<double> mF   = readElementMatrix("./inputs/pheno/13TeVCMSF.txt");

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

std::vector<double> SME::generateMatrixPerMassBin(Wilson wilson_p, int bin) const
{
    
    std::string cutName = "13TeVCMSnew";
    std::string suffix = "_" + std::to_string(bin);
    
    std::string File_qqbar = "./inputs/pheno/" + cutName + "Pqqbar" + suffix + ".txt";
    std::string File_gg = "./inputs/pheno/" + cutName + "P2g" + suffix + ".txt";
    std::string File_F = "./inputs/pheno/" + cutName + "F" + suffix + ".txt";
    
    double nPqq              = readNumberOfEvents(File_qqbar);
    std::vector<double> mPqq = readElementMatrix(File_qqbar);
    double nP2g              = readNumberOfEvents(File_gg);
    std::vector<double> mP2g = readElementMatrix(File_gg);
    double nF                = readNumberOfEvents(File_F);
    std::vector<double> mF   = readElementMatrix(File_F);
    
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

double SME::fXX_hours(double* x, double* par) const 
{
    double T = x[0];
    double t0 = par[0];
    double t = t0+T;
    double part1 = ((a1()-a2())/2)*cos(2*OMEGA_GMST*t);
    double part2 = a3()*sin(2*OMEGA_GMST*t);
    return 2*(part1+part2);
}



/////////////////////////////////
// Public methods
/////////////////////////////////

void SME::generateModulation(int t0, int nBin)
{

    std::cout << "t0="<< t0 << "nBin="<<nBin<<std::endl;

    std::string wilsonName;
    if(wilson == Wilson::L)
        wilsonName = "cL";
    else if(wilson == Wilson::R)
        wilsonName = "cR";
    else if(wilson == Wilson::C)
        wilsonName = "c";
    else if(wilson == Wilson::D)
        wilsonName = "d";

    TH1F *hXX = new TH1F((wilsonName+"XX").c_str(), (wilsonName+"XX").c_str(), nBin, 0, nBin);
    TH1F *hXY = new TH1F((wilsonName+"XY").c_str(), (wilsonName+"XY").c_str(), nBin, 0, nBin);
    TH1F *hXZ = new TH1F((wilsonName+"XZ").c_str(), (wilsonName+"XZ").c_str(), nBin, 0, nBin);
    TH1F *hYZ = new TH1F((wilsonName+"YZ").c_str(), (wilsonName+"YZ").c_str(), nBin, 0, nBin);
    for(int i = 0; i < nBin; ++i){
        hXX->SetBinContent(i+1, fXX((t0+i)%24*3600));
        hXY->SetBinContent(i+1, fXY((t0+i)%24*3600));
        hXZ->SetBinContent(i+1, fXZ((t0+i)%24*3600));
        hYZ->SetBinContent(i+1, fYZ((t0+i)%24*3600));
    }

    //TF1* fXX = new TF1(("function_"+wilsonName+"XX").c_str(), fXX_hours,0,24);
    //fXX->SetParameter(0, t0);
    TH1F *hXX_details = new TH1F((wilsonName+"XX_details").c_str(), (wilsonName+"XX_details").c_str(), nBin*100, 0, nBin);
    TH1F *hXY_details = new TH1F((wilsonName+"XY_details").c_str(), (wilsonName+"XY_details").c_str(), nBin*100, 0, nBin);
    TH1F *hXZ_details = new TH1F((wilsonName+"XZ_details").c_str(), (wilsonName+"XZ_details").c_str(), nBin*100, 0, nBin);
    TH1F *hYZ_details = new TH1F((wilsonName+"YZ_details").c_str(), (wilsonName+"YZ_details").c_str(), nBin*100, 0, nBin);

    for (int i=0; i<nBin*3600; i++){
	hXX_details->SetBinContent(i+1, fXX(t0+i*3600/100));
        hXY_details->SetBinContent(i+1, fXY(t0+i*3600/100));
        hXZ_details->SetBinContent(i+1, fXZ(t0+i*3600/100));
        hYZ_details->SetBinContent(i+1, fYZ(t0+i*3600/100));
    }

    hXX->Write();
    hXY->Write();
    hXZ->Write();
    hYZ->Write();

    hXX_details->Write();
    hXY_details->Write();
    hXZ_details->Write();
    hYZ_details->Write();

}

void SME::generateModulationPerMassBin(int t0, int nBin, int binMass)
{
    std::string wilsonName;
    if(wilson == Wilson::L)
        wilsonName = "cL";
    else if(wilson == Wilson::R)
        wilsonName = "cR";
    else if(wilson == Wilson::C)
        wilsonName = "c";
    else if(wilson == Wilson::D)
        wilsonName = "d";

    TH1F *hXX = new TH1F((wilsonName+"XX_"+std::to_string(binMass)).c_str(), (wilsonName+"XX_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    TH1F *hXY = new TH1F((wilsonName+"XY_"+std::to_string(binMass)).c_str(), (wilsonName+"XY_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    TH1F *hXZ = new TH1F((wilsonName+"XZ_"+std::to_string(binMass)).c_str(), (wilsonName+"XZ_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    TH1F *hYZ = new TH1F((wilsonName+"YZ_"+std::to_string(binMass)).c_str(), (wilsonName+"YZ_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    for(int i = 0; i < nBin; ++i){
        hXX->SetBinContent(i+1, fXX((t0+i)%24*3600));
        hXY->SetBinContent(i+1, fXY((t0+i)%24*3600));
        hXZ->SetBinContent(i+1, fXZ((t0+i)%24*3600));
        hYZ->SetBinContent(i+1, fYZ((t0+i)%24*3600));
    }

    TH1F *hXX_details = new TH1F((wilsonName+"XX_details"+std::to_string(binMass)).c_str(), (wilsonName+"XX_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);
    TH1F *hXY_details = new TH1F((wilsonName+"XY_details"+std::to_string(binMass)).c_str(), (wilsonName+"XY_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);
    TH1F *hXZ_details = new TH1F((wilsonName+"XZ_details"+std::to_string(binMass)).c_str(), (wilsonName+"XZ_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);
    TH1F *hYZ_details = new TH1F((wilsonName+"YZ_details"+std::to_string(binMass)).c_str(), (wilsonName+"YZ_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);

    for (int i=0; i<nBin*3600; i++){
        hXX_details->SetBinContent(i+1, fXX(t0+i*3600/100));
        hXY_details->SetBinContent(i+1, fXY(t0+i*3600/100));
        hXZ_details->SetBinContent(i+1, fXZ(t0+i*3600/100));
        hYZ_details->SetBinContent(i+1, fYZ(t0+i*3600/100));
    }
    
    hXX->Write();
    hXY->Write();
    hXZ->Write();
    hYZ->Write();

    hXX_details->Write();
    hXY_details->Write();
    hXZ_details->Write();
    hYZ_details->Write();

}

