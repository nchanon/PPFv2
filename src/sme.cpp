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

SME::SME(Wilson wilson_p, std::string observable, bool doSingleTop)
    : wilson(wilson_p)
{
    if(wilson != Wilson::L and wilson != Wilson::R and 
       wilson != Wilson::C and wilson != Wilson::D){
        Log("Error with Wilson coefficent choice");
        exit(0);
    }

    std::vector<double> matrix = generateMatrix(wilson, observable, doSingleTop);
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


SME::SME(Wilson wilson_p, int bin, std::string observable, bool doSingleTop)
{

    wilson = wilson_p;

    if(wilson != Wilson::L and wilson != Wilson::R and
       wilson != Wilson::C and wilson != Wilson::D){
        Log("Error with Wilson coefficent choice");
        exit(0);
    }

    std::vector<double> matrix = generateMatrixPerMassBin(wilson, bin, observable, doSingleTop);
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

std::vector<double> SME::generateMatrix(Wilson wilson_p, std::string observable, bool doSingleTop) const
{

    std::string cutName = "13TeVCMSnanoGEN";
    std::string suffix = "_inc_particle";

    std::string File_qqbar = "./inputs/pheno/signal_" + observable + "_" + cutName + "Pqqbar" + suffix + ".txt";
    std::string File_gg = "./inputs/pheno/signal_" + observable + "_" + cutName + "P2g" + suffix + ".txt";
    std::string File_F = "./inputs/pheno/signal_" + observable + "_" + cutName + "F" + suffix + ".txt";

    std::string File_F_singletop = "./inputs/pheno/singletop_" + observable + "_" + cutName + "F" + suffix + ".txt";
    std::cout << "doSingleTop=" << doSingleTop << std::endl;


    double nPqq              = readNumberOfEvents(File_qqbar);
    std::vector<double> mPqq = readElementMatrix(File_qqbar);
    double nP2g              = readNumberOfEvents(File_gg);
    std::vector<double> mP2g = readElementMatrix(File_gg);
    double nF                = readNumberOfEvents(File_F);
    std::vector<double> mF   = readElementMatrix(File_F);

    double nF_singletop                = readNumberOfEvents(File_F_singletop);
    std::vector<double> mF_singletop   = readElementMatrix(File_F_singletop);

    if (doSingleTop){
	for (unsigned int i=0; i<mPqq.size(); i++){
	    mPqq[i] = 0;
	    mP2g[i] = 0;
	    mF[i] = mF_singletop[i];
	    nF = nF_singletop;
	}
    }

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

std::vector<double> SME::generateMatrixPerMassBin(Wilson wilson_p, int bin, std::string observable, bool doSingleTop) const
{
   
    std::string cutName = "13TeVCMSnanoGEN";
    std::string suffix = "_" + std::to_string(bin) + "_particle"; 
    
    std::string File_qqbar = "./inputs/pheno/signal_" + observable + "_" + cutName + "Pqqbar" + suffix + ".txt";
    std::string File_gg = "./inputs/pheno/signal_" + observable + "_" + cutName + "P2g" + suffix + ".txt";
    std::string File_F = "./inputs/pheno/signal_" + observable + "_" + cutName + "F" + suffix + ".txt";
  
    std::string File_F_singletop = "./inputs/pheno/singletop_" + observable + "_" + cutName + "F" + suffix + ".txt";
    std::cout << "doSingleTop=" << doSingleTop << std::endl;

    double nPqq              = readNumberOfEvents(File_qqbar);
    std::vector<double> mPqq = readElementMatrix(File_qqbar);
    double nP2g              = readNumberOfEvents(File_gg);
    std::vector<double> mP2g = readElementMatrix(File_gg);
    double nF                = readNumberOfEvents(File_F);
    std::vector<double> mF   = readElementMatrix(File_F);
  
    double nF_singletop                = readNumberOfEvents(File_F_singletop);
    std::vector<double> mF_singletop   = readElementMatrix(File_F_singletop);
 
    if (doSingleTop){
        for (unsigned int i=0; i<mPqq.size(); i++){
            mPqq[i] = 0;
            mP2g[i] = 0;
            mF[i] = mF_singletop[i];
            nF = nF_singletop;
        }
    }
 
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

double SME::fXX_primitive(double t) const
{
    double part1 = ((a1()-a2())/2)*sin(2*OMEGA_GMST*t)/(2*OMEGA_GMST);
    double part2 = a3()*cos(2*OMEGA_GMST*t)/(-2*OMEGA_GMST);
    return 2*(part1+part2);
}

double SME::fXY(double t) const
{
    double part1 = ((a1()-a2())/2)*sin(2*OMEGA_GMST*t); 
    double part2 = a3()*cos(2*OMEGA_GMST*t);
    return 2*(part1-part2);
}

double SME::fXY_primitive(double t) const
{   
    double part1 = ((a1()-a2())/2)*cos(2*OMEGA_GMST*t)/(-2*OMEGA_GMST);
    double part2 = a3()*sin(2*OMEGA_GMST*t)/(2*OMEGA_GMST);
    return 2*(part1-part2);
}

double SME::fXZ(double t) const{
    return 2*(a4()*cos(OMEGA_GMST*t) + a5()*sin(OMEGA_GMST*t));
}

double SME::fXZ_primitive(double t) const{
    return 2*(a4()*sin(OMEGA_GMST*t)/OMEGA_GMST + a5()*cos(OMEGA_GMST*t)/(-OMEGA_GMST));
}

double SME::fYZ(double t) const{
    return 2*(a4()*sin(OMEGA_GMST*t) - a5()*cos(OMEGA_GMST*t));
}

double SME::fYZ_primitive(double t) const{
    return 2*(a4()*cos(OMEGA_GMST*t)/(-OMEGA_GMST) - a5()*sin(OMEGA_GMST*t)/OMEGA_GMST);
}



/*
double SME::fXX_hours(double* x, double* par)  
{
    double T = x[0];
    double t0 = par[0];
    double t = t0+T;
    //return SME::fXX(t);
    double part1 = ((a1()-a2())/2)*cos(2*OMEGA_GMST*t);
    double part2 = a3()*sin(2*OMEGA_GMST*t);
    return 2*(part1+part2);
}
*/


/////////////////////////////////
// Public methods
/////////////////////////////////

void SME::generateModulation(int t0, int nBin, bool doSingleTop)
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

    std::string prefix = "";
    if (doSingleTop) prefix = "singletop_";


    TH1F *hXX = new TH1F((prefix+wilsonName+"XX").c_str(), (prefix+wilsonName+"XX").c_str(), nBin, 0, nBin);
    TH1F *hXY = new TH1F((prefix+wilsonName+"XY").c_str(), (prefix+wilsonName+"XY").c_str(), nBin, 0, nBin);
    TH1F *hXZ = new TH1F((prefix+wilsonName+"XZ").c_str(), (prefix+wilsonName+"XZ").c_str(), nBin, 0, nBin);
    TH1F *hYZ = new TH1F((prefix+wilsonName+"YZ").c_str(), (prefix+wilsonName+"YZ").c_str(), nBin, 0, nBin);
/*
    for(int i = 0; i < nBin; ++i){
        hXX->SetBinContent(i+1, SME::fXX((t0+i)%24*3600));
        hXY->SetBinContent(i+1, fXY((t0+i)%24*3600));
        hXZ->SetBinContent(i+1, fXZ((t0+i)%24*3600));
        hYZ->SetBinContent(i+1, fYZ((t0+i)%24*3600));
    }
*/
    //TF1* fXX = new TF1(("function_"+wilsonName+"XX").c_str(), fXX_hours,0,24);
    //fXX->SetParameter(0, t0);
    TH1F *hXX_details = new TH1F((prefix+wilsonName+"XX_details").c_str(), (prefix+wilsonName+"XX_details").c_str(), nBin*100, 0, nBin);
    TH1F *hXY_details = new TH1F((prefix+wilsonName+"XY_details").c_str(), (prefix+wilsonName+"XY_details").c_str(), nBin*100, 0, nBin);
    TH1F *hXZ_details = new TH1F((prefix+wilsonName+"XZ_details").c_str(), (prefix+wilsonName+"XZ_details").c_str(), nBin*100, 0, nBin);
    TH1F *hYZ_details = new TH1F((prefix+wilsonName+"YZ_details").c_str(), (prefix+wilsonName+"YZ_details").c_str(), nBin*100, 0, nBin);

    for (int i=0; i<nBin*3600; i++){
	hXX_details->SetBinContent(i+1, SME::fXX(t0+i*3600/100));
        hXY_details->SetBinContent(i+1, fXY(t0+i*3600/100));
        hXZ_details->SetBinContent(i+1, fXZ(t0+i*3600/100));
        hYZ_details->SetBinContent(i+1, fYZ(t0+i*3600/100));
    }

    for (int j=0; j<nBin; j++){
        hXX->SetBinContent(j+1,(fXX_primitive(t0+(j+1)*3600)-fXX_primitive(t0+j*3600))/3600.);
        hXY->SetBinContent(j+1,(fXY_primitive(t0+(j+1)*3600)-fXY_primitive(t0+j*3600))/3600.);
        hXZ->SetBinContent(j+1,(fXZ_primitive(t0+(j+1)*3600)-fXZ_primitive(t0+j*3600))/3600.);
        hYZ->SetBinContent(j+1,(fYZ_primitive(t0+(j+1)*3600)-fYZ_primitive(t0+j*3600))/3600.);
        hXX->SetBinError(j+1,0);
        hXY->SetBinError(j+1,0);
        hXZ->SetBinError(j+1,0);
        hYZ->SetBinError(j+1,0);
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

void SME::generateModulationPerMassBin(int t0, int nBin, int binMass, bool doSingleTop)
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

    std::string prefix = "";
    if (doSingleTop) prefix = "singletop_";

    TH1F *hXX = new TH1F((prefix+wilsonName+"XX_"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"XX_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    TH1F *hXY = new TH1F((prefix+wilsonName+"XY_"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"XY_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    TH1F *hXZ = new TH1F((prefix+wilsonName+"XZ_"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"XZ_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
    TH1F *hYZ = new TH1F((prefix+wilsonName+"YZ_"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"YZ_"+std::to_string(binMass)).c_str(), nBin, 0, nBin);
//    for(int i = 0; i < nBin; ++i){
//        hXX->SetBinContent(i+1, fXX((t0+i)%24*3600));
//        hXY->SetBinContent(i+1, fXY((t0+i)%24*3600));
//        hXZ->SetBinContent(i+1, fXZ((t0+i)%24*3600));
//        hYZ->SetBinContent(i+1, fYZ((t0+i)%24*3600));
//    }

    TH1F *hXX_details = new TH1F((prefix+wilsonName+"XX_details"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"XX_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);
    TH1F *hXY_details = new TH1F((prefix+wilsonName+"XY_details"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"XY_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);
    TH1F *hXZ_details = new TH1F((prefix+wilsonName+"XZ_details"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"XZ_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);
    TH1F *hYZ_details = new TH1F((prefix+wilsonName+"YZ_details"+std::to_string(binMass)).c_str(), (prefix+wilsonName+"YZ_details"+std::to_string(binMass)).c_str(), nBin*100, 0, nBin);

    //std::vector<double> vec_avg_XX;
    //std::vector<double> vec_avg_XY;
    //std::vector<double> vec_avg_XZ;
    //std::vector<double> vec_avg_YZ;
    //avg_XX[j]=0; 
    //avg_XY[j]=0, avg_XZ=0, avg_YZ=0;
    //int j=-1;
    for (int i=0; i<nBin*3600; i++){
	//if (i%100==0) j++;
	//avg_XX[j] += fXX(t0+i*3600/100);
	//avg_XY[j] += fXY(t0+i*3600/100);
	//avg_XZ[j] += fXZ(t0+i*3600/100);
	//avg_YZ[j] += fYZ(t0+i*3600/100);
	//hXX->Fill(0.5+j, fXX(t0+i*3600/100));
        //hXY->Fill(0.5+j, fXY(t0+i*3600/100));
        //hXZ->Fill(0.5+j, fXZ(t0+i*3600/100));
        //hYZ->Fill(0.5+j, fXZ(t0+i*3600/100));
        hXX_details->SetBinContent(i+1, fXX(t0+i*3600/100));
        hXY_details->SetBinContent(i+1, fXY(t0+i*3600/100));
        hXZ_details->SetBinContent(i+1, fXZ(t0+i*3600/100));
        hYZ_details->SetBinContent(i+1, fYZ(t0+i*3600/100));
    }

    for (int j=0; j<nBin; j++){
	hXX->SetBinContent(j+1,(fXX_primitive(t0+(j+1)*3600)-fXX_primitive(t0+j*3600))/3600.);
        hXY->SetBinContent(j+1,(fXY_primitive(t0+(j+1)*3600)-fXY_primitive(t0+j*3600))/3600.);
        hXZ->SetBinContent(j+1,(fXZ_primitive(t0+(j+1)*3600)-fXZ_primitive(t0+j*3600))/3600.);
        hYZ->SetBinContent(j+1,(fYZ_primitive(t0+(j+1)*3600)-fYZ_primitive(t0+j*3600))/3600.);
        hXX->SetBinError(j+1,0);
        hXY->SetBinError(j+1,0);
        hXZ->SetBinError(j+1,0);
        hYZ->SetBinError(j+1,0);
    }
/*
    for (int j=0; j<nBin; j++){
	int bin1 = 1+j*100;
	int bin2 = (j+1)*100;
	//hXX->SetBinContent(j+1, hXX_details->Integral(bin1, bin2)/100.);
        hXY->SetBinContent(j+1, hXY_details->Integral(bin1, bin2)/100.);
        hXZ->SetBinContent(j+1, hXZ_details->Integral(bin1, bin2)/100.);
        hYZ->SetBinContent(j+1, hYZ_details->Integral(bin1, bin2)/100.);
        hXX->SetBinError(j+1,0);
        hXY->SetBinError(j+1,0);
        hXZ->SetBinError(j+1,0);
        hYZ->SetBinError(j+1,0);
    }
*/
    hXX->Write();
    hXY->Write();
    hXZ->Write();
    hYZ->Write();

    hXX_details->Write();
    hXY_details->Write();
    hXZ_details->Write();
    hYZ_details->Write();

}

