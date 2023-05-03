#pragma once 

#include <vector>
#include <string>
#include <cmath>

constexpr double LATITUDE = 46.309/180*M_PI;
constexpr double AZIMUT   = 101.2790/180*M_PI;
constexpr double OMEGA_GMST = 7.2722e-5;
constexpr double TILT = -atan(1.23 / 100);

enum class Wilson{
    L,
    R,
    C,
    D
};

class SME{

    private:

        Wilson wilson;
        double Axx;
        double Azz;

        double a1() const;
        double a2() const;
        double a3() const;
        double a4() const;
        double a5() const;

        std::vector<double> generateMatrix(Wilson wilson_p, std::string observable, bool doSingleTop) const;
        std::vector<double> generateMatrixPerMassBin(Wilson wilson_p, int bin, std::string observable, bool doSingleTop, std::string analysis_level) const;
        double readNumberOfEvents(std::string const& path_p) const;
        std::vector<double> readElementMatrix(std::string const& path_p) const;

        double fXX(double t) const;
        double fXY(double t) const;
        double fXZ(double t) const;
        double fYZ(double t) const;

	double fXX_primitive(double t) const;
        double fXY_primitive(double t) const;
        double fXZ_primitive(double t) const;
        double fYZ_primitive(double t) const;

	//static double fXX_hours(double* x, double* par);

    public:

        SME();
        SME(Wilson wilson_p, std::string observable, bool doSingleTop);
        SME(SME const& other);
        SME(Wilson wilson_p, int bin, std::string observable, bool doSingleTop,std::string analysis_level);

        SME &operator=(SME const& other);

        void generateModulation(int t0,int nBin = 24, bool doSingleTop = false);
	void generateModulationPerMassBin(int t0, int nBin, int binMass, bool doSingleTop = false, std::string analysis_level = "particle");

	//bool doSingleTop = false;
        //std::string observable;

        //double fXX_hours(double* x, double* par);

};
