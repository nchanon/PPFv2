#include "../src/generator.hpp"
#include "../src/debug.h"
#include "../src/sample_2017.hpp"
#include "../src/sample_2016.hpp"
#include <ctime>


int main()
{

   std::vector<int> binning(3);
        binning[0] = 25; binning[1] = 0; binning[2] = 500;
    std::string observable = "m_dilep"; 
    std::string year = "2017"; 
    Generator gen(observable,binning, year);

    return 0;
}