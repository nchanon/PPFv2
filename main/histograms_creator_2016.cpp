#include "../src/generator.hpp"
#include "../src/debug.h"
#include "../src/sample_2016.hpp"

#include <ctime>


int main(int argc, char** argv){
    
    std::string observable; 

    std::clock_t t0 = std::clock();

    if(argc != 2){
        Log("Wrong number of arguments !!!");
        Log(" missing arguments : observable");
        return 0;
    }
    else{
         observable = argv[1];
    }
    Log(argc, argv);
    Log(observable);

    std::vector<int> binning(3);
    if(observable == "m_dilep"){
        binning[0] = 25; binning[1] = 0; binning[2] = 500;
    }
    else if(observable == "n_bjets"){
        binning[0] = 5; binning[1] = 0; binning[2] = 5;
    }


// ------------------------------------------------------------- //


    Generator gen(observable, binning, "2016");

    gen.generateMC(sampleList_MC_2016, triggerList, ttbarList, systematicList, mc_rescale_2016,
                   "RECREATE");
    gen.generateData(sampleList_DATA_2016, triggerList, data_2016, succedJobs_2016,
                     "UPDATE");
    gen.generateDataTimmed(sampleList_DATA_2016, triggerList, data_2016, succedJobs_2016,
                           24);

    std::clock_t t1 = std::clock();
    Log(elapsedTime(t0,t1));

    return 0;
}