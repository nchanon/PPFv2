#include "../src/generator.hpp"
#include "../src/debug.h"
#include "../src/sample_2017.hpp"

#include <ctime>


int main(int argc, char** argv){
    
    std::string observable; 
    bool isClean = true;
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

    std::string launch;
    do{
        std::cout << "Run on ? All, MC, Timed : ";
        std::cin >> launch;
    }
    while(launch != "All" and launch != "Timed" and launch != "MC");

    Generator gen(observable, binning, "2017");

    if(launch == "MC"){
        gen.generateMC(sampleList_MC_2017, triggerList, ttbarList, systematicList,
                       mc_rescale_2017, "RECREATE", isClean);
        gen.generateData(sampleList_DATA_2017, triggerList, data_2017, 
                         succedJobs_2017, "UPDATE", isClean);      
    }
    else if(launch == "Timed"){
        gen.generateDataTimmed(sampleList_DATA_2017, triggerList, data_2017, 
                               succedJobs_2017, 24, isClean);
    }
    else{
        gen.generateMC(sampleList_MC_2017, triggerList, ttbarList, systematicList,
                       mc_rescale_2017, "RECREATE", isClean);
        gen.generateData(sampleList_DATA_2017, triggerList, data_2017, 
                         succedJobs_2017, "UPDATE", isClean);           
        gen.generateDataTimmed(sampleList_DATA_2017, triggerList, data_2017, 
                               succedJobs_2017, 24, isClean);
    }


    std::clock_t t1 = std::clock();
    Log(elapsedTime(t0,t1));

    return 0;
}