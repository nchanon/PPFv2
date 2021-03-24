#include "../src/generator.hpp"
#include "../src/debug.h"
#include "../src/sample_2017.hpp"
#include "../src/sample_2016.hpp"
#include <ctime>

void VectorCheck(namelist v){
    for(size_t i = 0; i < v.size(); ++i){
        std::cout << v[i] << std::endl;
    }
    std::cout << std::endl;
}
void VectorCheck(std::vector<double> v){
    for(size_t i = 0; i < v.size(); ++i){
        std::cout << v[i] << std::endl;
    }
    std::cout << std::endl;
}





int main(int argc, char** argv){
    
    std::string observable; 
    std::string year; 
    bool isClean = true;
    std::clock_t t0 = std::clock();

    if(argc != 3){
        Log("Wrong number of arguments !!!");
        Log(" missing arguments : observable, year");
        return 0;
    }
    else{
         observable = argv[1];
         year = argv[2];
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
    else if(observable == "pt_lead"){
        binning[0] = 25; binning[1] = 0; binning[2] = 250;
    }
    else if(observable == "pt_sublead"){
        binning[0] = 25; binning[1] = 0; binning[2] = 250;
    }

// ------------------------------------------------------------- //

    namelist triggerList;
    // mc
    namelist sampleList_MC;
    std::vector<double> mc_rescale;
    // data
    namelist sampleList_DATA;
    namelist data;
    std::vector<double> succedJobs;
    
    if(year == "2016"){
        triggerList     = triggerList_2016;
        sampleList_MC   = sampleList_MC_2016;
        mc_rescale      = mc_rescale_2016;
        sampleList_DATA = sampleList_DATA_2016;
        data            = data_2016;
        succedJobs      = succedJobs_2016;
    } 
    else if(year == "2017"){        
        triggerList     = triggerList_2017;
        sampleList_MC   = sampleList_MC_2017;
        mc_rescale      = mc_rescale_2017;
        sampleList_DATA = sampleList_DATA_2017;
        data            = data_2017;
        succedJobs      = succedJobs_2017;
    }
    else{
        std::cout << "Bad year input (2016 or 2017)" << std::endl;
        return 0;
    }

    /*
    VectorCheck(sampleList_MC);
    VectorCheck(mc_rescale);
    VectorCheck(sampleList_DATA);
    VectorCheck(data);
    VectorCheck(succedJobs);
    VectorCheck(triggerList);
    VectorCheck(ttbarList);
    VectorCheck(systematicList);
    VectorCheck(systematicTimeList);
    VectorCheck(systematicRate);
    */

// ------------------------------------------------------------- //


    std::string launch;
    do{
        std::cout << "Run on ? All, MC, Timed : ";
        std::cin >> launch;
    }
    while(launch != "All" and launch != "Timed" and launch != "MC");

    Generator gen(observable, binning, year);

    if(launch == "MC"){
        gen.generateMC(sampleList_MC, triggerList, ttbarList, 
                       systematicList, systematicTimeList,
                       mc_rescale, "RECREATE", isClean);
        gen.generateData(sampleList_DATA, triggerList, data, 
                         succedJobs, "UPDATE", isClean);      
    }
    else if(launch == "Timed"){
        gen.generateDataTimed(sampleList_DATA, triggerList, data, 
                               succedJobs, 24, isClean);
    }
    else{
        gen.generateMC(sampleList_MC, triggerList, ttbarList, 
                       systematicList,systematicTimeList,
                       mc_rescale, "RECREATE", isClean);
        gen.generateData(sampleList_DATA, triggerList, data, 
                         succedJobs, "UPDATE", isClean);           
        gen.generateDataTimed(sampleList_DATA, triggerList, data, 
                               succedJobs, 24, isClean);
    }


    std::clock_t t1 = std::clock();
    Log(elapsedTime(t0,t1));



    return 0;
}