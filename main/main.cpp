#include "../src/generator.hpp"
#include "../src/debug.h"
#include "../src/sample_2017.hpp"

#include <ctime>


int main(int argc, char** argv){
    
    std::string observable; 
    std::clock_t t0 = std::clock();

    if(argc != 2){
        Log("Wrong number of arguments !!!");
        return 0;
    }
    else
         observable = argv[1];

    Log(argc, argv);
    Log(observable);

    Generator gen("2017");
    gen.generate("DATA", observable, 
                 sampleList_DATA_2017, triggerList, data_2017, systematicList, succedJobs_2017,
                 "RECREATE");
    gen.generate("MC", observable, 
                 sampleList_MC_2017, triggerList, ttbarList, systematicList, mc_rescale_2017,
                 "UPDATE");


    std::clock_t t1 = std::clock();
    Log(elapsedTime(t0,t1));

    return 0;
}