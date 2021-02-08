#include "debug.h"

#include <iostream>

void Log(double value){
    std::cout << "Time : " << value << "s" << std::endl;
}

void Log(std::string const& message){
    std::cout << message << std::endl;
}

void Log(int argc, char** argv){
    std::cout << std::endl;
    std::cout << "#####################" << std::endl;
    std::cout << "####  Aguments  #####" << std::endl;
    std::cout << "#####################" << std::endl;
    std::cout << "argc : " << argc << std::endl;
    for(int i = 0; i < argc; ++i){
            std::cout << "argc " << i << " : " << argv[i] << std::endl;
    }
    std::cout << "#####################" << std::endl;
    std::cout << std::endl;
}

double elapsedTime(std::clock_t t0, std::clock_t t1){
    return (float)(t1-t0)/CLOCKS_PER_SEC;
}