#include "../src/card.hpp"
#include "../src/debug.h"
#include "../src/sample_2017.hpp"
#include "../src/sample_2016.hpp"

#include <iostream>
#include <fstream>

int main(int argc, char** argv){

    std::string year;
    std::string observable;
    std::string process;
    std::string wilson;

    if(argc < 4 or argc >5){
        Log("Wrong number of arguments !!!");
        Log(" missing arguments : observable, year, combine process (OneBin, Unolled), [wilson] (fXX_L, fXZ_C, ...");
        return 0;
    }
    else{
         observable = argv[1];
         year       = argv[2];
         process    = argv[3];
         if(argc == 5)
            wilson = argv[4];
    }

    int nBin = 24;
 //   std::cout << "number of bin : ";
 //   std::cin >> nBin;

// ------------------------------------------------------------- //

    if(process == "OneBin"){

        std::vector<double> numberOfEvents;
        double number = 0;
        std::ifstream f("./combine/"+year+"/"+observable+"_noe_data_timed.txt");
        for(int i = 0; i < nBin; ++i){
            f >> number;
            numberOfEvents.push_back(number);
        }

        for(int i = 0; i < nBin; ++i){
            Card datacard;
            std::string name = observable+"_"+std::to_string(nBin)+"_"+std::to_string(i);

            datacard.addGlobalParameter(ttbarList);
            datacard.addSeparator();
            datacard.addInputsProcess("./inputs/"+year+"/", name+".root");
            datacard.addSeparator();
            datacard.addChanels(observable, numberOfEvents[i]);
            datacard.addSeparator();

            datacard.addProcToCard(observable, ttbarList);
            datacard.addSeparator();
            datacard.addRateToCard(ttbarList, systematicRate);
            for(std::string const& syst : systematicList)
                datacard.addSystToCard(syst, "shape", ttbarList);
            //datacard.addSystToCard("lumi", "lnN", groupList_p, "1.023");
            for(std::string const& syst : systematicTimeList)
                datacard.addSystToCard(syst, "shape", ttbarList);
            datacard.addSystToCard_alternative();
            datacard.addSeparator();
            //datacard.addLine("* autoMCStats 5");

            datacard.saveCard("./combine/"+year+"/one_bin/inputs/"+name+"_datacard.txt");
            if(i == 7){
                std::cout << " >> Card " << i << std::endl;
                datacard.printCard();
            }
        }
    }
    else if(process == "Unrolled"){

        double numberOfEvents;
        std::ifstream f("./combine/"+year+"/"+observable+"_noe_data.txt");
        f >> numberOfEvents;
        std::cout << "Number of events in data : " << numberOfEvents << std::endl;

        Card datacard;        
        std::string name = observable;

        datacard.addGlobalParameter(ttbarList);
        datacard.addSeparator();
        datacard.addInputsProcess("./inputs/"+year+"/", name+".root");
        datacard.addSeparator();
        datacard.addChanels(observable, numberOfEvents);
        datacard.addSeparator();

        datacard.addProcToCard(observable, ttbarList);
        datacard.addSeparator();
        datacard.addRateToCard(ttbarList, systematicRate);
        for(std::string const& syst : systematicList)
            datacard.addSystToCard(syst, "shape", ttbarList);
        //datacard.addSystToCard("lumi", "lnN", ttbarList, "1.023");
        for(std::string const& syst : systematicTimeList)
            datacard.addSystToCard(syst, "shape", ttbarList);
        datacard.saveCard("./combine/"+year+"/unrolled/inputs/"+name+"_datacard.txt");
    }

    else if(process == "SME"){

        ttbarList.push_back("");
        systematicRate.push_back("");
        for(size_t i = ttbarList.size(); i > 1 ; --i){
            ttbarList[i-1] = ttbarList[i-2];
        }
        for(size_t i = systematicRate.size(); i > 1 ; --i){
            systematicRate[i-1] = systematicRate[i-2];
        }
        ttbarList[0] =  wilson;
        systematicRate[0] =  "1.0";

        double numberOfEvents;
        std::ifstream f("./combine/"+year+"/"+observable+"_noe_data.txt");
        f >> numberOfEvents;
        std::cout << "Number of events in data : " << numberOfEvents << std::endl;

        Card datacard;        
        std::string name = observable;
        datacard.addGlobalParameter(ttbarList);
        datacard.addSeparator();
        datacard.addInputsProcess("./inputs/"+year+"/", name+"_"+wilson+".root");
        datacard.addSeparator();
        datacard.addChanels(observable, numberOfEvents);
        datacard.addSeparator();

        datacard.addProcToCard(observable, ttbarList);
        datacard.addSeparator();
        datacard.addRateToCard(ttbarList, systematicRate);
        for(std::string const& syst : systematicList)
            datacard.addSystToCard(syst, "shape", ttbarList);
        //datacard.addSystToCard("lumi", "lnN", ttbarList, "1.023");
        for(std::string const& syst : systematicTimeList)
            datacard.addSystToCard(syst, "shape", ttbarList);
        datacard.addSystToCard_alternative(true);
        datacard.addSeparator();
        //datacard.addLine("* autoMCStats 5");

        datacard.saveCard("./combine/"+year+"/sme/inputs/"+name+"_"+wilson+"_datacard.txt");
    }



    std::cout << "Finished !!!" << std::endl;
    return 0;
}

