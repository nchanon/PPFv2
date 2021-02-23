#pragma once 

#include <vector>
#include <string>

using  namelist = std::vector<std::string>;

class Card
{

    private:

        std::string datacard;

        std::string completeBlock(std::string  const& word_p,
                                  unsigned int        blockSize_p
                                 );


    public:


        Card();  

        void addSeparator();

        void addGlobalParameter(namelist const& groupList_p,
                                int             numberOfBins = 1
                               );

        void addInputsProcess(std::string const& directory_p, 
                              std::string const& rootfile_p
                             );

        void addChanels(std::string const& observable_p,
                        double             numberOfEvents_p
                       );

        void addProcToCard(std::string const& observable_p,
                           namelist    const& groupList_p
                          );

        void addSystToCard(std::string const& systName_p,
                           std::string const& shape_p,
                           namelist    const& groupList_p,
                           std::string const& value_p = "1"
                          );

        void addRateToCard(namelist    const& groupList_p,
                           namelist    const& systematicsRate_p
                          );


        void printCard();
        void saveCard(std::string const& name_p);

};