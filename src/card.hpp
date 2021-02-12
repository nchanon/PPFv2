#pragma once 

#include <vector>
#include <string>

using  namelist = std::vector<std::string>;

class Card
{

    private:

        std::string completeBlock(std::string  const& word_p,
                                  unsigned int        blockSize_p
                                 );

        void addProcToCard(std::string const& observable_p,
                           namelist    const& groupList_p,
                           std::string      & card_p
                          );

        void addSystToCard(std::string const& systName_p,
                           std::string const& shape_p,
                           namelist    const& groupList_p,
                           std::string      & card_p,
                           std::string const& value_p = "1"
                          );

        void addRateToCard(namelist    const& groupList_p,
                           namelist    const& systematicsRate_p,
                           std::string      & card_p
                          );

        std::string cardGenerator(std::string const& rootfile_p,
                                  std::string const& observable_p,
                                  double             numberOfEvents_p,
                                  namelist    const& groupList_p,
                                  namelist    const& systematicsList_p,
                                  namelist    const& systematicsRate_p
                                 );


    public:


        Card();  


    void generateCard(std::string const& name_p,
                      std::string const& rootfile_p,
                      std::string const& observable_p,
                      double             numberOfEvents_p,
                      namelist    const& groupList_p,
                      namelist    const& systematicsList_p,
                      namelist    const& systematicsRate_p
                     );

};