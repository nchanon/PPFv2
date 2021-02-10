#pragma once 

#include <vector>
#include <string>

using  namelist = std::vector<std::string>;

class Combine
{

    private:



        std::string cardGenerator(std::string const& rootfile_p,
                                  std::string const& observable_p,
                                  double             numberOfEvents_p,
                                  namelist    const& groupList_p
                                 );

    public:


        Combine();  





};