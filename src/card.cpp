#include "card.hpp"
#include "debug.h"

#include <fstream>

/////////////////////////////////
// Constructor
/////////////////////////////////

Card::Card(){}

/////////////////////////////////
// Private methods
/////////////////////////////////

std::string Card::completeBlock(std::string  const& word_p,
                                   unsigned int        blockSize_p
                                  )
{ 
    unsigned int size = blockSize_p - word_p.size();
    if(size == 0){
        Log("Error with "+word_p+" block in datacard");
        return "";
    }
    std::string spaces;
    for(size_t i = 0; i < size; ++i){
        spaces += " ";
    }
    return word_p+spaces;
}


void Card::addProcToCard(std::string const& observable_p,
                            namelist    const& groupList_p,
                            std::string      & card_p
                           )
{
    std::string line0 = completeBlock("bin", 30);
    std::string line1 = completeBlock("process", 30);
    std::string line2 = line1;
    std::string line3 = completeBlock("rate", 30);

    for(size_t i = 0; i < groupList_p.size(); ++i){
        line0 += completeBlock(observable_p, 12);
        line1 += completeBlock(groupList_p[i], 12);
        line2 += completeBlock(std::to_string(i), 12);
        line3 += completeBlock("-1", 12);
    }
    card_p += line0+'\n'+line1+'\n'+line2+'\n'+line3+'\n'; 
}

void Card::addSystToCard(std::string const& systName_p,
                            std::string const& shape_p,
                            namelist    const& groupList_p,
                            std::string      & card_p,
                            std::string const& value_p
                           )
{
    card_p += completeBlock(systName_p, 24) + completeBlock(shape_p, 6);
    for(size_t i = 0; i < groupList_p.size(); ++i)
        card_p += completeBlock(value_p, 12);
    card_p += '\n';
}

std::string Card::cardGenerator(std::string const& rootfile_p,
                                   std::string const& observable_p,
                                   double             numberOfEvents_p,
                                   namelist    const& groupList_p,
                                   namelist    const& systematicsList_p
                                  )
{
    std::string output;
    output += "imax 1 number of bins\n";
    output += "jmax "+std::to_string(groupList_p.size()-1)+" number of background processes\n";
    output += "kmax * number of nuisance parameters\n";
    output += "-------------------------------------------------------------\n";

    output += "shapes * * ./inputs/"+rootfile_p+" $PROCESS $PROCESS_$SYSTEMATIC\n";
    output += "shapes sig * ./inputs/"+rootfile_p+" $PROCESS $PROCESS_$SYSTEMATIC\n";
    output += "--------------------------------------------------------------\n";

    output += "bin "+observable_p+"\n";
    output += "observation "+std::to_string(numberOfEvents_p)+"\n";
    output += "--------------------------------------------------------------------------------------------- \n";

    addProcToCard(observable_p, groupList_p, output);
    output += "--------------------------------------------------------------------------------------------- \n";

    for(std::string const& syst : systematicsList_p){
        addSystToCard(syst, "shape", groupList_p, output);
    }
    addSystToCard("lumi", "lnN", groupList_p, output, "1.023");

    return output;
}


/////////////////////////////////
// Public methods
/////////////////////////////////


void Card::printCard(std::string const& rootfile_p,
                        std::string const& observable_p,
                        double             numberOfEvents_p,
                        namelist    const& groupList_p,
                        namelist    const& systematicsList_p
                       )
{
    std::cout << cardGenerator(rootfile_p, 
                               observable_p, 
                               numberOfEvents_p, 
                               groupList_p,
                               systematicsList_p
                              ); 
    std::cout << std::endl;
}

void Card::generateCard(std::string const& name_p,
                           std::string const& rootfile_p,
                           std::string const& observable_p,
                           double             numberOfEvents_p,
                           namelist    const& groupList_p,
                           namelist    const& systematicsList_p
                          )
{

    std::ofstream file(name_p);
    file << cardGenerator(rootfile_p, 
                          observable_p, 
                          numberOfEvents_p, 
                          groupList_p,
                          systematicsList_p
                         ); 
    file.close();
}