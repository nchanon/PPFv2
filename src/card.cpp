#include "card.hpp"
#include "debug.h"

#include <fstream>

unsigned int block_proc = 30;
unsigned int block_syst = 24;
unsigned int block_grp  = 15;

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

/////////////////////////////////
// Public methods
/////////////////////////////////

void Card::addSeparator()
{
    datacard += "--------------------------------------------------------------------------------------------- \n";
}

void Card::addLine(std::string const& line)
{
    datacard += line + '\n';
}


void Card::addGlobalParameter(namelist const& groupList_p,
                              int             numberOfBins)
{
    datacard += "imax "+std::to_string(numberOfBins)+" number of bins\n";
    datacard += "jmax "+std::to_string(groupList_p.size()-1)+" number of background processes\n";
    datacard += "kmax * number of nuisance parameters\n";
}

void Card::addInputsProcess(std::string const& directory_p,
                            std::string const& rootfile_p
                           )
{
    datacard += "shapes * * "+directory_p+rootfile_p+" $PROCESS $PROCESS_$SYSTEMATIC\n";
    datacard += "shapes sig * "+directory_p+rootfile_p+" $PROCESS $PROCESS_$SYSTEMATIC\n";
}

void Card::addChanels(std::string const& observable_p,
                      double             numberOfEvents_p
                     )
{
    datacard += "bin "+observable_p+"\n";
    datacard += "observation "+std::to_string(numberOfEvents_p)+"\n";
}


void Card::addProcToCard(std::string const& observable_p,
                         namelist    const& groupList_p
                        )
{

    std::string line0 = completeBlock("bin", block_proc);
    std::string line1 = completeBlock("process", block_proc);
    std::string line2 = line1;
    std::string line3 = completeBlock("rate", block_proc);

    for(int i = 0; i < int(groupList_p.size()); ++i){
        line0 += completeBlock(observable_p, block_grp);
        line1 += completeBlock(groupList_p[i], block_grp);
        if(groupList_p[0] == "signal")
            line2 += completeBlock(std::to_string(i), block_grp);
        else
            line2 += completeBlock(std::to_string(i-1), block_grp);
        line3 += completeBlock("-1", block_grp);
    }
    datacard += line0+'\n'+line1+'\n'+line2+'\n'+line3+'\n'; 
}

void Card::addSystToCard_alternative(bool isSME)
{
    if(!isSME){
      datacard += "CP5                     shape 1              0              0              0              0              \n";
      datacard += "hdamp                   shape 1              0              0              0              0              \n";
      datacard += "erdOn                   shape 1              0              0              0              0              \n";
      datacard += "QCDinspired             shape 1              0              0              0              0              \n";
      datacard += "GluonMove               shape 1              0              0              0              0              \n";
      datacard += "mtop                    shape 1              0              0              0              0              \n";
      //datacard += "color_reco              shape 1              0              0              0              0              \n";
      datacard += "jec              	   shape 1              1              1              1              1              \n";
    }
    else{
      //datacard += "CP5                     shape - 1              -              -              -              -              \n";
      //datacard += "hdamp                   shape 0 1              -              -              -              -              \n";
      //datacard += "color_reco              shape 0 1              -              -              -              -              \n";
      datacard += "CP5                     shape 1              1              0              0              0              0              \n";
      datacard += "hdamp                   shape 1              1              0              0              0              0              \n";
      datacard += "erdOn                   shape 1              1              0              0              0              0              \n";
      datacard += "QCDinspired             shape 1              1              0              0              0              0              \n";
      datacard += "GluonMove               shape 1              1              0              0              0              0              \n";
      datacard += "mtop                    shape 1              1              0              0              0              0              \n";
      datacard += "jec              	  shape 1              1              1              1              1              1	\n";
    }

}

void Card::addProcSystToCard(std::string const& systName_p,
                             std::string const& shape_p,
                             namelist    const& groupList_p,
                             std::string const& process_p,
			     bool isSME
                             )
{
    datacard += completeBlock(systName_p, block_syst) 
            + completeBlock(shape_p, block_proc-block_syst);
    for(size_t i = 0; i < groupList_p.size(); ++i)
    {
	if (!isSME){
          if(groupList_p[i] == process_p)
              datacard += completeBlock("1", block_grp);
          else 
              datacard += completeBlock("0", block_grp);
	}
	else if (isSME){
          if(groupList_p[i] == process_p || i==0)
              datacard += completeBlock("1", block_grp);
          else 
              datacard += completeBlock("0", block_grp);
	}
    }
    datacard += '\n';
}


void Card::addSystToCard(std::string const& systName_p,
                            std::string const& shape_p,
                            namelist    const& groupList_p,
                            std::string const& value_p
                           )
{
    datacard += completeBlock(systName_p, block_syst) 
             + completeBlock(shape_p, block_proc-block_syst);
    for(size_t i = 0; i < groupList_p.size(); ++i)
        datacard += completeBlock(value_p, block_grp);
    datacard += '\n';
}

void Card::addRateToCard(namelist    const& groupList_p,
                         namelist    const& systematicsRate_p,
			 bool isSME
                        )
{
    for(size_t i = 1; i < groupList_p.size(); ++i){
	if (isSME && i==1) continue;
        std::string line = completeBlock("r"+groupList_p[i], block_syst) 
                         + completeBlock("lnN", block_proc-block_syst);
        for(size_t j = 0; j < groupList_p.size(); ++j){
            if(i == j)
                line += completeBlock(systematicsRate_p[j], block_grp);
	    else if (j == 0 && i==1 && isSME)
		line += completeBlock(systematicsRate_p[0], block_grp);
            else
                line += completeBlock("0", block_grp);
        }  
        datacard += line + '\n';
    }
}

void Card::printCard()
{
    std::cout << " >> datacard :" << std::endl;
    std::cout << datacard << std::endl;
}

void Card::saveCard(std::string const& name_p)
{
    std::ofstream file(name_p);
    file << datacard;
    file.close();
}
