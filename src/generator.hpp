#pragma once 

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

using  namelist = std::vector<std::string>;


class Generator{

    private:

        std::string year;

        double generateWeight(TTree *tree_p);

        double generateSystematics(TTree            * tree_p,
                                   std::string const& systematicName,
                                   bool               isUp);

        
        bool isTriggerPassed(TTree         * tree_p,
                             namelist const& triggerList_p);
        
        void write(std::string       const& filename,
                   std::vector<TH1F>      & listObject,
                   std::string       const& option_p);

        void grouping(std::vector<TH1F>      & list,
                      namelist          const& groupList_p,
                      bool                     isData = false);

        void groupingSystematics(std::vector<TH1F>      & list,
                                 namelist          const& groupList_p,
                                 namelist          const& systematicsList_p,
                                 bool                     isUp);

        


    public:

        Generator(std::string const& year_p);
        ~Generator(){};

        void generate(std::string         const& nature_p, 
                      std::string         const& observable_p,
                      namelist            const& sampleList_p,
                      namelist            const& triggerList_p,
                      namelist            const& groupList_p,
                      namelist            const& systematicsList_p,
                      std::vector<double> const& correction_p,
                      std::string         const& option_p);


};