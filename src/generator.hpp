#pragma once 

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

using  namelist = std::vector<std::string>;


class Generator{

    private:

        std::string observable;
        int nBin, minBin, maxBin;
        std::string year;

        double siderealHour(double time_p);

        double generateWeight(TTree *tree_p);

        double generateSystematics(TTree            * tree_p,
                                   std::string const& systematicName,
                                   bool               isUp
                                  );

        void generateTimeSystematics(std::vector<double>      & weightsUp,
                                     std::vector<double>      & weightsDown
                                    );

        
        bool isTriggerPassed(TTree         * tree_p,
                             namelist const& triggerList_p
                            );
        
        void write(std::string       const& filename,
                   std::vector<TH1F>      & listObject,
                   std::string       const& option_p
                  );

        void groupingMC(std::vector<TH1F>      & list,
                        namelist          const& groupList_p,
                        bool                     clean
                       );

        void groupingData(std::vector<TH1F>      & list,
                          namelist          const& groupList_p,                        bool                     clean
                         );


        void groupingDataTimed(std::vector<TH1F>      & list,
                               namelist          const& groupList_p,
                               int                      bin,
                               int                      nBin,
                               bool                     clean
                              );


        void groupingSystematics(std::vector<TH1F>      & list,
                                 namelist          const& groupList_p,
                                 namelist          const& systematicsList_p,
                                 bool                     isUp,
                                 bool                     clean
                                );

        


    public:

        Generator(std::string     const& observable_p,
                  std::vector<int> const& binning_p,
                  std::string     const& year_p
                 );

        ~Generator(){};

        void generateMC(namelist            const& sampleList_p,
                        namelist            const& triggerList_p,
                        namelist            const& groupList_p,
                        namelist            const& systematicsList_p,
                        namelist            const& systematicsTimeList_p,
                        std::vector<double> const& correction_p,
                        std::string         const& rootOption_p,
                        bool                       clean_p = true
                       );

        void generateData(namelist            const& sampleList_p,
                          namelist            const& triggerList_p,
                          namelist            const& groupList_p,
                          std::vector<double> const& correction_p,
                          std::string         const& rootOption_p,
                          bool                       clean_p = true
                         );

        void generateDataTimed(namelist            const& sampleList_p,
                               namelist            const& triggerList_p,
                               namelist            const& groupList_p,
                               std::vector<double> const& correction_p,
                               int                        nBin_p,
                               bool                       clean_p = true
                              );

};