#pragma once 

#include <string>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>


using  namelist = std::vector<std::string>;


class Generator{

    private:

        std::string observable;
        int nBin, minBin, maxBin;
        std::string year;

        std::vector<time_t> timestamp;
        std::vector<double> instLumi;

        time_t utcConverter(std::string const& time);

        double siderealHour(double time_p);

	double luminositySumOfWeight(TTree *tree_p);
        double luminosityCorrection(TTree *tree_p, double lumiavg);

	void drawHisto1D(TTree* tree, std::string obs, std::string string_eventSelection, std::string string_weight, std::string string_triggered, TH1F* hist);
        void drawHisto2D(TTree* tree, std::string obs1, std::string obs2, std::string string_eventSelection, std::string string_weight, std::string string_triggered, TH2F* hist);

        double generateWeight(TTree *tree_p, bool isTimed=true);
        std::string generateWeightString(bool isTimed=true);

        double generateSystematics(TTree            * tree_p,
                                   std::string const& systematicName,
                                   bool               isUp
                                  );
	std::string generateSystematicsString(std::string const& systematicName,
		                              bool               isUp
                                  	     );

	void generateLHEweightSystematics(TTree            * tree_p,
                                          std::string const& systematicName,
                                          double           * systList);
	void generateLHEweightSystematicsStrings(std::string const& systematicName,
                                                    std::string        * systList);

        void generateTimeSystematics(std::vector<double>      & weightsUp,
                                     std::vector<double>      & weightsDown
                                    );

        
        bool isTriggerPassed(TTree         * tree_p,
                             namelist const& triggerList_p,
                             bool            is2016H = false
                            );
	std::string isTriggerPassedString(namelist const& triggerList_p,
                             bool            is2016H = false
                            );


	bool eventSelection(TTree           * tree_p);
        bool eventSelection(TTree           * tree_p,
                                std::string     jecName);
	std::string eventSelectionString();
	std::string eventSelectionString(std::string     jecName);
 
        float getGenObservableValue(TTree         * tree_p); 
	float getObservableValue(TTree         * tree_p);
 
        void write(std::string       const& filename,
                   std::vector<TH1F>      & listObject,
                   std::string       const& option_p
                  );

        void write(std::string       const& filename,
                   std::vector<TH2F>      & listObject,
                   std::string       const& option_p
                  );


        void write(std::string       const& filename,
                   std::vector<std::vector<TH1F>> & listObject,
                   std::string       const& option_p
                  );


        void groupingMC(std::vector<TH1F>      & list,
                        namelist          const& groupList_p,
                        bool                     clean
                       );

        void groupingMC(std::vector<TH2F>      & list,
                        namelist          const& groupList_p,
                        std::string       const& name,
                        bool                     clean
                       );

        void groupingMC(std::vector<TH1F>      & list,
                        namelist          const& groupList_p,
                        std::string       const& name,
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

        void groupingSystematics(std::vector<TH2F>      & list,
                                 namelist          const& groupList_p,
				 std::string       const& name,
                                 namelist          const& systematicsList_p,
                                 bool                     isUp,
                                 bool                     clean
                                );


	void groupingLHEweightSystematics(std::vector<TH1F>      & listLHE,
                                          std::vector<TH1F>      & list,
                                    	  namelist          const& groupList_p,
                                          std::string           syst,
                                    	  bool                     clean
                                   	);
        
        void groupingLHEweightSystematics(std::vector<TH2F>      & listLHE,
                                          std::vector<TH2F>      & list,
                                          namelist          const& groupList_p,
 	                                  std::string       const& name,
                                          std::string           syst,
                                          bool                     clean
                                        );


    public:

	bool doLoop = false;

        Generator(std::string     const& observable_p,
                  std::vector<int> const& binning_p,
                  std::string     const& year_p
                 );

        ~Generator(){};

        void generateAltMC(namelist            const& sampleList_p,
                           namelist            const& groupList_p,
                           namelist            const& triggerList_p,
                           std::vector<double> const& correction_p,
                           bool                clean_p = true,
                           bool                isTimed_p = true
                          );

        void generateJecMC(namelist            const& sampleList_p,
                           namelist            const& jecList_p,
                           namelist            const& groupList_p,
                           namelist            const& triggerList_p,
                           //std::vector<std::vector<double>> const& correction_p,
                           std::vector<double> const& correction_p,
                           bool                clean_p = true,
                           bool                isTimed_p = true
                          );

        void generateMC(namelist            const& sampleList_p,
                        namelist            const& triggerList_p,
                        namelist            const& groupList_p,
                        namelist            const& systematicsList_p,
                        namelist            const& systematicsTimeList_p,
			std::vector<double>    const& numberofevents_p,
                        std::vector<double> const& correction_p,
                        std::string         const& rootOption_p,
                        bool                       clean_p = true,
                        bool                       isTimed_p = true

                       );

        void generateMCforComp(namelist            const& sampleList_p,
                        namelist            const& triggerList_p,
                        namelist            const& groupList_p,
                        std::vector<double> const& correction_p,
                        std::string         const& rootOption_p,
                        bool                       clean_p = true
                       );

        void generateData(namelist            const& sampleList_p,
                          namelist            const& triggerList_p,
                          namelist            const& groupList_p,
                          std::vector<double> const& correction_p,
                          std::string         const& rootOption_p,
                          bool                       correctedLumi,
                          bool                       clean_p = true
                         );

        void generateDataTimed(namelist            const& sampleList_p,
                               namelist            const& triggerList_p,
                               namelist            const& groupList_p,
                               std::vector<double> const& correction_p,
                               int                        nBin_p,
                               bool                       correctedLumi,
                               bool                       clean_p = true
                              );

};
