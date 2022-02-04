#include "../src/card.hpp"
#include "../src/debug.h"
#include "../src/sample_2017.hpp"
#include "../src/sample_2016.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include <TH1F.h>
#include <TFile.h>

using namespace std;

int main(int argc, char** argv){

    std::string year;
    std::string observable;
    std::string process;
    std::string wilson;

    std::string debug_syst = "...";

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

    int iyear = -1;
    if (year=="2016") iyear = 0;
    if (year=="2017") iyear = 1;

    //Lumi flat uncertainties
    std::string lumi_flat_uncorr[2] = {"1.009", "1.014"};
    std::string lumi_flat_corr[2] = {"1.006", "1.009"};
    std::string lumi_flat_uncorr_inc[2] = {"1.01", "1.02"};
    std::string lumi_flat_corr_inc[2] = {"1.006", "1.009"};

    //Alternative theory sample uncertainties
    std::string alt_hdamp[2] = {"1", "1"};
    std::string alt_CP5[2] =  {"1", "1"};
    std::string alt_erdOn[2] = {"1", "1"};
    std::string alt_GluonMove[2] = {"1", "1"};
    std::string alt_QCDinspired[2] = {"1", "1"};

    std::string sinc="";
    if (process == "Inclusive") sinc = "_inclusive";

    string f_nom_name = "./results/" + year + "/flattree/" + observable + sinc + ".root";
    TFile* f_nom = new TFile(f_nom_name.c_str(), "READ");
    string f_colorreco_name = "./results/" + year + "/flattree/" + observable + "_color_reco"+sinc+".root";
    TFile* f_colorreco = new TFile(f_colorreco_name.c_str(), "READ");
    TH1F* h_nom = (TH1F*) f_nom->Get("signal");
    TH1F* h_hdampUp = (TH1F*) f_colorreco->Get("signal_hdampUp");
    TH1F* h_hdampDown = (TH1F*) f_colorreco->Get("signal_hdampDown");
    TH1F* h_CP5Up = (TH1F*) f_colorreco->Get("signal_CP5Up");
    TH1F* h_CP5Down = (TH1F*) f_colorreco->Get("signal_CP5Down");
    TH1F* h_erdOnUp = (TH1F*) f_colorreco->Get("signal_erdOnUp");
    TH1F* h_GluonMoveUp = (TH1F*) f_colorreco->Get("signal_GluonMoveUp");
    TH1F* h_QCDinspiredUp = (TH1F*) f_colorreco->Get("signal_QCDinspiredUp");

    string tmp1 = std::to_string(h_hdampUp->GetBinContent(1)/h_nom->GetBinContent(1));
    string tmp2 = std::to_string(h_hdampDown->GetBinContent(1)/h_nom->GetBinContent(1));
    alt_hdamp[iyear] = tmp1.substr(0, tmp1.find(".")+4) + "/" + tmp2.substr(0, tmp2.find(".")+4);
    tmp1 = std::to_string(h_CP5Up->GetBinContent(1)/h_nom->GetBinContent(1));
    tmp2 = std::to_string(h_CP5Down->GetBinContent(1)/h_nom->GetBinContent(1));
    alt_CP5[iyear] = tmp1.substr(0, tmp1.find(".")+4) + "/" + tmp2.substr(0, tmp2.find(".")+4);
    tmp1 = std::to_string(h_erdOnUp->GetBinContent(1)/h_nom->GetBinContent(1));
    alt_erdOn[iyear] =  tmp1.substr(0, tmp1.find(".")+4);
    tmp1 = std::to_string(h_GluonMoveUp->GetBinContent(1)/h_nom->GetBinContent(1));
    alt_GluonMove[iyear] =  tmp1.substr(0, tmp1.find(".")+4);
    tmp1 = std::to_string(h_QCDinspiredUp->GetBinContent(1)/h_nom->GetBinContent(1));
    alt_QCDinspired[iyear] = tmp1.substr(0, tmp1.find(".")+4);



    //std::cout << alt_hdamp[iyear] << std::endl; 

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
            datacard.addRateToCard(ttbarList, systematicRate,false);
            for(std::string const& syst : systematicList){
                if(syst == "syst_pt_top")
                    datacard.addProcSystToCard(syst, "shape", ttbarList, "signal",false);
                else if(syst == debug_syst)
                    continue;
                else
                    datacard.addSystToCard(syst, "shape", ttbarList);
            }
            datacard.addSystToCard("lumi_flat_uncor_"+year, "lnN", ttbarList, lumi_flat_uncorr[iyear]);
            datacard.addSystToCard("lumi_flat_cor", "lnN", ttbarList, lumi_flat_corr[iyear]);
            datacard.addProcSystToCard("hdamp", "lnN", ttbarList, "signal", false, alt_hdamp[iyear]);
            datacard.addProcSystToCard("CP5", "lnN", ttbarList, "signal", false, alt_CP5[iyear]);
            datacard.addProcSystToCard("erdOn", "lnN", ttbarList, "signal", false, alt_erdOn[iyear]);
            datacard.addProcSystToCard("GluonMove", "lnN", ttbarList, "signal", false, alt_GluonMove[iyear]);
            datacard.addProcSystToCard("QCDinspired", "lnN", ttbarList, "signal", false, alt_QCDinspired[iyear]);

            for(std::string const& syst : systematicTimeList){
                if (syst == "lumi_flat") continue;
                if (syst != "emu_trig" && syst != "lumi_stability" && syst != "lumi_linearity") datacard.addSystToCard(syst, "shape", ttbarList);
                else datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
            }

            //datacard.addSystToCard("lumi", "lnN", groupList_p, "1.023");
            //for(std::string const& syst : systematicTimeList){
            //    if(syst == debug_syst)
            //        continue;
            //    else
            //        datacard.addSystToCard(syst, "shape", ttbarList);
            //}
            datacard.addSystToCard_alternative(false);
            datacard.addSeparator();
            datacard.addLine("* autoMCStats 0");

            datacard.saveCard("./combine/"+year+"/one_bin/inputs/"+name+"_datacard.txt");
            if(i == 7){
                std::cout << " >> Card " << i << std::endl;
                datacard.printCard();
            }
        }
    }
    else if(process == "Inclusive"){

        double numberOfEvents;
        std::ifstream f("./combine/"+year+"/"+observable+"_noe_data.txt");
        f >> numberOfEvents;
        std::cout << "Number of events in data : " << numberOfEvents << std::endl;

        Card datacard;
        std::string name = observable;

        datacard.addGlobalParameter(ttbarList);
        datacard.addSeparator();
        datacard.addInputsProcess("./inputs/"+year+"/", name+"_inclusive.root");
        datacard.addSeparator();
        datacard.addChanels(observable, numberOfEvents);
        datacard.addSeparator();

        datacard.addProcToCard(observable, ttbarList);
        datacard.addSeparator();
        datacard.addRateToCard(ttbarList, systematicRate,false);
        for(std::string const& syst : systematicList){
            if(syst == "syst_pt_top")
                datacard.addProcSystToCard(syst, "shape", ttbarList, "signal",false);
            else if (syst == "syst_em_trig")
                datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
            else
                datacard.addSystToCard(syst, "shape", ttbarList);
        }

        datacard.addSystToCard("lumi_uncor_"+year, "lnN", ttbarList, lumi_flat_uncorr_inc[iyear]);
        datacard.addSystToCard("lumi_cor", "lnN", ttbarList, lumi_flat_corr_inc[iyear]);
        datacard.addProcSystToCard("hdamp", "lnN", ttbarList, "signal", false, alt_hdamp[iyear]);
        datacard.addProcSystToCard("CP5", "lnN", ttbarList, "signal", false, alt_CP5[iyear]);
        datacard.addProcSystToCard("erdOn", "lnN", ttbarList, "signal", false, alt_erdOn[iyear]);
        datacard.addProcSystToCard("GluonMove", "lnN", ttbarList, "signal", false, alt_GluonMove[iyear]);
        datacard.addProcSystToCard("QCDinspired", "lnN", ttbarList, "signal", false, alt_QCDinspired[iyear]);

        //for(std::string const& syst : systematicTimeList){
        //    if (syst == "lumi_flat") continue;
        //    if (syst != "emu_trig" && syst != "lumi_stability" && syst != "lumi_linearity") datacard.addSystToCard(syst, "shape", ttbarList);
        //    else datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
        //}
        datacard.addSystToCard_alternative(false);
        datacard.addSeparator();
        datacard.addLine("* autoMCStats 0");

        datacard.saveCard("./combine/"+year+"/inclusive/inputs/"+name+"_datacard.txt");

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
        datacard.addRateToCard(ttbarList, systematicRate,false);
        for(std::string const& syst : systematicList){
            if(syst == "syst_pt_top")
                datacard.addProcSystToCard(syst, "shape", ttbarList, "signal",false);
            else if (syst == "syst_em_trig") continue;
            else
                datacard.addSystToCard(syst, "shape", ttbarList);
        }
        //datacard.addSystToCard("lumi", "lnN", ttbarList, "1.023");
        for(std::string const& syst : systematicTimeList)
            datacard.addSystToCard(syst, "shape", ttbarList);
	datacard.addSystToCard_alternative(false);

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
        systematicRate[0] =  "1.05";

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
        datacard.addRateToCard(ttbarList, systematicRate,true);
        for(std::string const& syst : systematicList){
            if(syst == "syst_pt_top")
                datacard.addProcSystToCard(syst, "shape", ttbarList, "signal",true);
            else if (syst == "syst_em_trig") continue;
            else
                datacard.addSystToCard(syst, "shape", ttbarList);
        }
        //datacard.addSystToCard("lumi", "lnN", ttbarList, "1.023");
        datacard.addSystToCard("lumi_flat_uncor_"+year, "lnN", ttbarList, lumi_flat_uncorr[iyear]);
        datacard.addSystToCard("lumi_flat_cor", "lnN", ttbarList, lumi_flat_corr[iyear]);
	datacard.addProcSystToCard("hdamp", "lnN", ttbarList, "signal", true, alt_hdamp[iyear]);
        datacard.addProcSystToCard("CP5", "lnN", ttbarList, "signal", true, alt_CP5[iyear]);
        datacard.addProcSystToCard("erdOn", "lnN", ttbarList, "signal", true, alt_erdOn[iyear]);
        datacard.addProcSystToCard("GluonMove", "lnN", ttbarList, "signal", true, alt_GluonMove[iyear]);
        datacard.addProcSystToCard("QCDinspired", "lnN", ttbarList, "signal", true, alt_QCDinspired[iyear]);

        for(std::string const& syst : systematicTimeList){
            if (syst == "lumi_flat") continue;
	    if (syst != "emu_trig" && syst != "lumi_stability" && syst != "lumi_linearity") datacard.addSystToCard(syst, "shape", ttbarList);
	    else datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
	}
        datacard.addSystToCard_alternative(true);
        datacard.addSeparator();
        datacard.addLine("* autoMCStats 0");

        datacard.saveCard("./combine/"+year+"/sme/inputs/"+name+"_"+wilson+"_datacard.txt");
    }



    std::cout << "Finished !!!" << std::endl;
    return 0;
}

