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

    bool doExpTimeNuisance = true;
    int triggerOption = 0; //Full trigger syst uncertainties including Nvtx partition
    //int triggerOption = 1; //Trigger syst uncertainties without Nvtx partition
    //int triggerOption = 2; //Trigger syst uncertainties treated as uncorrelated in time

    bool doPuTime = true;
    //bool doPuTime = false;
    
    bool doAllWilson = false;
    if (wilson=="sme_all") doAllWilson = true;
    std::string wilsonList[16] = {"cLXX","cLXY","cLXZ","cLYZ","cRXX","cRXY","cRXZ","cRYZ","cXX","cXY","cXZ","cYZ","dXX","dXY","dXZ","dYZ"};

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

    //std::string stimebin = "_puold";
    std::string stimebin = "_punew";
    //std::string stimebin = "_puinc";
    //if (timebin>=0) stimebin = "_put"+std::to_string(timebin);    

    string f_nom_name = "./results/" + year + "/flattree/" + observable + sinc + stimebin + ".root";
    //std::cout << "Histogram file: "<<f_nom_name<<std::endl;

    TFile* f_nom = new TFile(f_nom_name.c_str(), "READ");
    string f_colorreco_name = "./results/" + year + "/flattree/" + observable + "_color_reco"+sinc+stimebin+".root";
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

	string timebin="";

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

            if (doExpTimeNuisance){
		timebin = "_t" + std::to_string(i);
            }   

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
                else if (syst == "syst_em_trig") {
		    if (triggerOption==0 || triggerOption==1) continue;
		    else if (triggerOption==2) datacard.addSystToCard(syst + "_" + year + timebin, "shape", ttbarList);
		}
		else if (syst == "syst_b_uncorrelated" || syst == "syst_l_uncorrelated")
		    datacard.addSystToCard(syst + "_" + year + timebin, "shape", ttbarList);
		else if (syst == "syst_qcdscale"){
		    datacard.addProcSystToCard(syst + "_signal", "shape", ttbarList, "signal",false);
		    datacard.addProcSystToCard(syst + "_singletop", "shape", ttbarList, "singletop",false);
                    datacard.addProcSystToCard(syst + "_ttx", "shape", ttbarList, "ttx",false);
                    datacard.addProcSystToCard(syst + "_vjets", "shape", ttbarList, "vjets",false);
		}
                else if (syst == "syst_pdfas"){
		    datacard.addSystToCard("syst_pdfas", "shape", ttbarList, "1", "dibosons");
                }
		else if (syst == "syst_ps_isr"){
		     datacard.addProcSystToCard(syst + "_signal", "shape", ttbarList, "signal",false);
                     datacard.addProcSystToCard(syst + "_singletop", "shape", ttbarList, "singletop",false);
		}
		else if (syst == "syst_ps_fsr"){
		     datacard.addProcSystToCard(syst, "shape", ttbarList, "signal", false, "1", "singletop");
		}
                else if(syst == debug_syst)
                    continue;
		else if (syst == "syst_pu" && doPuTime){
		    datacard.addSystToCard(syst, "shape", ttbarList);
		}
                else
                    datacard.addSystToCard(syst + timebin, "shape", ttbarList);
            }
            datacard.addSystToCard("lumi_flat_uncor_"+year, "lnN", ttbarList, lumi_flat_uncorr[iyear]);
            datacard.addSystToCard("lumi_flat_cor", "lnN", ttbarList, lumi_flat_corr[iyear]);
            datacard.addProcSystToCard("hdamp", "lnN", ttbarList, "signal", false, alt_hdamp[iyear]);
            datacard.addProcSystToCard("CP5", "lnN", ttbarList, "signal", false, alt_CP5[iyear]);
            datacard.addProcSystToCard("erdOn", "lnN", ttbarList, "signal", false, alt_erdOn[iyear]);
            datacard.addProcSystToCard("GluonMove", "lnN", ttbarList, "signal", false, alt_GluonMove[iyear]);
            datacard.addProcSystToCard("QCDinspired", "lnN", ttbarList, "signal", false, alt_QCDinspired[iyear]);
	    datacard.addProcSystToCard("mtop", "shape", ttbarList, "signal",false);

  	    //datacard.addSystToCard("jec" + timebin, "shape", ttbarList);
            datacard.addSystToCard("Absolute_jec"+ timebin, "shape", ttbarList);
            datacard.addSystToCard("Absolute_"+year+"_jec"+ timebin, "shape", ttbarList);
            datacard.addSystToCard("FlavorQCD_jec"+ timebin, "shape", ttbarList);
            datacard.addSystToCard("BBEC1_jec"+ timebin, "shape", ttbarList);
            datacard.addSystToCard("BBEC1_"+year+"_jec"+ timebin, "shape", ttbarList);
            datacard.addSystToCard("RelativeBal_jec"+ timebin, "shape", ttbarList);
            datacard.addSystToCard("RelativeSample_"+year+"_jec"+ timebin, "shape", ttbarList);


            for(std::string const& syst : systematicTimeList){
                if (syst == "lumi_flat") continue;
                if (syst != "emu_trig" && syst != "lumi_stability" && syst != "lumi_linearity") datacard.addSystToCard(syst, "shape", ttbarList);
                else if (syst != "emu_trig") datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
		else if (syst == "emu_trig") {
		    if (triggerOption==0 || triggerOption==1) datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
		    if (triggerOption==2) continue;
		}
            }

            //datacard.addSystToCard_alternative(false);
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

        string f_data_name = "./results/" + year + "/flattree/" + observable + "_data.root";
        TFile* f_data = new TFile(f_data_name.c_str(), "READ");
	TH1F* h_data = (TH1F*) f_data->Get("data_obs");
        double numberOfEvents = h_data->Integral();
        //std::ifstream f("./combine/"+year+"/"+observable+"_noe_data.txt");
        //f >> numberOfEvents;
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
            else if (syst == "syst_b_uncorrelated" || syst == "syst_l_uncorrelated")
                datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
            else if (syst == "syst_qcdscale"){
                datacard.addProcSystToCard(syst + "_signal", "shape", ttbarList, "signal",false);
                datacard.addProcSystToCard(syst + "_singletop", "shape", ttbarList, "singletop",false);
                datacard.addProcSystToCard(syst + "_ttx", "shape", ttbarList, "ttx",false);
                datacard.addProcSystToCard(syst + "_vjets", "shape", ttbarList, "vjets",false);
            }
            else if (syst == "syst_pdfas"){
                datacard.addSystToCard("syst_pdfas", "shape", ttbarList, "1", "dibosons");
            }
            else if (syst == "syst_ps_isr"){
                datacard.addProcSystToCard(syst + "_signal", "shape", ttbarList, "signal",false);
                datacard.addProcSystToCard(syst + "_singletop", "shape", ttbarList, "singletop",false);
            }
            else if (syst == "syst_ps_fsr"){
                datacard.addProcSystToCard(syst, "shape", ttbarList, "signal", false, "1", "singletop");
            }
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
        datacard.addProcSystToCard("mtop", "shape", ttbarList, "signal",false);

        datacard.addSystToCard("Absolute_jec", "shape", ttbarList);
        datacard.addSystToCard("Absolute_"+year+"_jec", "shape", ttbarList);
        datacard.addSystToCard("FlavorQCD_jec", "shape", ttbarList);
        datacard.addSystToCard("BBEC1_jec", "shape", ttbarList);
        datacard.addSystToCard("BBEC1_"+year+"_jec", "shape", ttbarList);
        datacard.addSystToCard("RelativeBal_jec", "shape", ttbarList);
        datacard.addSystToCard("RelativeSample_"+year+"_jec", "shape", ttbarList);


        //datacard.addSystToCard_alternative(false);
        datacard.addSeparator();
        datacard.addLine("* autoMCStats 0");

        datacard.saveCard("./combine/"+year+"/inclusive/inputs/"+name+"_datacard.txt");
	datacard.printCard();

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

        string timebin[24]; //="";
        if (doExpTimeNuisance){
            for (int k=0; k<24; k++) timebin[k] = "_t" + std::to_string(k);
        }

	if (!doAllWilson){
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
	}
	else if (doAllWilson){
	    for (int iw=0; iw<16; iw++){
		ttbarList.insert(ttbarList.begin(), wilsonList[15-iw]);
		systematicRate.insert(systematicRate.begin(), "1.05");
	    }
	}

        double numberOfEvents;
        std::ifstream f("./combine/"+year+"/"+observable+"_noe_data.txt");
        f >> numberOfEvents;
        std::cout << "Number of events in data : " << numberOfEvents << std::endl;

        Card datacard;        
        std::string name = observable;
        datacard.addGlobalParameter(ttbarList);
        datacard.addSeparator();
        if (!doAllWilson) datacard.addInputsProcess("./inputs/"+year+"/", name+"_"+wilson+".root");
	if (doAllWilson) datacard.addInputsProcess("./inputs/"+year+"/", name+"_sme_all.root");
        datacard.addSeparator();
        datacard.addChanels(observable, numberOfEvents);
        datacard.addSeparator();

        datacard.addProcToCard(observable, ttbarList);
        datacard.addSeparator();
        datacard.addRateToCard(ttbarList, systematicRate,true);
        for(std::string const& syst : systematicList){
            if(syst == "syst_pt_top")
                datacard.addProcSystToCard(syst, "shape", ttbarList, "signal",true);
            else if (syst == "syst_em_trig") {
		if (triggerOption==0 || triggerOption==1) continue;
		if (triggerOption==2) {
		    for (int k=0; k<24; k++) datacard.addSystToCard(syst + "_" + year + timebin[k], "shape", ttbarList);
		}
	    }
            else if (syst == "syst_b_uncorrelated" || syst == "syst_l_uncorrelated") {
                if (doExpTimeNuisance) {
		    for (int k=0; k<24; k++) datacard.addSystToCard(syst + "_" + year + timebin[k], "shape", ttbarList);
		}
		else datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
	    }
            else if (syst == "syst_qcdscale"){
                datacard.addProcSystToCard(syst + "_signal", "shape", ttbarList, "signal",true);
                datacard.addProcSystToCard(syst + "_singletop", "shape", ttbarList, "singletop",false);
                datacard.addProcSystToCard(syst + "_ttx", "shape", ttbarList, "ttx",false);
                datacard.addProcSystToCard(syst + "_vjets", "shape", ttbarList, "vjets",false);
            }
            else if (syst == "syst_pdfas"){
                datacard.addSystToCard("syst_pdfas", "shape", ttbarList, "1", "dibosons");
            }
            else if (syst == "syst_ps_isr"){
                datacard.addProcSystToCard(syst + "_signal", "shape", ttbarList, "signal",true);
                datacard.addProcSystToCard(syst + "_singletop", "shape", ttbarList, "singletop",false);
            }
            else if (syst == "syst_ps_fsr"){
                datacard.addProcSystToCard(syst, "shape", ttbarList, "signal", false, "1", "singletop");
            }
            else {
		if (doExpTimeNuisance) {
		    for (int k=0; k<24; k++) datacard.addSystToCard(syst + timebin[k], "shape", ttbarList);
		}
                else datacard.addSystToCard(syst, "shape", ttbarList);
	    }
        }
        //datacard.addSystToCard("lumi", "lnN", ttbarList, "1.023");
        datacard.addSystToCard("lumi_flat_uncor_"+year, "lnN", ttbarList, lumi_flat_uncorr[iyear]);
        datacard.addSystToCard("lumi_flat_cor", "lnN", ttbarList, lumi_flat_corr[iyear]);
	datacard.addProcSystToCard("hdamp", "lnN", ttbarList, "signal", true, alt_hdamp[iyear]);
        datacard.addProcSystToCard("CP5", "lnN", ttbarList, "signal", true, alt_CP5[iyear]);
        datacard.addProcSystToCard("erdOn", "lnN", ttbarList, "signal", true, alt_erdOn[iyear]);
        datacard.addProcSystToCard("GluonMove", "lnN", ttbarList, "signal", true, alt_GluonMove[iyear]);
        datacard.addProcSystToCard("QCDinspired", "lnN", ttbarList, "signal", true, alt_QCDinspired[iyear]);
        datacard.addProcSystToCard("mtop", "shape", ttbarList, "signal",true);
	if (doExpTimeNuisance) {
            for (int k=0; k<24; k++) {
		//datacard.addSystToCard("jec" + timebin[k], "shape", ttbarList);
                datacard.addSystToCard("Absolute_jec"+ timebin[k], "shape", ttbarList);
                datacard.addSystToCard("Absolute_"+year+"_jec"+ timebin[k], "shape", ttbarList);
                datacard.addSystToCard("FlavorQCD_jec"+ timebin[k], "shape", ttbarList);
                datacard.addSystToCard("BBEC1_jec"+ timebin[k], "shape", ttbarList);
                datacard.addSystToCard("BBEC1_"+year+"_jec"+ timebin[k], "shape", ttbarList);
                datacard.addSystToCard("RelativeBal_jec"+ timebin[k], "shape", ttbarList);
                datacard.addSystToCard("RelativeSample_"+year+"_jec"+ timebin[k], "shape", ttbarList);
	    }
	}
        //else datacard.addSystToCard("jec", "shape", ttbarList);

        for(std::string const& syst : systematicTimeList){
            if (syst == "lumi_flat") continue;
	    if (syst != "emu_trig" && syst != "lumi_stability" && syst != "lumi_linearity") datacard.addSystToCard(syst, "shape", ttbarList);
	    else if (syst != "emu_trig") datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
	    else if (syst=="emu_trig") {
		if (triggerOption==0 || triggerOption==1) datacard.addSystToCard(syst + "_" + year, "shape", ttbarList);
		else if (triggerOption==2) continue;
	    }
	}

        datacard.addProcSystToCard("sme_decay", "shape", ttbarList, "singletop",false);

        //datacard.addSystToCard_alternative(true);
        datacard.addSeparator();
        datacard.addLine("* autoMCStats 0");

        datacard.saveCard("./combine/"+year+"/sme/inputs/"+name+"_"+wilson+"_datacard.txt");
    }



    std::cout << "Finished !!!" << std::endl;
    return 0;
}

