
##   PyPlotFramework

Current version 2.0

This framework is plotter based on pyROOT-CERN library. It was build for ttbar analysis.
PyPlotFramework is optimized for the CMS HEPPY output.

This framework is built with python2.7 !

# Installation :

In your terminal type following commands : 

    > git clone https://github.com/Arc-Pintade/PPFv2.git

In the pyplotframework directory :

    # will create the directory tree 
    > bash scripts/install.sh

    # will generate all the sample utils from tools/sample_manager.py in a C++ header file in src directory. 
    NB : a compilation will be launched.
    > python scripts/sample_update.py

    # this file need to be update with your path to sourcing ROOT
    > source rootenv

NB: Caution, you need to have ROOT-CERN lirbary on your computer. And the rootenv 
need to be update (if needed) in function of your own ROOT installation path.


# Presentation of the framework 

This PyPlotFramework, is built in 2 part: 

  -> Build histograms from Heppy output rootfile, with reweight for MC samples and 
  triggers for Data. 

  These histograms are stored in :

    ./results/"year"/flattree

  -> Creation of comparaison data/mc plots, ...

  These results are stored in :

    ./results/"year"/comparaison
    ./results/"year"/systematics
    ...

  -> Creation of combine inputs and datacards, ...
    
    ./combine/"year"/inclusive
    ./combine/"year"/unrolled
    ...

# Generate histograms, reweighted and triggered

This histogram creator will be executed for each observable. This c++ code will ask you a option to run on all type of sample or just specifique one.

The list of option is : [All, mc, data, jec, alt, timed, sme]
'sme' create a rootfile with f function in './results/"year"/flattree' directory, usefull for combine inputs creation.

    > ./bin/histograms_creator "observable" "year" "option"

example : 

    > ./bin/histograms_creator m_dilep 2017 All

NB : The particular case of alternative samples is for now in progress.
To make the treatment of the flattree produced by 'histograms_creator' for color reconection for exemple, you need call the command :

    >  python bin/color_reco.py "observable" "year"

example :

    >  python bin/color_reco.py m_dilep 2017


If your not sure of what you need, type option "All".

# Comparaison Data/Monte-Carlo

Data/mc comparaison are created for a given observable

    > python ./bin/data_mc_comparaison.py "observable" "year" "title"

example : 

    > python ./bin/data_mc_comparaison.py m_dilep 2017 "Mass dileptonic"

# Check Systematics

In PPFv2 there is 2 scripts to check systematics, the first is to get an particular systematics, the other to get them all in one time. The results are stored in './results/'year'/systematics'

It's possible to check each systematics : 

    > python ./bin/systematics_observable.py "observable" "year" "systematic" "title"

example : 

    > python ./bin/systematics_observable.py m_dilep 2017 syst_elec_id "electron id"

In case you don't sure what systematics are implemented just type :

    > python ./bin/systematics_observable.py

If you want you can run for all systematics in one time with the scripts :

    > python scripts/control_systematics.py

You can change the content of this scripts as you wish to run with other observable.

##

In the particular case of alternative sample systematics and jec systematics, the scripts is work in progress but you can draw plots with : 

    > python bin/check_systematics.py "observable" "year" "systematics"

example: 

    > python bin/check_systematics.py m_dilep 2017 jec

The result are stored in './results/'year'/other directory'. For now the accessible systematics are : 'jec', 'hdamp', 'CP5'

# Combine input file creation 

For combine work, you will need to create combine inputs and datacard separatly.
All will be stored in the ./combine/year/ directory.

    # combine inputs file
    > python ./bin/combine_"methode".py "observable" "year" 
    # datacards
    > ./bin/card_creator "observable" "year" " "methode"

example :

    > python ./bin/combine_unrolled.py m_dilep 2017 
    > ./bin/card_creator m_dilep 2017 Unrolled

    or 

    > python ./bin/combine_one_bin.py n_bjets 2016 
    > ./bin/card_creator n_bjets 2016 OneBin

Some usefull scripts can be used in the self-named directory. 

example : 

    This command will create combine inputs then datacards then export in the directory of your servers (adress in ./scripts/export_combine.py)
    > bash scripts/launch_unrolled_stuff.sh m_dilep 2017

# Example of full analysis.

For example, to make an analysis with the dilepton mass observable for 2017 with cLXX Wilson coefficient, let's type : 

    # produce flattree correctred :
    > ./bin/histograms_creator m_dilep 2017 All
    >  python bin/color_reco.py m_dilep 2017

    # produce inputs for combine, differential part :
    > python ./bin/combine_one_bin.py m_dilep 2017 
    > ./bin/card_creator m_dilep 2017 OneBin
    
    # produce inputs for combine, cmunu measurement part :
    > python ./bin/combine_unrolled.py m_dilep 2017 
    > ./bin/card_creator m_dilep 2017 SME cLXX

