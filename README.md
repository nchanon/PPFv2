
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
    > bash install.sh

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

This histogram creator will be executed for each observable. This c++ code will ask you if you want to run on MC, time binned observable or on All (MC+data)

    > ./bin/histograms_creator "observable" "year"

example : 

    > ./bin/histograms_creator m_dilep 2017


# Comparaison Data/Monte-Carlo

Data/mc comparaison are created for a given observable

    > python ./bin/data_mc_comparaison.py "observable" "year" "title"

example : 

    > python ./bin/data_mc_comparaison.py m_dilep 2017 "Mass dileptonic"

# Combine input file creation 

For combine work, you will need to create combine inputs and datacard separatly.
All will be stored in the ./combine/year/ directory.

    # combine inputs file
    > python ./bin/combine_"methode".py "observable" "year" 
    # datacards
    > ./bin/card_creator "observable" "year" " "methode"

example :

    > python ./bin/combine_unrolled.py m_dilep 2017 
    > ./bin/card_creator.py m_dilep 2017 Unrolled

    or 

    > python ./bin/combine_one_bin.py n_bjets 2016 
    > ./bin/card_creator.py n_bjets 2016 OneBin

Some usefull scripts can be used in the self-named directory. 

example : 

    This command will create combine inputs then datacards then export in the directory of your servers (adress in ./scripts/export_combine.py)
    > bash scripts/launch_unrolled_stuff.sh m_dilep 2017
