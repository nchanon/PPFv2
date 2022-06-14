
##   PyPlotFramework

Current version 2.0

This framework is plotter based on pyROOT-CERN library. It was build for ttbar analysis.
PyPlotFramework is optimized for the CMS HEPPY output.

This framework is built with python2.7 !

# Installation :

In your terminal type following commands : 

    > git clone https://github.com/nchanon/PPFv2.git

In the pyplotframework directory :

    # will create the directory tree 
    > bash scripts/install.sh

    # will generate all the sample utils from tools/sample_manager.py in a C++ header file in src directory. 
    # NB : a compilation will be launched: to be redone (or type make) if a C++ file is modified
    # Relaunch this command each time new input rootfiles are available
    > python scripts/sample_update.py

    # this file need to be update with your path to sourcing ROOT (local machine)
    # NB: Caution, you need to have ROOT-CERN lirbary on your computer. And the rootenv
    # need to be update (if needed) in function of your own ROOT installation path.
    > source rootenv
    # if working on lyoserv or lyoui, use instead (loading devtoolset and an appropriate version of root):
    > source setup.sh

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

The list of option is : [All, mc, data, jec, alt, timed, sme, inclusive]
'sme' create a rootfile with f function in './results/"year"/flattree' directory, usefull for combine inputs creation.

    > ./bin/histograms_creator "observable" "year" "option"

To produce all histograms (takes some times...):

    > ./bin/histograms_creator m_dilep 2017 All

Producing histograms specifically for data/mc comparisons:

    > ./bin/histograms_creator m_dilep 2017 forComp  

Producing nominal MC histograms + weight uncertainties (needed for differential/sme measurement):

    > ./bin/histograms_creator m_dilep 2017 mc

Producing JEC uncertainties histograms (needed for differential/sme measurement):

    > ./bin/histograms_creator m_dilep 2017 jec

Producing alternative uncertainties histograms (needed for differential/sme measurement):

    > ./bin/histograms_creator m_dilep 2017 alt

Alternative uncertainties need an additional layer (requiring the histograms_creator alt command is already done):

    >  python bin/color_reco.py m_dilep 2017 timed

If you want to run the inclusive analysis:

    > ./bin/histograms_creator m_dilep 2017 inclusive

    >  python bin/color_reco.py m_dilep 2017 inclusive

You also need to run on data:

    > ./bin/histograms_creator m_dilep 2017 data


# Comparaison Data/Monte-Carlo

Once the histograms produced with the "forComp" option (same as inclusive setup), data/mc comparaison can be created for a given observable:

    > python ./bin/data_mc_comparaison.py "observable" "year" "title"

Actually, one should now use the data/mc comparisons by splitting in MC processes, using:

    > python ./bin/data_mc_comparaison_splitted.py "observable" "year" "title"


example: 

    > python ./bin/data_mc_comparaison_splitted.py m_dilep 2017 "Dilepton mass (GeV)"

or:

    > python ./bin/data_mc_comparaison_splitted.py n_bjets 2017 "number of b jets"


If willing to produce all data/mc comparisons in one shot:

    > python scripts/launch_hist_creator.py
    > source scripts/launch_datamccomp.bash


# Comparaison signal / backrgounds in the Monte-Carlo

To compare signal and backgrounds (normalized to 1), one can use bin/mc_comparaison.py script.

Example:

    > python bin/mc_comparaison.py n_bjets 2017 "number of b jets"


# Checking Systematics

Checking all systematics at a time (weight uncertainties, ALT and JEC):

    > python scripts/control_systematics.py "year" "timed"

Please check the scripts for command line options, you need to add "timed" or "inclusive"

As an example:

    > python scripts/control_systematics.py 2017 timed

or:

    > python scripts/control_systematics.py 2017 inclusive


By modifying the script, you can choose to plot some systematics and some other not.

It's also possible to check each systematics individually (function called by the control_systematics.py): 

    > python ./bin/systematics_observable.py "observable" "year" "timed or inclusive" "systematic" "title"

example : 

    > python ./bin/systematics_observable.py m_dilep 2017 inclusive syst_elec_id "electron id"

In case you don't sure what systematics are implemented just type :

    > python ./bin/systematics_observable.py

In case you want to check specifically the JEC or ALT systematics (with lower quality plots), you can also use:

    > python bin/check_systematics.py  "observable" "year" "systematic"

For instance:

    > python bin/check_systematics.py m_dilep 2017 mtop

For the special case of mtop, this command line will check into the color_reco.py output the top mass up/down variation where the difference to nominal has been divided by 3 (at the moment, the plots produced by systematics_observable.py or control_systematics.py do not include this).

The result are stored in './results/'year'/other directory'. For now the accessible systematics are : 'jec', 'hdamp', 'CP5'

# Combine input file creation 

For combine work, you will need to create combine inputs and datacards.
All will be stored in the ./combine/year/ directory.

    # combine inputs file
    > python ./bin/combine_"methode".py "observable" "year" 
    # datacards
    > ./bin/card_creator "observable" "year" " "methode"

example :

    > python ./bin/combine_inclusive.py m_dilep 2017
    > ./bin/card_creator m_dilep 2017 Inclusive

    or

    > python ./bin/combine_unrolled.py m_dilep 2017 
    > ./bin/card_creator m_dilep 2017 Unrolled

    or 

    > python ./bin/combine_one_bin.py m_dilep 2016 
    > ./bin/card_creator m_dilep 2016 OneBin

    or 

    > python ./bin/combine_unrolled.py m_dilep 2017
    > ./bin/card_creator m_dilep 2016 SME cLXX

Some usefull scripts can be used in the self-named directory. 

example : 

    This command will create combine inputs then datacards then export in the directory of your servers (adress in ./scripts/export_combine.py)
    > bash scripts/launch_unrolled_stuff.sh m_dilep 2017

# Exporting the combine inputs to combine area

First modify the line 19 of file script/scripts/export_combine.py, to set the path of your combine area, for instance:
 
    # path_out = '/gridgroup/cms/nchanon/CMSSW_10_2_13/src/combine-ttbar/'+directory+'/inputs/'+year+'/'

Export the inputs to inclusive measurement combine area:

    > python scripts/export_combine.py inclusive 2017

Export the inputs to differential measurement combine area:

    > python scripts/export_combine.py one_bin 2017

Export the inputs to SME measurement combine area:

    > python scripts/export_combine.py sme 2017

# Example of full analysis.

Produce data/mc comparisons:

    > python scripts/sample_update.py
    > ./bin/histograms_creator m_dilep 2017 forComp
    > python ./bin/data_mc_comparaison.py m_dilep 2017 "Dilepton mass (GeV)"

Prepare histograms for combine:

    > ./bin/histograms_creator m_dilep 2017 mc
    > ./bin/histograms_creator m_dilep 2017 alt
    > ./bin/histograms_creator m_dilep 2017 jec
    > ./bin/histograms_creator m_dilep 2017 timed
    > python bin/color_reco.py m_dilep 2017

Check systematics:

    > python scripts/control_systematics.py 2017

Prepare combine inputs for inclusive measurement and export it to the combine area:

    > python ./bin/combine_inclusive.py m_dilep 2017
    > ./bin/card_creator m_dilep 2017 Inclusive
    > python scripts/export_combine.py inclusive 2017

Prepare combine inputs for differential measurement and export it to the combine area:

    > python ./bin/combine_one_bin.py m_dilep 2017
    > ./bin/card_creator m_dilep 2017 OneBin
    > python scripts/export_combine.py one_bin 2017

Prepare combine inputs for SME fit and export it to the combine area:

    > python ./bin/combine_unrolled.py m_dilep 2017 
    > ./bin/card_creator m_dilep 2017 SME cLXX
    > python scripts/export_combine.py sme 2017


