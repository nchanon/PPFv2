import sys, os
import argparse

sys.path.append('./')

from tools.sample_manager import *

################################################################################
## Functions stuff
################################################################################

def listing(cpp_type, name, content):
    if(cpp_type.find('namelist') != -1):
        binding = '"'
    else:
        binding = '' 
    foo = cpp_type+' '+name+' {\n'
    for i in range(len(content)):
        foo += '    '+binding+str(content[i])+binding
        if(i != len(content)-1):
            foo +=','
        foo += '\n'
    foo += '};\n\n'
    return foo


def multilisting(cpp_type, name, content):
    if(cpp_type.find('namelist') != -1):
        binding = '"'
    else:
        binding = '' 
    foo = 'std::vector<'+cpp_type+' > '+name+' {\n'
    for l in range(len(content)):
        foo += '    {\n'
        for i in range(len(content[l])):
            foo += '        '+binding+str(content[l][i])+binding
            if(i != len(content[l])-1):
                foo +=','
            foo += '\n'
        foo += '    }'
        if(l != len(content)-1):
            foo += ','
        foo += '\n'        
    foo += '};\n\n'
    return foo


################################################################################
## Code body
################################################################################

years = [
    #'2016',
    '2017'
]

for year in years:

    numberofevents_N0 = []
    effective_N0 = []
    mc_rescale = []
    numberofevents_N0 = generate_numberofevents(year, sample_MC[year]) 
    effective_N0 = generate_eventN0(year, sample_MC[year])
    mc_rescale   = rescaling(year, effective_N0)

    print effective_N0
    print mc_rescale

    #jec_effective_N0 = []
    #jec_mc_rescale = []
    #for l in jec_list:
    #    jec_effective_N0.append(generate_eventN0(year, sample_MC[year],'JEC'+'/'+l))
    #    jec_mc_rescale.append(rescaling(year, effective_N0,'JEC'+'/'+l))

    jec_mc_rescale = []
    for l in jec_list:
         jec_mc_rescale.append(mc_rescale)

    alt_effective_N0 = generate_eventN0(year, sample_ALT[year],'ALT')
    alt_mc_rescale   = rescaling(year, alt_effective_N0,'ALT')

    core = "// "+year+" samples : \n\n"
    core += '#include <vector> \n'
    core += '#include <string> \n'
    core += ' \n'
    core += 'using  namelist = std::vector<std::string>; \n'
    core += ' \n'
    core += '#ifndef COMMON_LIST \n'
    core += '#define COMMON_LIST \n'
    core += listing('namelist', 'ttbarList', ttbar_list)
    core += listing('namelist', 'systematicList', systematic_list)
    core += listing('namelist', 'systematicAltList', alt_syst_list)
    core += listing('namelist', 'systematicTimeList', systematic_time_list)
    core += listing('namelist', 'systematicRate', systematic_rate[year])
    core += '#endif \n\n'
    core += listing('namelist', 'triggerList_'+year, triggers[year])
    core += listing('namelist', 'jecList_'+year, jec_list[year])
    core += listing('namelist', 'sampleList_MC_'+year, sample_MC[year])
    core += listing('namelist', 'sampleList_ALT_'+year, sample_ALT[year])
    core += listing('std::vector<double>', 'number_of_events_'+year, numberofevents_N0)
    core += listing('std::vector<double>', 'mc_rescale_'+year, mc_rescale)
    core += listing('std::vector<double>', 'alt_mc_rescale_'+year, alt_mc_rescale)
    core += multilisting('std::vector<double>', 'jec_mc_rescale_'+year, jec_mc_rescale)
    core += listing('namelist', 'sampleList_DATA_'+year, sample_DATA[year])
    core += listing('std::vector<double>', 'succedJobs_'+year, effective_data_event[year])
    core += listing('namelist', 'data_'+year, data_name[year])

    
    file = open("./src/sample_"+year+".hpp",'w') 
    file.write(core) 
    file.close() 

################################################################################
## Rebuild project
################################################################################

os.system('make clean')
os.system('make')
