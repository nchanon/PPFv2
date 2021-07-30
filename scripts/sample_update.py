import sys, os
import argparse

sys.path.append('./')

from tools.sample_manager import *

################################################################################
## Functions stuff
################################################################################

def eventfilter_input(year, nature, sample, filter):
    return './inputs/'+year+'/'+nature+'/'+sample+'/'+filter+'/SkimReport.txt'

def is_same_sample(name1, name2, year):
    foo = False
    if name1.find('TTW3') != -1 or name2.find('TTW3')  != -1 or name1.find('TTZ4') != -1 or name2.find('TTZ4')  != -1:
	return foo
    if year == '2017':
	if name1.find('TTZ3') != -1 or name2.find('TTZ3')  != -1:
	    return foo
    for c1,c2 in zip(name1, name2):
        if c1 != c2:
            try:
                if(int(c1)-int(c2) == 1):
                    return True
            except:
                pass
            if foo:
                return False
        else:
            foo = True
    return foo

def sum_of_weight(year, sample, sampletype='MC'):
    filterfile = eventfilter_input(year, sampletype, sample, 'MCWeighter')
    event = []
    f = open(filterfile, 'r')
    for i in f.read().split():
        try:
            event.append(float(i.strip()))
        except:
            pass
    return float(event[3])

def generate_eventN0(year, sample_list, sampletype='MC'):
    foo = []
    foo.append(sum_of_weight(year,sample_list[0],sampletype))
    for i in range(1,len(sample_list)):
        foo.append(sum_of_weight(year,sample_list[i],sampletype))
        if is_same_sample(sample_list[i], sample_list[i-1], year):
            foo[i-1] += foo[i]
            foo[i] = foo[i-1]
            if is_same_sample(sample_list[i], sample_list[i-2], year):
                foo[i-2] = foo[i]
            if is_same_sample(sample_list[i], sample_list[i-3], year):
                foo[i-3] = foo[i]
                foo[i-2] = foo[i]

    return foo

def rescaling(year, N0_list, sampletype='MC'):
    foo = []
    cross_sec_alt = {
    '2016' : 89.0482256,
    '2017' : 89.0482256
    #    '2018' : trig_2018    
    }
    for n in range(len(N0_list)):
        if sampletype=='MC':
            foo.append(1000.*luminosity[year]*cross_sec[year][n]/N0_list[n]) 
        else:
            foo.append(1000.*luminosity[year]*cross_sec_alt[year]/N0_list[n]) 
    return foo


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
    '2016',
    '2017'
]

for year in years:

    effective_N0 = generate_eventN0(year, sample_MC[year])
    mc_rescale   = rescaling(year, effective_N0)

    jec_effective_N0 = []
    jec_mc_rescale = []
    for l in jec_list:
        jec_effective_N0.append(generate_eventN0(year, sample_MC[year],'JEC'+'/'+l))
        jec_mc_rescale.append(rescaling(year, effective_N0,'JEC'+'/'+l))

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
    core += listing('namelist', 'jecList', jec_list)
    core += listing('namelist', 'systematicList', systematic_list)
    core += listing('namelist', 'systematicAltList', alt_syst_list)
    core += listing('namelist', 'systematicTimeList', systematic_time_list)
    core += listing('namelist', 'systematicRate', systematic_rate[year])
    core += '#endif \n\n'
    core += listing('namelist', 'triggerList_'+year, triggers[year])
    core += listing('namelist', 'sampleList_MC_'+year, sample_MC[year])
    core += listing('namelist', 'sampleList_ALT_'+year, sample_ALT[year])
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