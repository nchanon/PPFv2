import sys, os
import argparse

sys.path.append('./')

from tools.sample_manager import *

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year


################################################################################
## Functions stuff
################################################################################

def inputs_name(path, year, nature):
    foo = os.listdir('./'+path+'/'+year+'/'+nature+'/')
    foo.sort()
    if(not foo):
        print 'Error : dataset list is empty'
    return foo

def eventfilter_input(year, nature, sample, filter):
    return './inputs/'+year+'/'+nature+'/'+sample+'/0000/'+sample+'/'+filter+'/SkimReport.txt'

def is_same_sample(name1, name2):
    foo = False;
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

def sum_of_weight(year, sample):
    filterfile = eventfilter_input(year, 'MC', sample, 'MCWeighter')
    event = []
    f = open(filterfile, 'r')
    for i in f.read().split():
        try:
            event.append(float(i.strip()))
        except:
            pass
    return float(event[3])

def generate_eventN0(year, sample_list):
    foo = []
    foo.append(sum_of_weight(year,sample_list[0]))
    for i in range(1,len(cross_sec[year])):
        foo.append(sum_of_weight(year,sample_list[i]))
        if is_same_sample(sample_list[i], sample_list[i-1]):
            foo[i-1] += foo[i]
            foo[i] = foo[i-1]
            if is_same_sample(sample_list[i], sample_list[i-2]):
                foo[i-2] = foo[i]
            if is_same_sample(sample_list[i], sample_list[i-3]):
                foo[i-3] = foo[i]
                foo[i-2] = foo[i]

    return foo

def rescaling(year, N0_list):
    foo = []
    for n in range(len(N0_list)):
        foo.append(1000.*luminosity[year]*cross_sec[year][n]/N0_list[n]) 
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


################################################################################
## Code body
################################################################################

sample_MC    = inputs_name('inputs', year, 'MC')
sample_DATA  = inputs_name('inputs', year, 'DATA')

effective_N0 = generate_eventN0(year, sample_MC)
mc_rescale   = rescaling(year, effective_N0)

core = "// "+year+" samples : \n\n"
core += '#include <vector> \n'
core += '#include <string> \n'
core += ' \n'
core += 'using  namelist = std::vector<std::string>; \n'
core += ' \n'
core += '#ifndef COMMON_LIST \n'
core += '#define COMMON_LIST \n'
core += listing('namelist', 'triggerList', triggers[year])
core += listing('namelist', 'ttbarList', ttbar_list)
core += listing('namelist', 'systematicList', systematic_list)
core += '#endif \n\n'

core += listing('namelist', 'systematicRate', systematic_rate[year])
core += listing('namelist', 'sampleList_MC_'+year, sample_MC)
core += listing('std::vector<double>', 'mc_rescale_'+year, mc_rescale)
core += listing('namelist', 'sampleList_DATA_'+year, sample_DATA)
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