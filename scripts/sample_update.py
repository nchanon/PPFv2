import sys, os
import argparse

################################################################################
## Initialisation stuff
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')

args = parser.parse_args()
year = args.year

##########

def succeeds_job(percent):
    return 100./(percent)

################################################################################
## 2016 
################################################################################

cross_sec_2016 = [ # alphabetic order
    118.7,                # dibosons WW          
    118.7,                # dibosons WW2         
    47.13,                # dibosons WZ          
    47.13,                # dibosons WZ2         
    16.523,               # dibosons ZZ          
    16.523,               # dibosons ZZ2         

    89.048226,            # signal dilep       
    377.96006,            # signal hadronic    
    366.91429,            # signal semilep     

    10.32,                # singletop STs
    80.95,                # singletop STt antitop
    136.02,               # singletop STt top
    35.85,                # singletop tW  antitop
    35.85,                # singletop tW  top

    0.2043,               # TTX TTW         
    0.2043,               # TTX TTW2        
    0.4062,               # TTX TTW3        
    0.2529,               # TTX TTZ         
    0.2529,               # TTX TTZ2        
    0.2529,               # TTX TTZ3        
    0.5297,               # TTX TTZ4        
    #3.697,                # TTX  TTG        
    #3.697,                # TTX  TTG2        

    61526.7,              # wjets WJets       
    61526.7,              # wjets WJets2      

    6225.4,               # zjets DY        
    #22635.14              # zjets DY 10-50   
    #22635.14              # zjets DY 10-50   
    #22635.14              # zjets DY 10-50   
]

effective_data_event_2016 = [
    succeeds_job(100.),    #MuonEG_B
    succeeds_job(100.),    #MuonEG_C
    succeeds_job(100.),    #MuonEG_D
    succeeds_job(95.8),    #MuonEG_E
    succeeds_job(100.),    #MuonEG_F
    succeeds_job(87.6),    #SingleElectron_B
    succeeds_job(100.),    #SingleElectron_C
    succeeds_job(100.),    #SingleElectron_D
    succeeds_job(83.3),    #SingleElectron_E
    succeeds_job(76.5),    #SingleElectron_F
    succeeds_job(100.),    #SingleMuon_B
    succeeds_job(100.),    #SingleMuon_C
    succeeds_job(100.),    #SingleMuon_D
    succeeds_job(100.),    #SingleMuon_E
    succeeds_job(100.)     #SingleMuon_F
]

################################################################################
## 2017
################################################################################

cross_sec_2017 = [
    118.7,    # dibosons WW
    47.13,    # dibosons WZ
    16.523,   # dibosons ZZ 

    89.05,    # signal dilep
    380.11,   # signal hadronic
    364.31,   # signal semilep

    10.32,    # singletop STs
    10.32,    # singletop STs2
    35.5,     # singletop tW antitop
    35.5,     # singletop tw2 antitop
    35.5,     # singletop tW  top
    35.5,     # singletop tW2 top
    80.95,    # singletop STt antitop
    136.02,   # singletop STt top

    0.2043,   # TTX TTW         
    0.2043,   # TTX TTW2        
    0.4062,   # TTX TTW3        
    0.2529,   # TTX TTZ         
    0.2529,   # TTX TTZ2        
    0.5297,   # TTX TTZ3        

    0.4062,   # wjets WJets       
    0.4062,   # wjets Wjets2      

    22635.1,  # zjets DY 10-50     
    6225.4,   # zjets DY        
    6225.4    # zjets DY2         
]

effective_data_event_2017 = [
    succeeds_job(100.),    #MuonEG_B
    succeeds_job(100.),    #MuonEG_C
    succeeds_job(100.),    #MuonEG_D
    succeeds_job(95.8),    #MuonEG_E
    succeeds_job(100.),    #MuonEG_F
    succeeds_job(87.6),    #SingleElectron_B
    succeeds_job(100.),    #SingleElectron_C
    succeeds_job(100.),    #SingleElectron_D
    succeeds_job(83.3),    #SingleElectron_E
    succeeds_job(76.5),    #SingleElectron_F
    succeeds_job(100.),    #SingleMuon_B
    succeeds_job(100.),    #SingleMuon_C
    succeeds_job(100.),    #SingleMuon_D
    succeeds_job(100.),    #SingleMuon_E
    succeeds_job(100.)     #SingleMuon_F
]


################################################################################
## 2018 
################################################################################

#
#

################################################################################
################################################################################

effective_data_event = {
    '2016' : effective_data_event_2016,
    '2017' : effective_data_event_2017,
    '2018' : 35.9,
}

data_name = {
    '2016' : ['Run2016'],
    '2017' : ['Run2017'],
    '2018' : 0
}

luminosity = {
    '2016' : 35.9,
    '2017' : 41.53,
    '2018' : 0
}

cross_sec = {
    '2016' : cross_sec_2016,
    '2017' : cross_sec_2017,
    '2018' : 0
}


ttbar_list = [ 
    'signal',
    'ttx',
    'singletop',
    'dibosons',
    'wjets',
    'zjets'
]

systematic_list = [
    'syst_elec_reco',
    'syst_elec_id',
    'syst_muon_id',
    'syst_muon_iso',
    'syst_em_trig',
    'syst_pu'
]

trig_2016 = [
    'trg_muon_electron_mu23ele12_fired',
    'trg_muon_electron_mu23ele12DZ_fired',
    'trg_muon_electron_mu8ele23_fired',
    'trg_muon_electron_mu8ele23DZ_fired',
    'trg_muon_mu24_fired',
    'trg_muon_mutk24_fired',
    'trg_electron_ele27_fired'
]

trig_2017 = [
    'trg_muon_electron_mu8ele23DZ_fired',
    'trg_muon_electron_mu23ele12_fired',
    'trg_muon_mu27_fired',
    'trg_electron_ele35_fired'
]

triggers = {
    '2016' : trig_2016,
    '2017' : trig_2017,
#    '2018' : trig_2018
}

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