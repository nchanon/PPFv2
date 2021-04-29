from directory_manager import *

from manager2016 import *
from manager2017 import *


################################################################################
## Analysis values
################################################################################

# constants
omega_UT1  = 7.2921e-5
omega_GMST = 7.2722e-5
t0         = 1483228800.
phase      = 3.2830


ttbar_list = [ 
    'signal',
    'ttx',
    'singletop',
    'dibosons',
    'wjets',
    'zjets'
]

variables = [
    'm_dilep',
    'n_bjets'
    #'unix_time'
]

systematic_list = [
    'syst_elec_reco',
    'syst_elec_id',
    'syst_muon_id',
    'syst_muon_iso',
    'syst_em_trig',
    'syst_pu'
]


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
################################################################################
################################################################################


sample_list_groups_2016 = inputs_name('results', '2016', 'groups/MC')
sample_list_groups_2017 = inputs_name('results', '2017', 'groups/MC')
#
sample_list_groups_SYST_2016 = inputs_name('results', '2016', 'groups/SYST')
sample_list_groups_SYST_2017 = inputs_name('results', '2017', 'groups/SYST')
#
sample_list_MC_2016 = inputs_name('inputs', '2016', 'MC')
sample_list_MC_2017 = inputs_name('inputs', '2017', 'MC')
#sample_list_MC_2018 = inputs_name('inputs', '2018', 'MC')
sample_list_DATA_2016 = inputs_name('inputs', '2016', 'DATA')
sample_list_DATA_2017 = inputs_name('inputs', '2017', 'DATA')
#sample_list_DATA_2018 = inputs_name('inputs', '2018', 'DATA')

sample_list_groups = {
    '2016' : sample_list_groups_2016,
    '2017' : sample_list_groups_2017,
    #'2018' :
}

sample_list_groups_SYST = {
    '2016' : sample_list_groups_SYST_2016,
    '2017' : sample_list_groups_SYST_2017,
    #'2018' :
}

sample_list_MC = {
    '2016' : sample_list_MC_2016,
    '2017' : sample_list_MC_2017,
#    '2018' : sample_list_MC_2018
}

sample_list_DATA = {
    '2016' : sample_list_DATA_2016,
    '2017' : sample_list_DATA_2017,
#    '2018' : sample_list_DATA_2018
}

triggers = {
    '2016' : trig_2016,
    '2017' : trig_2017,
#    '2018' : trig_2018
}

elecmu_trig =  {
    '2016' : elecmu_trig_2016,
    '2017' : elecmu_trig_2017,
#    '2018' : elecmu_trig_2018
}

mu_trig = {
    '2016' : mu_trig_2016,
    '2017' : mu_trig_2017,
#    '2018' : mu_trig_2018
}

elec_trig = {
    '2016' : elec_trig_2016,
    '2017' : elec_trig_2017,
#    '2018' : elec_trig_2018
}

# To call sample_list do sample_list[nature][year]
sample_list = {
    'MC'   : sample_list_MC,
    'DATA' : sample_list_DATA
}


################################################################################
# Utils
################################################################################

## time ones

def sideral_time(time_user):
    return (omega_UT1 * (time_user - t0))/omega_GMST
##

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

def index(year, sample, dataset='MC'):
    foo = 0
    for i in sample_list[dataset][year]:
        if(i == sample):
            return foo
        foo += 1 

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

def generate_eventN0(year):
    foo = []
    foo.append(sum_of_weight(year,sample_list['MC'][year][0]))
    for i in range(1,len(cross_sec[year])):
        foo.append(sum_of_weight(year,sample_list['MC'][year][i]))
        if is_same_sample(sample_list['MC'][year][i], sample_list['MC'][year][i-1]):
            foo[i-1] += foo[i]
            foo[i] = foo[i-1]
            if is_same_sample(sample_list['MC'][year][i], sample_list['MC'][year][i-2]):
                foo[i-2] = foo[i]
            if is_same_sample(sample_list['MC'][year][i], sample_list['MC'][year][i-3]):
                foo[i-3] = foo[i]
                foo[i-2] = foo[i]

    return foo

def percent(x, total):
    return round(100.* x/total, 2)

################################################################################
# Variables
################################################################################

#events_N0_2016 = []
#for i in range(len(cross_sec_2016)):
#    events_N0_2016.append(sum_of_weight('2016',sample_list['MC']['2016']#[i]))
#
#events_N0_2017 = []
#for i in range(len(cross_sec_2017)):
#    events_N0_2017.append(sum_of_weight('2017',sample_list['MC']['2017']#[i]))

effective_data_event = {
    '2016' : effective_data_event_2016,
    '2017' : effective_data_event_2017,
    '2018' : 35.9,
}


events_N0 = {
    '2016' : generate_eventN0('2016'),
    '2017' : generate_eventN0('2017'),
    '2018' : 0
}


def mc_rescale(year, sample):
    foo = 1000.*luminosity[year]*cross_sec[year][index(year, sample)]
    return foo/(events_N0[year][index(year, sample)])
