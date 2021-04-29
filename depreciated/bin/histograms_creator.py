import sys
sys.path.append('./')
from tools.generator_manager import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('year', help='year of samples')
parser.add_argument('nature', help='DATA, MC, SYST, time')

args = parser.parse_args()
year = args.year
nature = args.nature

################################################################################
# Generate pure TH1 without fancy style 
################################################################################


#MC
if nature == "MC":
    print 'Start Monte Carlo'

    for i in range(len(sample_list['MC'][year])):
        '''
        generate_TH1('pt_lead',      50, 0, 250, year, 'MC', sample_list['MC'][year][0])
        generate_TH1('pt_sublead',   50, 0, 250, year, 'MC', sample_list['MC'][year][0])
        generate_TH1('j1_pt',        50, 0, 250, year, 'MC', sample_list['MC'][year][0])
        generate_TH1('j2_pt',        50, 0, 250, year, 'MC', sample_list['MC'][year][0])
        generate_TH1('pt_elec',      50, 0, 250, year, 'MC', sample_list['MC'][year][0])
        generate_TH1('pt_muon',      50, 0, 250, year, 'MC', sample_list['MC'][year][0])

        generate_TH1('eta_elec',     50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('eta_muon',     50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('j1_eta',       50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('j2_eta',       50, -3, 3,  year, 'MC', sample_list['MC'][year][0])

        generate_TH1('phi_elec',     50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('phi_muon',     50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('j1_phi',       50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('j2_phi',       50, -3, 3,  year, 'MC', sample_list['MC'][year][0])

        generate_TH1('n_bjets',      5, 0, 5,    year, 'MC', sample_list['MC'][year][0])
        generate_TH1('n_jets_pt30',  7, 0, 7,    year, 'MC', sample_list['MC'][year][0])
        generate_TH1('rho',          50, 0, 80,  year, 'MC', sample_list['MC'][year][0])
        generate_TH1('met',          50, 0, 600, year, 'MC', sample_list['MC'][year][0])
        generate_TH1('metphi',       50, -3, 3,  year, 'MC', sample_list['MC'][year][0])
        '''
        generate_TH1('m_dilep',      50, 0, 800, year, 'MC', sample_list['MC'][year][i])

#        print '----'

    print ''


# DATA
elif nature == "DATA":
    print 'Start Data'

    for i in range(len(sample_list['DATA'][year])):
        '''
         generate_TH1('pt_lead',     50, 0, 250, year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('pt_sublead',  50, 0, 250, year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('pt_elec',     50, 0, 250, year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('pt_muon',     50, 0, 250, year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('j1_pt',       50, 0, 250, year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('j2_pt',       50, 0, 250, year, 'DATA', sample_list['DATA'][year][i])
         
         generate_TH1('eta_muon',    50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('eta_elec',    50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('j1_eta',      50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('j2_eta',      50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         
         generate_TH1('phi_elec',    50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('phi_muon',    50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('j1_phi',      50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('j2_phi',      50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         
         generate_TH1('n_bjets',     5, 0, 5,    year, 'DATA', sample_list['DATA'][year][i])
         #generate_TH1('n_jets_pt30', 7, 0, 7,    year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('rho',         50, 0, 80,  year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('met',         50, 0, 600, year, 'DATA', sample_list['DATA'][year][i])
         generate_TH1('metphi',      50, -3, 3,  year, 'DATA', sample_list['DATA'][year][i])
         '''
        generate_TH1('n_bjets', 5, 0, 5, year, 'DATA', sample_list['DATA'][year][i])
        generate_TH1('m_dilep', 25, 0, 300, year, 'DATA', sample_list['DATA'][year][i])

        print '----'
         
    print ''

# DATA
elif nature == "time":

    os.system("rm "+results_path(year,'number_of_events','data.txt'))
    nbin = int(raw_input("number of bin : "))
    print 'Start Data'

    for i in range(len(sample_list['DATA'][year])):
        #generate_DATA_unrolled('m_dilep', nbin, year, sample_list['DATA'][year][i])
        generate_DATA_timed('m_dilep', nbin, year, sample_list['DATA'][year][i])
        print '----'
    print ''


#SYST
elif nature == "SYST":
    print 'Start Systematics'

    for i in range(len(sample_list['MC'][year])):
        for j in range(len(systematic_list)):
            generate_TH1_systematic('m_dilep', 50, 0, 800, year, systematic_list[j], 'Up', sample_list['MC'][year][i])
            generate_TH1_systematic('m_dilep', 50, 0, 800, year, systematic_list[j], 'Down',sample_list['MC'][year][i])
        print '----'

    print ''

else:
    print "Wrong nature of sample (DATA, MC, SYST)"
