from directory_manager import *
from sample_manager import *
from tools.style_manager import *

import numpy as np

from ROOT import TFile, TH1F

################################################################################
# General
################################################################################

# for MC
def generate_weight(tree):
    foo = 1.
    foo *= tree.GetLeaf('weight_pu').GetValue()
    foo *= tree.GetLeaf('weight_generator').GetValue()
    foo *= tree.GetLeaf('weight_top').GetValue()
    #lepton
    foo *= tree.GetLeaf('weight_sfe_id').GetValue()
    foo *= tree.GetLeaf('weight_sfe_reco').GetValue()
    foo *= tree.GetLeaf('weight_sfm_id').GetValue()
    foo *= tree.GetLeaf('weight_sfm_iso').GetValue()
    #hadron
    for i in range(int(tree.GetLeaf('n_bjets').GetValue())):
        foo *= tree.GetLeaf('weight_sfb').GetValue()
    #trigger
    foo *= tree.GetLeaf('weight_sf_em_trig').GetValue()
    return foo

# for DATA
def trigger_passed(tree, trigger_list):
    trig = 0
    for l in trigger_list:
        if(tree.GetLeaf(l).GetValue() == 1):
            trig += 1
    if(trig>1):
        return True
    return False


def generate_trigger(year, tree, hist, variable, file_name):
    if(trigger_passed(tree, triggers[year])):
        if (file_name.find('MuonEG')!=-1):
            hist.Fill(tree.GetLeaf(variable).GetValue())
            return
            if (file_name.find('SingleMuon')!=-1):
                hist.Fill(tree.GetLeaf(variable).GetValue())
                return
                if (file_name.find('SingleElectron')!=-1):
                    hist.Fill(tree.GetLeaf(variable).GetValue())
                    return
    return

def generate_trigger_timed(year, tree, variable, file_name, nbin):
    list = [ [ 0 for i in range(nbin)],[ 0 for i in range(nbin)] ]
    if(trigger_passed(tree, triggers[year])):
        if (file_name.find('MuonEG')!=-1):
            bin_foo = int(np.floor(sideral_time(tree.GetLeaf('unix_time').GetValue())%nbin))
            list[0][bin_foo] += tree.GetLeaf(variable).GetValue()
            list[1][bin_foo] +=1
            return list
            if (file_name.find('SingleMuon')!=-1):
                bin_foo = int(np.floor(sideral_time(tree.GetLeaf('unix_time').GetValue())%nbin))
                list[0][bin_foo] += tree.GetLeaf(variable).GetValue()
                list[1][bin_foo] +=1                
                return list
                if (file_name.find('SingleElectron')!=-1):
                    bin_foo = int(np.floor(sideral_time(tree.GetLeaf('unix_time').GetValue())%nbin))
                    list[0][bin_foo] += tree.GetLeaf(variable).GetValue()
                    list[1][bin_foo] +=1
                    return list
    return list

def read_DATA_variable(year, tree, variable, file_name):
    number = 0
    if(trigger_passed(tree, triggers[year])):
        if (file_name.find('MuonEG')!=-1):
            number = tree.GetLeaf(variable).GetValue()
            return number
            if (file_name.find('SingleMuon')!=-1):
                number = tree.GetLeaf(variable).GetValue()
                return number
                if (file_name.find('SingleElectron')!=-1):
                    number = tree.GetLeaf(variable).GetValue()
                    return number
    return number


# for syst
def generate_syst(tree, systematic_name, up_down):
    if(systematic_name == 'syst_pu'):
        if(up_down == 'Up'):
            return tree.GetLeaf('weight_pu_up').GetValue()
        elif(up_down == 'Down'):
            return tree.GetLeaf('weight_pu_down').GetValue()
    else:
        foo = tree.GetLeaf(systematic_name).GetValue()
        if(up_down == 'Up'):
            return 1+foo
        elif(up_down == 'Down'):
            return 1-foo
        else:
            print 'up_down error'

def write(hist, variable, year, where, nature, sample):
    newfile = TFile(file_inout('results', year, where, nature , sample+'_'+variable+'.root'), 'RECREATE')
    hist.Write()
    newfile.Close()


################################################################################
# For groups
################################################################################

def grouping(n_bin, bin_min, bin_max, year, hist, grp):
    foo = [[] for x in xrange(len(analysis_list[year]))]
    for h in hist:
        for i in range(len(h)):
            h[i].append(TH1F(grp+var, var, n_bin, bin_min, bin_max))
            for index, grp in enumerate(analysis_list):
                if(h[i].GetName().find(grp)!=-1):
                    foo[i][index] += h[i][index]
        return False


################################################################################
################################################################################
################################################################################

def generate_MC_groups(year, variables_list, analysis_list):

   #hist = [[] for x in xrange(len(sample_list_DATA[year]))]
    hist = []
    for ind, grp in enumerate(analysis_list):
        hist.append([])
        for var in variables_list:
            binning = observable_values(var)[0]
            hist[ind].append(TH1F(grp+'_'+var, var, binning[0], binning[1], binning[2]))

    for sample_index, sample in enumerate(sample_list['MC'][year]):
        rfile  = TFile(heppy_tree(year, 'MC', sample))
        tree   = rfile.Get('events')
        h      = []
        for var in variables_list:
            binning = observable_values(var)[0]
            h.append(TH1F(sample+'_'+var, var, binning[0], binning[1], binning[2]))
        for i in range(tree.GetEntriesFast()):
            tree.GetEntry(i)
            sf_weight   = generate_weight(tree)
            if(trigger_passed(tree, triggers[year])):
                for ind, var in enumerate(variables_list):
                    h[ind].Fill(tree.GetLeaf(var).GetValue(),sf_weight)
            if(i%100000 == 0):
                print '100 000 events passed'
        for i in range(len(variables_list)):
            h[i].Scale(mc_rescale(year, sample))
            for grp_index, grp in enumerate(analysis_list):
                if(h[i].GetName().find(grp) != -1):
                    hist[grp_index][i].Add(h[i])
        print 'next sample'

    for grp_i, grp in enumerate(analysis_list):
        for var_i, var in enumerate(variables_list):
            write(hist[grp_i][var_i], var, year, 'groups', 'MC', grp)

def generate_SYST_groups(year, variables_list, systematic, analysis_list):

   #hist = [[] for x in xrange(len(sample_list_DATA[year]))]
    hist_syst_up = []
    hist_syst_down = []
    for ind, grp in enumerate(analysis_list):
        hist_syst_up.append([])
        hist_syst_down.append([])
        for var in variables_list:
            binning = observable_values(var)[0]
            hist_syst_up[ind].append(TH1F(grp+'_'+var+'_'+systematic+'Up', var, binning[0], binning[1], binning[2]))
            hist_syst_down[ind].append(TH1F(grp+'_'+var+'_'+systematic+'Down', var, binning[0], binning[1], binning[2]))

    for sample_index, sample in enumerate(sample_list['MC'][year]):
        rfile  = TFile(heppy_tree(year, 'MC', sample))
        tree   = rfile.Get('events')
        h_up   = []
        h_down = []
        for var in variables_list:
            binning = observable_values(var)[0]
            h_up.append(TH1F(sample+'_'+var+'_'+systematic+'Up', var, binning[0], binning[1], binning[2]))
            h_down.append(TH1F(sample+'_'+var+'_'+systematic+'Down', var, binning[0], binning[1], binning[2]))
        for i in range(tree.GetEntriesFast()):
            tree.GetEntry(i)
            sf_weight = generate_weight(tree)
            syst_up   = generate_syst(tree, systematic, 'Up')
            syst_down = generate_syst(tree, systematic, 'Down')
            if(trigger_passed(tree, triggers[year])):
                for ind, var in enumerate(variables_list):
                    h_up[ind].Fill(tree.GetLeaf(var).GetValue(),sf_weight*syst_up)
                    h_down[ind].Fill(tree.GetLeaf(var).GetValue(),sf_weight*syst_down)
            if(i%100000 == 0):
                print '100 000 events passed'
        for i in range(len(variables_list)):
            h_up[i].Scale(mc_rescale(year, sample))
            h_down[i].Scale(mc_rescale(year, sample))

            for grp_index, grp in enumerate(analysis_list):
                if(h_up[i].GetName().find(grp) != -1):
                    hist_syst_up[grp_index][i].Add(h_up[i])
                    hist_syst_down[grp_index][i].Add(h_down[i])
        print 'next sample'

    for grp_i, grp in enumerate(analysis_list):
        for var_i, var in enumerate(variables_list):
            write(hist_syst_up[grp_i][var_i], var+'_'+systematic+'Up', year, 'groups', 'SYST', grp)
            write(hist_syst_down[grp_i][var_i], var+'_'+systematic+'Down', year, 'groups', 'SYST', grp)


def generate_DATA_timed(variable, n_bin, year, sample):
    rfile = TFile(heppy_tree(year, 'DATA', sample))
    hist  = TH1F(variable+'_time'+str(n_bin), variable+'_time'+str(n_bin), n_bin, 0, n_bin)
    tree  = rfile.Get('events')

    content = [ [ 0 for i in range(n_bin)],[ 0 for i in range(n_bin)] ] 

    #######################
    for i in range(tree.GetEntriesFast()):
        tree.GetEntry(i)

        bin_foo = int(np.floor(sideral_time(tree.GetLeaf('unix_time').GetValue())%n_bin))        
        if(trigger_passed(tree, triggers[year])):
            if (sample.find('MuonEG')!=-1):
                content[0][bin_foo] += tree.GetLeaf(variable).GetValue()
                content[1][bin_foo] +=1
                continue
                if (sample.find('SingleMuon')!=-1):
                    content[0][bin_foo] += tree.GetLeaf(variable).GetValue()
                    content[1][bin_foo] +=1                
                    continue
                    if (sample.find('SingleElectron')!=-1):
                        content[0][bin_foo] += tree.GetLeaf(variable).GetValue()
                        content[1][bin_foo] +=1
                        continue
    print sample, content
    total = 0
    for i in range(n_bin):
        if(content[0][i] != 0 and content[1][i] != 0 ):
            hist.SetBinContent(i+1, content[0][i]/content[1][i])
        else:
            hist.SetBinContent(i+1, 0)
        total += content[1][i]

    hist.Scale(effective_data_event[year][index(year,sample,'DATA')])
    ########################
    file = open(results_path(year,'number_of_events','data.txt'),"a") 
    file.write(sample+'   '+str(total)+'\n')
    file.close()

    write(hist, variable+'_time'+str(n_bin), year, 'TH1', 'DATA', sample)
    rfile.Close()



def generate_DATA_unrolled(variable, n_time, year, sample):
    rfile = TFile(heppy_tree(year, 'DATA', sample))
    binning = observable_values(variable)[0]
    hist = []
    for i in range(n_time):
        hist.append(TH1F(variable+'_time'+str(n_time)+'_'+str(i), variable+'_time'+str(n_time)+'_'+str(i), binning[0], binning[1], binning[2]))
    tree  = rfile.Get('events')
    ########################
    for i in range(tree.GetEntriesFast()):
        tree.GetEntry(i)
        ibin = int(np.floor(sideral_time(tree.GetLeaf('unix_time').GetValue())%n_time))
        generate_trigger(year, tree, hist[ibin], variable, sample)
    ########################
    for i in range(n_time):
        hist[i].Scale(effective_data_event[year][index(year,sample,'DATA')])
        write(hist[i], variable+'_time'+str(n_time)+'_'+str(i), year, 'TH1', 'DATA/time', sample)
    rfile.Close()




def generate_TH1(variable, n_bin, bin_min, bin_max, 
                    year, nature, sample):
    rfile = TFile(heppy_tree(year, nature, sample))
    hist  = TH1F(variable, variable, n_bin, bin_min, bin_max)
    tree  = rfile.Get('events')

    ########################
    if(nature == 'MC'):
        
        for i in range(tree.GetEntriesFast()):
            tree.GetEntry(i)
            sf_weight   = generate_weight(tree)
            if(trigger_passed(tree, triggers[year])):
                hist.Fill(tree.GetLeaf(variable).GetValue(),sf_weight)
            if(i%100000 == 0):
                print '100 000 events passed'
        hist.Scale(effective_mc_event[year][index(year,sample,'MC')])
        hist.Scale(mc_rescale(year, sample))
    ########################
    elif(nature == 'DATA'):
        for i in range(tree.GetEntriesFast()):
            tree.GetEntry(i)
            generate_trigger(year, tree, hist, variable, sample)
        hist.Scale(effective_data_event[year][index(year,sample,'DATA')])
    ########################
    write(hist, variable, year, 'TH1', nature, sample)
    rfile.Close()


def generate_TH1_systematic(variable, n_bin, bin_min, bin_max, year, systematic_name, up_down, sample):
    rfile = TFile(heppy_tree(year, 'MC', sample))
    hist  = TH1F(variable, variable, n_bin, bin_min, bin_max)
    tree  = rfile.Get('events')

    for i in range(tree.GetEntriesFast()):
        tree.GetEntry(i)
        sf_weight   = generate_weight(tree)
        syst_weight = generate_syst(tree, systematic_name, up_down)
        hist.Fill(tree.GetLeaf(variable).GetValue(),sf_weight*syst_weight)
        if(i%100000 == 0):
            print '100 000 events passed'
    hist.Scale(mc_rescale(year, sample))

    newfile = TFile(file_inout('results', year, 'TH1', 'SYST' , sample+'_'+variable+'_'+systematic_name+up_down+'.root'), 'RECREATE')
    hist.Write()
    newfile.Close()
    rfile.Close()