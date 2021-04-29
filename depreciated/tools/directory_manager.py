import sys, os

def inputs_name(path, year, nature):
    foo = os.listdir('./'+path+'/'+year+'/'+nature+'/')
    foo.sort()
    if(not foo):
        print 'Error : dataset list is empty'
    return foo

def file_inout(path, year, nature, sample, file):
    return './'+path+'/'+year+'/'+nature+'/'+sample+'/'+file

def results_path(year, directory, file):
    return './results/'+year+'/'+directory+'/'+file

def heppy_tree(year, nature, sample):
    return file_inout('inputs', year, nature, sample, 'tree.root')

def eventfilter_input(year, nature, sample, filter):
    return './inputs/'+year+'/'+nature+'/'+sample+'/0000/'+sample+'/'+filter+'/SkimReport.txt'

def efficiency_input(year, nature, sample, filter):
    return './inputs/'+year+'/'+nature+'/'+sample+'/0000/'+sample+'/'+filter+'/efficiency.txt'

