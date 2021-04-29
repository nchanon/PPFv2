import sys
sys.path.append('./')
from tools.generator_manager import *

for index, i in enumerate(events_N0['2016']):
    print sample_list['MC']['2016'][index], i
print ''
for index, i in enumerate(events_N0['2017']):
    print sample_list['MC']['2017'][index], i