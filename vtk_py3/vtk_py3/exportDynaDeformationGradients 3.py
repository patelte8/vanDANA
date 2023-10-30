########################################################################

import os

########################################################################
#open d3plot "''' + d3plot_file_basename + '''.d3plot"

def exportDynaDeformationGradients(d3plot_file_basename, stateno):
    lspp_file_name = d3plot_file_basename+'.lspp'
    lspp_file = open(lspp_file_name, 'w')
    lspp_file.write('''\
open d3plot "d3plot"
selectpart shell off
fringe 91
output "''' + d3plot_file_basename + '''.history#11"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 92
output "''' + d3plot_file_basename + '''.history#12"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 93
output "''' + d3plot_file_basename + '''.history#13"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 94
output "''' + d3plot_file_basename + '''.history#14"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 95
output "''' + d3plot_file_basename + '''.history#15"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 96
output "''' + d3plot_file_basename + '''.history#16"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 97
output "''' + d3plot_file_basename + '''.history#17"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 98
output "''' + d3plot_file_basename + '''.history#18"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 99
output "''' + d3plot_file_basename + '''.history#19"''' + str(stateno) + ''' 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
exit
''')
    lspp_file.close()
#    os.system('Xvfb :2 -screen 0 1074x800x24 &')
    os.system('/opt/lsprepost4.1_centos6/lspp41 -nographics c='+lspp_file_name)
