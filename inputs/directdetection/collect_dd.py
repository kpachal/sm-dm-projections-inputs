from xml.dom import minidom
import ast
import numpy as np
import matplotlib.pyplot as plt

def parse_lims(infile) :
    openfile = minidom.parse(infile)
    limits = openfile.getElementsByTagName('limit')[0]
    vals = limits.getElementsByTagName('data-values')[0].firstChild.data
    vals_list = vals.strip("{[").strip("]}")
    vals_groups = vals_list.split(";")
    xvals = []
    yvals = []
    for vals_group in vals_groups :
        tokens = vals_group.split()
        xvals.append(eval(tokens[0]))
        yvals.append(eval(tokens[1]))
    return xvals, yvals

# All we need is PICO-60
def get_sd_proton() :
    # Read in and format
    file_pico = 'SD/PICO60_SDp_1902.04031v2.xml'
    pico_x, pico_y = parse_lims(file_pico)
    return {"PICO-60" : [np.array(pico_x), np.array(pico_y)]}

# Want LUX and XENON-1T
def get_sd_neutron() :

    # LUX
    file_lux = 'SD/LUX_SDn_1705.03380.xml'
    lux_x, lux_y = parse_lims(file_lux)
    # Lux limit is in pb not cm2, so need to convert.
    lux_y = np.array(lux_y)/1.0e36

    # XENON-1T
    file_x1T = 'SD/XENON1t_Neutron.xml'
    x1t_x, x1t_y = parse_lims(file_x1T)
    
    return {'LUX' : [np.array(lux_x), lux_y],
            'XENON-1T' : [np.array(x1t_x), np.array(x1t_y)]}

# Here:
def get_spin_independent() :
    pass

# validate
sd_proton = get_sd_proton()
print(sd_proton)
sd_neutron = get_sd_neutron()
print(sd_neutron)

fig,ax=plt.subplots(1,1)
plt.xlim(1, 2000)
plt.ylim(1e-46, 1e-37)
plt.xscale('log')
plt.yscale('log')
plt.plot(sd_proton['PICO-60'][0],sd_proton['PICO-60'][1])
plt.plot(sd_neutron['LUX'][0],sd_neutron['LUX'][1])
plt.plot(sd_neutron['XENON-1T'][0],sd_neutron['XENON-1T'][1])
plt.savefig('spin-dependent.pdf',bbox_inches='tight')

