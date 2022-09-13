from xml.dom import minidom
import ast
import numpy as np
import matplotlib.pyplot as plt
import os
path = '/Users/katherinepachal/Code/sm-dm-projections-inputs/inputs/directdetection'

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

def parse_lzlims(infile) :
    xvals = []
    yvals = []
    with open(infile) as inf :
        lines = inf.readlines()
        for line in lines[1:] :
            tokens = line.split()
            xvals.append(eval(tokens[0]))
            yvals.append(eval(tokens[1]))
    return xvals, yvals

# All we need is PICO-60
def get_sd_proton() :
    # Read in and format
    file_pico = path+'/SD/PICO60_SDp_1902.04031v2.xml'
    pico_x, pico_y = parse_lims(file_pico)
    file_lz = path+'/SD/Fig8_SpinDependentprotonLimitandSensitivity.txt'
    lz_x, lz_y = parse_lzlims(file_lz)
    return {"PICO-60" : [np.array(pico_x), np.array(pico_y)]} #,
            #"LZ" : [np.array(lz_x), np.array(lz_y)]} # too weak

# Want LUX and XENON1T
def get_sd_neutron() :

    # LUX
    file_lux = path+'/SD/LUX_SDn_1705.03380.xml'
    lux_x, lux_y = parse_lims(file_lux)
    # Lux limit is in pb not cm2, so need to convert.
    lux_y = np.array(lux_y)/1.0e36

    # XENON1T
    file_x1T = path+'/SD/XENON1t_Neutron.xml'
    x1t_x, x1t_y = parse_lims(file_x1T)

    # LZ
    file_lz = path+'/SD/Fig7_SpinDependentneutronLimitandSensitivity.txt'
    lz_x, lz_y = parse_lzlims(file_lz)

    return {'LUX' : [np.array(lux_x), lux_y],
            'XENON1T' : [np.array(x1t_x), np.array(x1t_y)],
            "LZ" : [np.array(lz_x), np.array(lz_y)]}

# Here, we want 3:
# XENON-1T (2018), XENON-1T low mass,
# DarkSide-50
def get_spin_independent() :

    # XENON1T
    file_x1T = path+'/SI/XENON1T_2018.xml'
    x1t_x, x1t_y = parse_lims(file_x1T)
    
    # XENON low mass
    file_x1T_low = path+'/SI/XENON1T_Migdal_LowMass_1907.12771.xml'
    x1t_low_x, x1t_low_y = parse_lims(file_x1T_low)
    
    # DarkSide
    file_darkside = path+'/SI/darkside_1023.xml'
    darkside_x, darkside_y = parse_lims(file_darkside)

    # LZ
    file_lz = path+'/SI/Fig5_SpinIndependentLimitandSensitivity.txt'
    lz_x, lz_y = parse_lzlims(file_lz)    
    
    return {'XENON1T' : [np.array(x1t_x),np.array(x1t_y)],
            'XENON1T MIGD' : [np.array(x1t_low_x),np.array(x1t_low_y)],
            'DarkSide-50' : [np.array(darkside_x),np.array(darkside_y)],
            "LZ" : [np.array(lz_x), np.array(lz_y)]}

# validate
sd_proton = get_sd_proton()
#print(sd_proton)
sd_neutron = get_sd_neutron()
#print(sd_neutron)
si = get_spin_independent()
#print(si)

fig,ax=plt.subplots(1,1)
plt.xlim(1, 2000)
plt.ylim(1e-46, 1e-37)
plt.xscale('log')
plt.yscale('log')
plt.plot(sd_proton['PICO-60'][0],sd_proton['PICO-60'][1])
plt.plot(sd_neutron['LUX'][0],sd_neutron['LUX'][1])
plt.plot(sd_neutron['XENON1T'][0],sd_neutron['XENON1T'][1])
plt.savefig('spin-dependent.pdf',bbox_inches='tight')

plt.clf()
plt.xlim(1, 2000)
plt.ylim(1e-48, 1e-37)
plt.xscale('log')
plt.yscale('log')
plt.plot(si['XENON1T'][0],si['XENON1T'][1])
plt.plot(si['XENON1T MIGD'][0],si['XENON1T MIGD'][1])
plt.plot(si['DarkSide-50'][0],si['DarkSide-50'][1])
plt.savefig('spin-independent.pdf',bbox_inches='tight')
