import pandas as pd
import numpy as np
from generator import get_flux, init_luminosity
from astropy.io import ascii
import astropy.units as u
import matplotlib.pyplot as plt
import os


def read_lum(direc, selected_objs, if_V):
    lum_list = []
    if if_V:
        for name in selected_objs:
            lum_table = pd.read_csv(f"{direc}/{name}_0.csv")
            lum_list.append((lum_table[['t', 'L_V', 'L_u', 'L_g', 'L_r', 'L_i', 'L_z', 'L_y']]))
    else:
        for name in selected_objs:
            lum_table = pd.read_csv(f"{direc}/{name}_0.csv")
            lum_list.append((lum_table[['t', 'L_u', 'L_g', 'L_r', 'L_i', 'L_z', 'L_y']]))
    return lum_list



def init_extin():
    dis_extin = ascii.read("analysis/extinction_1arcsec.csv")
    extinction_baye = dict(zip(dis_extin['name'], dis_extin['extinction_baye']))
    extinction_sfd = dict(zip(dis_extin['name'], dis_extin['extinction_sfd']))
    distances = dict(zip(dis_extin['name'], dis_extin['distance']))
    return extinction_baye, extinction_sfd, distances


extinction_baye, extinction_sfd, distances = init_extin()


def read_extin(objs):

    extin_list = []
    distance_list = []
    for obj_name in objs:
        if np.isnan(extinction_baye[obj_name]):
            extinction = extinction_sfd[obj_name]
        else:
            extinction = extinction_baye[obj_name]
        extin_list.append(extinction)
        distance_list.append(distances[obj_name] * u.parsec)
    return extin_list, distance_list


LSST_A_TO_EBV = {
    'V': 2.742,
    'u': 4.145,
    'g': 3.237,
    'r': 2.273,
    'i': 1.684,
    'z': 1.323,
    'y': 1.088,
}

coef = 2.742


############# compare with PS1 #####
tested_objs = ['OGLE BLG-DN-0036', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0826']
passbands = ['u', 'g', 'r', 'i', 'z', 'y']
modified_pass = [ 'g', 'r', 'i', 'z', 'y', 'V']
# lums = read_lum('analysis_Mdot', tested_objs)
lums = read_lum('analysis_Check_Model', tested_objs, True)
extinctions_list, distances_list = read_extin(tested_objs)
i = 45

for name, l, extinction, distance_parsec in zip(tested_objs, lums, extinctions_list, distances_list):
    t = l['t']
    l.rename(columns={'L_V':'V','L_u': 'u', 'L_g': 'g', 'L_r': 'r', 'L_i': 'i', 'L_z': 'z', 'L_y': 'y'}, inplace=True)
    l_no_t = l.drop(axis=1, labels='t')
    distance_cm = distance_parsec.to(u.cm).value
    flux = get_flux(L=l_no_t, d=distance_cm, i=i)
    mag_noe = -2.5 * np.log10(flux / 3.63e-20)

    plt.figure()
    for passband in modified_pass:
        mag_e = mag_noe[passband] + extinction /coef * LSST_A_TO_EBV[passband]
        plt.plot(t, mag_e, 'o', label=f'{passband}')

    plt.gca().invert_yaxis()
    plt.title(f'{name}')
    plt.legend()
    direc = 'pictures_check_PS1'
    os.makedirs(direc, exist_ok=True)
    plt.savefig(os.path.join(direc, f'{name}_no_normcoef.png'))
    plt.close()


"""
########Compare with OGLE itself#######
tested_objs = ['OGLE BLG-DN-0036', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0826']
modified_pass = ['V']
lums = read_lum('analysis_Check_Model', tested_objs, True)
extinctions_list, distances_list = read_extin(tested_objs)
i = 45

def data_loading(name):
    dff = np.genfromtxt(name, names='t,m,err')
    t = dff['t']
    m = dff['m']
    return t, m

data = ascii.read("analysis/fitting_info.csv")
t_st_dic = dict(zip(data['name'], data['t_start']))
t_ed_dic = dict(zip(data['name'], data['t_end']))
outburst_idx = dict(zip(data['name'], data['outburst_index']))



for name, l, extinction, distance_parsec in zip(tested_objs, lums, extinctions_list, distances_list):
    t = l['t']
    lum_V = l['L_V']
    # l.rename(columns={'L_V': 'V', 'L_u': 'u', 'L_g': 'g', 'L_r': 'r', 'L_i': 'i', 'L_z': 'z', 'L_y': 'y'}, inplace=True)
    # l_no_t = l.drop(axis=1, labels='t')
    distance_cm = distance_parsec.to(u.cm).value
    flux = get_flux(L=lum_V, d=distance_cm, i=i)
    mag_noe = -2.5 * np.log10(flux / 3.63e-20)
    # mag_e = mag_noe + extinction / coef
    mag_e = mag_noe + extinction

    k = name[-4:]
    file_name = "phot/OGLE-BLG-DN-" + k + ".dat"
    obj_name = "OGLE BLG-DN-" + k
    t_all, m = data_loading(file_name)
    t_st = t_st_dic[f'{obj_name}']
    t_ed = t_ed_dic[f'{obj_name}']

    idx = (t_all >= t_st) & (t_all <= t_ed)
    cand_t = t_all[idx]
    cand_m = m[idx]


    plt.figure()
    plt.plot(t, mag_e, 'o', label='derived_magnitude')
    plt.plot(cand_t, cand_m, 'x', label='original_mag')
    # for passband in modified_pass:
    #     # mag_e = mag_noe + extinction / coef * LSST_A_TO_EBV[passband]
    #     mag_e = mag_noe + extinction / coef
    #     plt.plot(t, mag_e[passband], 'o', label=f'{passband}')

    plt.gca().invert_yaxis()
    plt.title(f'{name}')
    plt.legend()
    direc = 'pictures_Check_Model'
    os.makedirs(direc, exist_ok=True)
    plt.savefig(os.path.join(direc, f'{name}_with_Ogle_no_norm.png'))
    plt.close()

"""

"""
### plot PS1 from Kostya#####
from itertools import count
from glob import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
def plot():   
    bands = dict(zip(count(1), 'grizy'))
    for f in glob('PS-1_*.csv'):
        df = pd.read_csv(f)
        df = df[df['psfFlux'] > 0]
        df['psfMag'] = -2.5 * np.log10(df['psfFlux']) + 8.90
        name = f.removesuffix('.csv').removeprefix('PS-1_')
        plt.figure()
        for filter_id in range(1, 6):
            lc = df[df['filterID'] == filter_id]
            plt.plot(lc['obsTime'], lc['psfMag'], 'o', label=bands[filter_id])
        plt.gca().invert_yaxis()
        plt.title(name)
        plt.xlabel('obsTime')
        plt.ylabel('PSF mag')
        plt.legend()
        plt.savefig(f'{name}.png')
        plt.close()
plot()
"""