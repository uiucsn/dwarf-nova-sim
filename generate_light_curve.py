import pandas as pd
import numpy as np
from generator import get_flux, init_luminosity
from astropy.io import ascii
import astropy.units as u
import matplotlib.pyplot as plt


def read_lum(direc, selected_objs):
    lum_list = []
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


tested_objs = ['OGLE BLG-DN-0036', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0826']
passbands = ['u', 'g', 'r', 'i', 'z', 'y']

lums = read_lum('analysis_Mdot', tested_objs)
extinctions_list, distances_list = read_extin(tested_objs)
i = 45

for name, l, extinction, distance_parsec in zip(tested_objs, lums, extinctions_list, distances_list):
    t = l['t']
    l.rename(columns={'L_u': 'u', 'L_g': 'g', 'L_r': 'r', 'L_i': 'i', 'L_z': 'z', 'L_y': 'y'}, inplace=True)
    l_no_t = l.drop(axis=1, labels='t')
    distance_cm = distance_parsec.to(u.cm).value
    flux = get_flux(L=l_no_t, d=distance_cm, i=i)
    mag_noe = -2.5 * np.log10(flux / 3.63e-20)
    #Extinction?

    for passband in passbands:
        plt.plot(t, mag_noe[passband], 'o', label=f'{name}_{passband}')
        plt.title(f'{name}_{passband}')
        plt.gca().invert_yaxis()
        plt.show()





"""
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