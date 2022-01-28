import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
from ch_vars.spatial_distr import MilkyWayDensityJuric2008 as MWDensity
from ch_vars.extinction import get_sfd_thin_disk_ebv as get_extinction
import os

passbands = ['u', 'g', 'r', 'i', 'z', 'y']


# selected_objs = ['OGLE BLG-DN-0001_3', 'OGLE BLG-DN-0001_4', 'OGLE BLG-DN-0002_0', 'OGLE BLG-DN-0002_1',
#                  'OGLE BLG-DN-0002_2', 'OGLE BLG-DN-0036_0', 'OGLE BLG-DN-0087_0', 'OGLE BLG-DN-0168_0',
#                  'OGLE BLG-DN-0174_0', 'OGLE BLG-DN-0233_1', 'OGLE BLG-DN-0286_0', 'OGLE BLG-DN-0305_2',
#                  'OGLE BLG-DN-0373_0', 'OGLE BLG-DN-0376_1', 'OGLE BLG-DN-0421_2', 'OGLE BLG-DN-0444_0',
#                  'OGLE BLG-DN-0531_0', 'OGLE BLG-DN-0588_0', 'OGLE BLG-DN-0595_0', 'OGLE BLG-DN-0690_2',
#                  'OGLE BLG-DN-0783_0', 'OGLE BLG-DN-0826_0', 'OGLE BLG-DN-0899_0']

selected_objs = ['OGLE BLG-DN-0001', 'OGLE BLG-DN-0001', 'OGLE BLG-DN-0002', 'OGLE BLG-DN-0002',
                 'OGLE BLG-DN-0002', 'OGLE BLG-DN-0036', 'OGLE BLG-DN-0087', 'OGLE BLG-DN-0168',
                 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0233', 'OGLE BLG-DN-0286', 'OGLE BLG-DN-0305',
                 'OGLE BLG-DN-0373', 'OGLE BLG-DN-0376', 'OGLE BLG-DN-0421', 'OGLE BLG-DN-0444',
                 'OGLE BLG-DN-0531', 'OGLE BLG-DN-0588', 'OGLE BLG-DN-0595', 'OGLE BLG-DN-0690',
                 'OGLE BLG-DN-0783', 'OGLE BLG-DN-0826', 'OGLE BLG-DN-0899']


"""
lum_list= []
names = ['OGLE BLG-DN-0001_3','OGLE BLG-DN-0001_4','OGLE BLG-DN-0002_0','OGLE BLG-DN-0002_1',\
         'OGLE BLG-DN-0002_2','OGLE BLG-DN-0036_0','OGLE BLG-DN-0087_0','OGLE BLG-DN-0168_0',\
         'OGLE BLG-DN-0174_0','OGLE BLG-DN-0233_1','OGLE BLG-DN-0286_0','OGLE BLG-DN-0305_2',\
         'OGLE BLG-DN-0373_0','OGLE BLG-DN-0376_1','OGLE BLG-DN-0421_2','OGLE BLG-DN-0444_0',\
         'OGLE BLG-DN-0531_0','OGLE BLG-DN-0588_0','OGLE BLG-DN-0595_0','OGLE BLG-DN-0690_2',\
         'OGLE BLG-DN-0783_0','OGLE BLG-DN-0826_0','OGLE BLG-DN-0899_0']

direc = 'analysis_Mdot'
for name in names:
    lum_table = pd.read_csv(f"{direc}/{name}.csv")
    lum_list.append((lum_table[['t', 'L_u','L_g','L_r','L_i','L_z', 'L_y']]))
luminosity_dict = dict(zip(names, lum_list))
"""


def init_luminosity(init_data_direc):
    lum_list = []
    # names = ['OGLE BLG-DN-0001_3','OGLE BLG-DN-0001_4','OGLE BLG-DN-0002_0','OGLE BLG-DN-0002_1',\
    #          'OGLE BLG-DN-0002_2','OGLE BLG-DN-0036_0','OGLE BLG-DN-0087_0','OGLE BLG-DN-0168_0',\
    #          'OGLE BLG-DN-0174_0','OGLE BLG-DN-0233_1','OGLE BLG-DN-0286_0','OGLE BLG-DN-0305_2',\
    #          'OGLE BLG-DN-0373_0','OGLE BLG-DN-0376_1','OGLE BLG-DN-0421_2','OGLE BLG-DN-0444_0',\
    #          'OGLE BLG-DN-0531_0','OGLE BLG-DN-0588_0','OGLE BLG-DN-0595_0','OGLE BLG-DN-0690_2',\
    #          'OGLE BLG-DN-0783_0','OGLE BLG-DN-0826_0','OGLE BLG-DN-0899_0']

    direc = init_data_direc
    for name in selected_objs:
        lum_table = pd.read_csv(f"{direc}/{name}.csv")
        lum_list.append((lum_table[['t', 'L_u','L_g','L_r','L_i','L_z', 'L_y']]))
    return dict(zip(selected_objs, lum_list))


# def get_luminosity(rng, count):
#     luminosities= []
#     names = ['OGLE BLG-DN-0001_3','OGLE BLG-DN-0001_4','OGLE BLG-DN-0002_0','OGLE BLG-DN-0002_1',\
#              'OGLE BLG-DN-0002_2','OGLE BLG-DN-0036_0','OGLE BLG-DN-0087_0','OGLE BLG-DN-0168_0',\
#              'OGLE BLG-DN-0174_0','OGLE BLG-DN-0233_1','OGLE BLG-DN-0286_0','OGLE BLG-DN-0305_2',\
#              'OGLE BLG-DN-0373_0','OGLE BLG-DN-0376_1','OGLE BLG-DN-0421_2','OGLE BLG-DN-0444_0',\
#              'OGLE BLG-DN-0531_0','OGLE BLG-DN-0588_0','OGLE BLG-DN-0595_0','OGLE BLG-DN-0690_2',\
#              'OGLE BLG-DN-0783_0','OGLE BLG-DN-0826_0','OGLE BLG-DN-0899_0']
#     obj_idxs = rng.integers(low=0,high=len(names)-1, size=count)
#     direc = 'analysis_Mdot'
#     for obj_idx in obj_idxs:
#         lum_table = pd.read_csv(f"{direc}/{names[obj_idx]}.csv")
#         luminosities.append((lum_table[['t', 'L_u','L_g','L_r','L_i','L_z', 'L_y']]))
#     return luminosities


luminosity_dict = init_luminosity('Mdot_test_Jan26')


def get_luminosity(rng, count):
    luminosities = []
    OGLE_id = []
    # names = ['OGLE BLG-DN-0001_3','OGLE BLG-DN-0001_4','OGLE BLG-DN-0002_0','OGLE BLG-DN-0002_1',
    #          'OGLE BLG-DN-0002_2','OGLE BLG-DN-0036_0','OGLE BLG-DN-0087_0','OGLE BLG-DN-0168_0',
    #          'OGLE BLG-DN-0174_0','OGLE BLG-DN-0233_1','OGLE BLG-DN-0286_0','OGLE BLG-DN-0305_2',
    #          'OGLE BLG-DN-0373_0','OGLE BLG-DN-0376_1','OGLE BLG-DN-0421_2','OGLE BLG-DN-0444_0',
    #          'OGLE BLG-DN-0531_0','OGLE BLG-DN-0588_0','OGLE BLG-DN-0595_0','OGLE BLG-DN-0690_2',
    #          'OGLE BLG-DN-0783_0','OGLE BLG-DN-0826_0','OGLE BLG-DN-0899_0']
    obj_idxs = rng.integers(low=0,high=len(selected_objs)-1, size=count)
    if count > 1:
        for obj_idx in obj_idxs:
            luminosities.append(luminosity_dict[selected_objs[obj_idx]])
            OGLE_id.append(selected_objs[obj_idx][12:16])
    else:
        luminosities = luminosity_dict[selected_objs[obj_idxs[0]]]
        OGLE_id = selected_objs[obj_idxs[0]][12:16]
    return luminosities, OGLE_id


MWDENSITY = MWDensity()


def get_coordinates(rng, count):
    mw_coords = MWDENSITY.sample_eq(shape=count, rng=rng)
    return mw_coords


def get_inclination(rng, count):
    return rng.uniform(low=0, high=89, size=count)


def get_inclination_single(rng, count=1):
    inclins = rng.uniform(low=0, high=89, size=count)
    return(inclins[0])


def get_flux(L, d, i):
    return L * np.cos(i * np.pi/180)/(2*np.pi * d**2)


LSST_A_TO_EBV = {
    'u': 4.145,
    'g': 3.237,
    'r': 2.273,
    'i': 1.684,
    'z': 1.323,
    'y': 1.088,
}


def hist_mpeak(ms, m_max):
    # passbands = ['u', 'g', 'r', 'i', 'z', 'y']
    fig, axs = plt.subplots(2, 3, figsize=(24,16))
    m_peaks = []

    for passband in passbands:
        m_peak_single_passband = []
        for m in ms:
            m_peak_single_passband.append(np.min(m[passband]))
        m_peaks.append(np.array(m_peak_single_passband))
    m_peak_dict = dict(zip(passbands, m_peaks))

    fig.suptitle("Histogram for Peak Magnitude in each passband")
    for ax, m_peak, passband in zip(axs.reshape(-1),m_peaks,passbands):
        ax.hist(m_peak[m_peak<m_max], bins=50)
        ax.set_title(f'passband_{passband}')
        ax.set_xlabel(f'magnitude')
        ax.set_xlim(5,99)

    direc = 'pictures_generated'
    os.makedirs(direc, exist_ok=True)
    plt.savefig(os.path.join(direc, f'Peak_magnitude_hist.png'))
    plt.close()

    return m_peak_dict


def write_header(f, event_num):
    # with open(file_name, 'w') as f:
#         f.write(f"""SURVEY: LSST
# FILTERS: ugrizY
# MODEL: m-Dwarf-Flare-Model
# RECUR_TYPE: NON-RECUR
# MODEL_PARNAMES: OGLE_ID,start_time,end_time,distance,inclination.
# NEVENT: {event_num}
#
# DOCUMENTATION:
#   PURPOSE: Supernovae outburst ligthtcurve using OGLE data and estimated distances from Gaia
#   REF:
#   - AUTHOR: Qifeng Cheng
#   USAGE_KEY: GENMODEL
#   NOTES:
#   - Lightcurve instances were taken from OGLE
#   - Distance data was taken from Gaia
#   - Extinction data was taken from SFD, Bayestar
#   PARAMS:
#   - OGLE_ID - OGLE Dwarf Nova Catalog object ID
#   - start_time - Start time of the reference outburst (in HJD-2450000)
#   - end_time - End time of the reference outburst (in HJD-2450000)
#   - distance - Distance to the supernovae (in pc)
#   - inclination - Inclination of the observation (in degree)
# DOCUMENTATION_END:
#
# #------------------------------
# """
#         )

    f.write(f"""DOCUMENTATION:
  PURPOSE: Supernovae outburst ligthtcurve using OGLE data and estimated distances from Gaia
  REF:
  - AUTHOR: Qifeng Cheng
  USAGE_KEY: GENMODEL
  NOTES:
  - Lightcurve instances were taken from OGLE
  - Distance data was taken from Gaia
  - Extinction data was taken from SFD, Bayestar
  PARAMS:
  - OGLE_ID - OGLE Dwarf Nova Catalog object ID
  - start_time - Start time of the reference outburst (in HJD-2450000)
  - end_time - End time of the reference outburst (in HJD-2450000)
  - distance - Distance to the supernovae (in pc)
  - inclination - Inclination of the observation (in degree)
DOCUMENTATION_END:

SURVEY: LSST
FILTERS: ugrizY
MODEL: m-Dwarf-Flare-Model
RECUR_TYPE: NON-RECUR
MODEL_PARNAMES: OGLE_ID,start_time,end_time,distance,inclination.
NEVENT: {event_num}

#------------------------------
"""
        )

def generate_outburst(start_index, event_num, file_name):

    i_event = start_index
    rng = np.random.default_rng(start_index)

    with open(file_name, 'w') as f:
        # write_header(f, event_num)

        mag_es_all = []
        mag_es_lclib = []

        while i_event < event_num + start_index:
            l, OGLE_id = get_luminosity(rng, 1)
            coord = get_coordinates(rng, 1)[0]
            i = get_inclination_single(rng, 1)
            extin = get_extinction(coord.ra.deg, coord.dec.deg, coord.distance.pc, cache_dir=None)

            distance_pc = coord.distance.to(u.pc).value
            distance_cm = coord.distance.to(u.cm).value
            ra = coord.ra.deg
            dec = coord.dec.deg
            coord_gal = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs').galactic

            l = l.copy()
            l.rename(columns={'L_u': 'u', 'L_g': 'g', 'L_r': 'r', 'L_i': 'i', 'L_z': 'z', 'L_y': 'y'}, inplace=True)
            l_no_t = l.drop(axis=1, labels='t')

            flux = get_flux(L=l_no_t, d=distance_cm, i=i)
            mag_noe = -2.5 * np.log10(flux / 3.63e-20)
            l.update(mag_noe)
            mag_e = l.copy()

            for passband in passbands:
                mag_e[passband] = mag_noe[passband] + extin * LSST_A_TO_EBV[passband]
            mag_es_all.append(mag_e)

            if any([np.any((mag_e[passband] > 99.0)) or np.any((mag_e[passband] < 5.0)) for passband in passbands]):
                continue

            anglematch_b = max(5, 0.5 * np.abs(coord_gal.b.deg))

            f.write(
                f'START_EVENT: {i_event}\n'
                f'NROW: {len(mag_e)+1} l: {coord_gal.l.value:.5f} b: {coord_gal.b.value:.5f}\n'
                f'PARVAL: {int(OGLE_id)},{mag_e["t"][0]},{mag_e["t"][len(mag_e) - 1]},{distance_pc:.2f},{i}\n'
                f'ANGLEMATCH_b: {anglematch_b:.1f}\n'
            )

            time = mag_e.loc[0]['t']
            f.write(f'T: {time:7.4f}')

            for passband in passbands:
                f.write(f' {mag_e.loc[0][passband]:.3f}')
            f.write(f'\n')

            for i_row in range(1, len(mag_e)):
                time = mag_e.loc[i_row]['t']
                f.write(f'S: {time:7.4f}')

                for passband in passbands:
                    f.write(f' {mag_e.loc[i_row][passband]:.3f}')
                f.write(f'\n')

            time = mag_e.loc[len(mag_e) - 1]['t'] + 0.01
            f.write(f'S: {time:7.4f}')

            for passband in passbands:
                f.write(f' {mag_e.loc[0][passband]:.3f}')
            f.write(f'\n')
            f.write(
                f'END_EVENT: {i_event}\n'
                '\n'
            )
            mag_es_lclib.append(mag_e)
            i_event = i_event + 1


def main():

    # directory = 'stitch_file'
    # os.makedirs(directory, exist_ok=True)
    #
    # header_filename = os.path.join(directory, f'header.txt')
    # gen_1 = os.path.join(directory, f'objects_1.txt')
    # gen_2 = os.path.join(directory, f'objects_2.txt')
    #
    # write_header(event_num=100, file_name=header_filename)
    # generate_outburst(start_index=0, event_num=50, file_name=gen_1)
    # generate_outburst(start_index=50, event_num=50, file_name=gen_2)

    # passbands = ['u', 'g', 'r', 'i', 'z', 'y']

    start_index = 42
    event_num = 100
    i_event = 0
    rng = np.random.default_rng(start_index)

    with open('LCLIB_test_jan_26.txt', 'w') as f:
        write_header(f, event_num)

        mag_es_all = []
        mag_es_lclib = []

        while i_event < event_num:
            l, OGLE_id = get_luminosity(rng, 1)
            coord = get_coordinates(rng, 1)[0]
            i = get_inclination_single(rng, 1)
            extin = get_extinction(coord.ra.deg, coord.dec.deg, coord.distance.pc, cache_dir=None)

            distance_pc = coord.distance.to(u.pc).value
            distance_cm = coord.distance.to(u.cm).value
            ra = coord.ra.deg
            dec = coord.dec.deg
            coord_gal = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs').galactic

            l = l.copy()
            l.rename(columns={'L_u': 'u', 'L_g': 'g', 'L_r': 'r', 'L_i': 'i', 'L_z': 'z', 'L_y': 'y'}, inplace=True)
            l_no_t = l.drop(axis=1, labels='t')

            flux = get_flux(L=l_no_t, d=distance_cm, i=i)
            mag_noe = -2.5 * np.log10(flux / 3.63e-20)
            l.update(mag_noe)
            mag_e = l.copy()

            for passband in passbands:
                mag_e[passband] = mag_noe[passband] + extin * LSST_A_TO_EBV[passband]
            mag_es_all.append(mag_e)

            if any([np.any((mag_e[passband] > 99.0)) or np.any((mag_e[passband] < 5.0)) for passband in passbands]):
                continue

            anglematch_b = max(5, 0.5 * np.abs(coord_gal.b.deg))

            f.write(
                f'START_EVENT: {i_event}\n'
                f'NROW: {len(mag_e)+1} l: {coord_gal.l.value:.5f} b: {coord_gal.b.value:.5f}\n'
                f'PARVAL: {int(OGLE_id)},{mag_e["t"][0]},{mag_e["t"][len(mag_e) - 1]},{distance_pc:.2f},{i}\n'
                f'ANGLEMATCH_b: {anglematch_b:.1f}\n'
            )

            time = mag_e.loc[0]['t']
            f.write(f'T: {time:7.4f}')

            for passband in passbands:
                f.write(f' {mag_e.loc[0][passband]:.3f}')
            f.write(f'\n')

            for i_row in range(1, len(mag_e)):
                time = mag_e.loc[i_row]['t']
                f.write(f'S: {time:7.4f}')

                for passband in passbands:
                    f.write(f' {mag_e.loc[i_row][passband]:.3f}')
                f.write(f'\n')

            time = mag_e.loc[len(mag_e) - 1]['t'] + 0.01
            f.write(f'S: {time:7.4f}')

            for passband in passbands:
                f.write(f' {mag_e.loc[0][passband]:.3f}')
            f.write(f'\n')
            f.write(
                f'END_EVENT: {i_event}\n'
                '\n'
            )
            mag_es_lclib.append(mag_e)
            i_event = i_event + 1
    mpeak_dict = hist_mpeak(mag_es_lclib, 30)


if __name__ == '__main__':
    main()