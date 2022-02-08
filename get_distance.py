from astropy.coordinates import SkyCoord
from pyvo.dal import TAPService
from astropy.table import QTable
import numpy as np
import os


def get_distance(obj_names, search_radius_arcsec, save_dir, file_name, print_position=False):
    table = QTable(names=('obj_name', 'd', 'r_med_geo', 'r_lo_geo' , 'r_hi_geo', 'r_med_photogeo', 'r_lo_photogeo', 'r_hi_photogeo', 'phot_g_mean_mag'),
                   dtype=('str', 'float64', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32', 'float32'))

    for name in obj_names:
        position = SkyCoord.from_name(name)
        ra = position.ra
        ra_deg = ra.deg
        dec = position.dec
        dec_deg = dec.deg

        string_query = f'''SELECT distance(ra, dec, {ra_deg}, {dec_deg}) as d, r_med_geo, r_lo_geo, r_hi_geo, r_med_photogeo,r_lo_photogeo, r_hi_photogeo, phot_g_mean_mag FROM gedr3dist.main JOIN gaia.edr3lite USING (source_id) WHERE distance(ra, dec, {ra_deg}, {dec_deg}) < {search_radius_arcsec / 3600.0}'''
        tap = TAPService('https://dc.zah.uni-heidelberg.de/tap')
        response = tap.search(string_query)
        r = response.to_table()

        if len(r) > 1:
            index = np.where(r['r_med_geo'] == np.amin(r['r_med_geo']))
            table.add_row((name, r[index]['d'], r[index]['r_med_geo'], r[index]['r_lo_geo'], r[index]['r_hi_geo'], r[index]['r_med_photogeo'], r[index]['r_lo_photogeo'], r[index]['r_hi_photogeo'], r[index]['phot_g_mean_mag']))
        elif len(r) == 0:
            continue
        else:
            table.add_row((name, r[0]['d'], r[0]['r_med_geo'], r[0]['r_lo_geo'], r[0]['r_hi_geo'], r[0]['r_med_photogeo'], r[0]['r_lo_photogeo'], r[0]['r_hi_photogeo'], r[0]['phot_g_mean_mag']))
        if print_position:
            print(position)

    os.makedirs(save_dir, exist_ok=True)
    table.write(os.path.join(save_dir , f'{file_name}.csv'), format='ascii.csv', overwrite=True)


def main():
    selected_obj = ['OGLE BLG-DN-0001', 'OGLE BLG-DN-0002', 'OGLE BLG-DN-0036', 'OGLE BLG-DN-0087',
                    'OGLE BLG-DN-0168', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0233', 'OGLE BLG-DN-0275',
                    'OGLE BLG-DN-0286', 'OGLE BLG-DN-0305', 'OGLE BLG-DN-0373', 'OGLE BLG-DN-0376',
                    'OGLE BLG-DN-0421', 'OGLE BLG-DN-0444', 'OGLE BLG-DN-0458', 'OGLE BLG-DN-0531',
                    'OGLE BLG-DN-0588', 'OGLE BLG-DN-0595', 'OGLE BLG-DN-0690', 'OGLE BLG-DN-0783',
                    'OGLE BLG-DN-0826', 'OGLE BLG-DN-0899']
    direc = 'analysis'
    file_name = 'distance_1arcsec'
    search_radius_arcsec = 1
    get_distance(selected_obj, search_radius_arcsec, direc, file_name, print_position=False)


if __name__ == '__main__':
    main()
