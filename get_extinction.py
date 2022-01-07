from astropy.table import Table
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import units as u
from dustmaps.bayestar import BayestarWebQuery
from dustmaps.sfd import SFDWebQuery
from dustmaps.marshall import MarshallQuery
import os

def main():
    data = ascii.read("analysis/distance_1arcsec.csv")
    selected_obj = data['obj_name']
    ra_deg = []
    dec_deg = []

    for name in selected_obj:
        position = SkyCoord.from_name(name)
        ra = position.ra
        ra_deg.append(ra.deg * u.deg)
        dec = position.dec
        dec_deg.append(dec.deg * u.deg)

    coef = 2.742
    distance = data['r_med_geo']
    d = []
    for dis in distance:
        d.append(dis * u.parsec)

    q = BayestarWebQuery(version='bayestar2019')
    sfd = SFDWebQuery()
    E_baye = []
    for i in range(len(d)):
        c_icrs = SkyCoord(ra=ra_deg[i], dec=dec_deg[i], distance=d[i], frame='icrs')
        c_gla = c_icrs.galactic
        E_baye.append(coef * q(c_gla, mode='median'))

    coords = SkyCoord(ra_deg, dec_deg, frame='icrs')
    E_sfd = coef * sfd(coords)

    coef = 2.742
    distance = data['r_med_geo']
    d = []
    for dis in distance:
        d.append(dis * u.parsec)

    q2 = MarshallQuery()
    E_marsh = []
    for i in range(len(d)):
        c_icrs = SkyCoord(ra=ra_deg[i], dec=dec_deg[i], distance=d[i], frame='icrs')
        c_gla = c_icrs.galactic
        E_marsh.append(coef * q2(c_gla))

    t_out = Table()
    t_out['name'] = selected_obj
    t_out['l, b'] = coords.galactic
    t_out['distance'] = distance
    t_out['extinction_baye'] = E_baye
    t_out['extinction_sfd'] = E_sfd
    t_out['extinction_marshall'] = E_marsh

    direc = 'analysis'
    os.makedirs(direc, exist_ok=True)
    t_out.write(os.path.join(direc, f'extinction_1arcsec.csv'), format = 'ascii.csv', overwrite=True)

if __name__ == '__main__':
    main()