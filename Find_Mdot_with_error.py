from scipy import integrate
import astropy.constants as cons
import astropy.units as u
from scipy.optimize import root
from functools import partial

from cand_all import data_loading, find_cand, plot_cand
import numpy as np
import matplotlib.pyplot as plt
import os
from astropy.io import ascii
from astropy.table import Table
from light_curve import BazinFit

h = cons.h.cgs.value
G = cons.G.cgs.value
sigma_sb = cons.sigma_sb.cgs.value
k_B = cons.k_B.cgs.value
R_sun = cons.R_sun.cgs.value
c = cons.c.cgs.value
M_sun = cons.M_sun.cgs.value
one_parsec = 1.0 * u.parsec
parsec_to_cm = one_parsec.to(u.cm).value
yr_to_sec = (1.0 * u.yr).to(u.s).value
wavelengths = {'V': 5400e-8,
         'u': 3694.25e-8,
         'g': 4840.83e-8,
         'r': 6257.74e-8,
         'i': 7560.48e-8,
         'z': 8701.29e-8,
         'y': 9749.32e-8
         }



def T_0(M, Mdot, r_in):
    return (2 ** (3 / 4) * (3 / 7) ** (7 / 4) * (G * M * Mdot / (np.pi * sigma_sb * r_in ** 3)) ** (1 / 4))


def r_0(r_in):
    return ((7 / 6) ** 2 * r_in)


def x_fun(nu, T0, r, r0):
    return (h * nu / (k_B * T0) * (r / r0) ** (3 / 4))


def f(Mdot, *, fv, nu, rin, rout, i, d, M):
    T0 = T_0(M, Mdot, rin)
    r0 = r_0(rin)
    xin = x_fun(nu, T0, rin, r0)
    xout = x_fun(nu, T0, rout, r0)
    fun_integr = lambda x: (x ** (5 / 3)) / (np.exp(x) - 1)
    integ, inte_err = integrate.quad(fun_integr, xin, xout)

    return (fv - (16 * np.pi) / (3 * d ** 2) * np.cos(i) * (k_B * T0 / h) ** (8 / 3) * h * (nu ** (1 / 3)) / (
                c ** 2) * (r0 ** 2) * integ)


def find_Mdot(fv, lamb, rin, rout, i, d, M, Mdot_previous):
    return (root(partial(f, fv=fv, nu=c / (lamb), rin=rin, rout=rout, i=i, d=d, M=M), Mdot_previous))


def find_flux(Mdot, nu, rin, rout, i, d, M):
    T0 = T_0(M, Mdot, rin)
    r0 = r_0(rin)
    xin = x_fun(nu, T0, rin, r0)
    xout = x_fun(nu, T0, rout, r0)
    fun_integr = lambda x: (x ** (5 / 3)) / (np.exp(x) - 1)
    integ, inte_err = integrate.quad(fun_integr, xin, xout)
    return ((16 * np.pi) / (3 * d ** 2) * np.cos(i) * (k_B * T0 / h) ** (8 / 3) * h * (nu ** (1 / 3)) / (c ** 2) * (
                r0 ** 2) * integ)


def norm_flux(F_old, F_peak):
    return (F_old * F_peak / np.max(F_old))


# basic initialization

data = ascii.read("analysis/fitting_info.csv")
dis_extin = ascii.read("analysis/extinction_1arcsec.csv")


distances_parsec = dict(zip(dis_extin['name'], dis_extin['distance']))
extinction_baye = dict(zip(dis_extin['name'], dis_extin['extinction_baye']))
extinction_sfd = dict(zip(dis_extin['name'], dis_extin['extinction_sfd']))
t_st_dic = dict(zip(data['name'], data['t_start']))
t_ed_dic = dict(zip(data['name'], data['t_end']))
outburst_idx = dict(zip(data['name'], data['outburst_index']))


def load_candidates_with_id(obj_name, data_direc):
    k = obj_name[-4:]
    file_name = f"{data_direc}/OGLE-BLG-DN-{k}.dat"
    # obj_name = "OGLE BLG-DN-" + k
    t, m, err = data_loading(file_name)
    t_st = t_st_dic[f'{obj_name}']
    t_ed = t_ed_dic[f'{obj_name}']
    idx = (t >= t_st) & (t <= t_ed)
    cand_t = t[idx]
    cand_m = m[idx]
    cand_m_err = err[idx]
    return cand_t, cand_m, cand_m_err


def mag_to_flux(mag, with_coef, mag_err=0.01):

    coef = 1
    if with_coef:
        coef = 3.63e-20
    # coef = (coef-1.) * with_coef + 1.
    # print(f'coef \n {coef}')
    flux = coef * 10 ** (-0.4 * mag)
    flux_err = 0.4 * np.log(10) * flux * mag_err
    return flux, flux_err


def flux_to_mag(flux, with_coef):
    coef = 1
    if with_coef:
        coef = 3.63e-20
    # coef = (coef-1.) * with_coef + 1.
    # print(f'coef \n {coef}')
    # print(f'coef \n {coef}')
    mag = np.log10(flux / coef) / (-0.4)
    return mag

def fit_light_curve(t, flux, flux_err):

    fit = BazinFit('mcmc-lmsder', mcmc_niter=100_000, lmsder_niter=20)
    parameters = fit(t, flux, flux_err)
    # print(f'parameters: \n {parameters}')
    delt_t = 0.1  # days
    t_start = t.min()
    t_end = t.max() + 2 * delt_t
    t_model = np.arange(t_start, t_end, delt_t)
    flux_model = fit.model(np.arange(t_start, t_end, delt_t), parameters)
    return t_model, flux_model


def plot_light_curve_model(t_data, m_data, t_model, m_model, obj_name):
    plt.figure()
    plt.scatter(t_data, m_data, label='original data')
    plt.plot(t_model, m_model, label='model data')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.xlabel('time')
    plt.ylabel('magnitude')
    plt.title(obj_name)

    direc = 'pictures_fittings_jan26'
    os.makedirs(direc, exist_ok=True)
    plt.savefig(os.path.join(direc, f'{obj_name}.png'))
    plt.close()


def find_Mdot_lum_table(obj_name, t_model, flux_model, mag_model, if_norm, mass=M_sun, MD_pre=1e18, conversion_coef=True):

    Mdot_peak = 4e-8 * mass / yr_to_sec  #question: mass or M_sun


    if obj_name in distances_parsec:  # ones that have distances
        print(f'Object {obj_name} have distance data from GAIA')

        distance = distances_parsec[obj_name] * parsec_to_cm
        # mass = M_sun
        # MD_pre = 1e18



        if np.isnan(extinction_baye[obj_name]):
            extinction = extinction_sfd[obj_name]
        else:
            extinction = extinction_baye[obj_name]
        mag_model_extins = mag_model - extinction
        fv_model_extins = mag_to_flux(mag_model_extins, with_coef=conversion_coef)[0]

        MD_al = np.zeros(len(fv_model_extins))

        # fv_model_extins = 3.63e-20 * 10 ** (-0.4 * mag_model_extins)

        # for fv_modle_extin in fv_model_extins:
        #     #                 fv = 3.63e-20 * 10 ** (-0.4 * m_model_extin)
        #     MD = find_Mdot(fv=fv_modle_extin, lamb=wavelengths['V'], rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=distance,
        #                    M=mass, Mdot_previous=MD_pre)
        #     MD_al.append(MD['x'][0])
        #     MD_pre = MD['x'][0]
        for flux_id in np.argsort(fv_model_extins)[::-1]:
            MD = find_Mdot(fv=fv_model_extins[flux_id], lamb=wavelengths['V'], rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4,
                           d=distance, M=mass, Mdot_previous=MD_pre)
            MD_al[flux_id] = MD['x'][0]
            MD_pre = MD['x'][0]

        if if_norm:
            if np.max(MD_al) > Mdot_peak:
                print('Mdot exceeds Mdot_peak, start to normalize...')
                flux_peak = find_flux(Mdot_peak, nu=c / wavelengths['V'], rin=0.0025 * R_sun, rout=10 ** 11, i=np.pi / 4,
                                      d=distance, M=M_sun)
                flux_norm = norm_flux(F_old=fv_model_extins, F_peak=flux_peak)
                print(f'old Mdot: Mdot_old = {MD_al}')

                MD_al = []
                # MD_pre = 1e18
                for f_norm in flux_norm:
                    #                 fv = 3.63e-20 * 10 ** (-0.4 * m_model_extin)
                    MD = find_Mdot(fv=f_norm, lamb=wavelengths['V'], rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=distance,
                                   M=mass, Mdot_previous=MD_pre)
                    MD_al.append(MD['x'][0])
                    MD_pre = MD['x'][0]
                print(f'after normalization: Mdot_norm = {MD_al}')

    else:  # the ones that don't have distances/
        print(f'No distance data found for {obj_name}, assume d = 1kpc and normalize...')

        distance = 1e3 * parsec_to_cm
        # mass = M_sun
        # MD_pre = 1e18
        # MD_al = []
        # mag_model_extins = mag_model
        flux_peak = find_flux(Mdot_peak, nu=c / 540e-7, rin=0.0025 * R_sun, rout=10 ** 11, i=np.pi / 4,
                              d=1e3 * parsec_to_cm, M=M_sun)


        flux_norm = norm_flux(F_old=flux_model, F_peak=flux_peak)
        MD_al = np.zeros(len(flux_norm))

        for flux_id in np.argsort(flux_norm)[::-1]:
            MD = find_Mdot(fv=flux_norm[flux_id], lamb=wavelengths['V'], rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4,
                           d=distance, M=mass, Mdot_previous=MD_pre)
            MD_al[flux_id] = MD['x'][0]
            MD_pre = MD['x'][0]

        # for f_norm in flux_norm:
        #     #                 fv = 3.63e-20 * 10 ** (-0.4 * m_model_extin)
        #     MD = find_Mdot(fv=f_norm, lamb=540e-7, rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=distance, M=mass,
        #                    Mdot_previous=MD_pre)
        #     MD_al.append(MD['x'][0])
        #     MD_pre = MD['x'][0]


    table = Table()
    table['t'] = t_model
    table['Mdot'] = MD_al

    # find luminosity at other passbands

        # lambs = [3694.25e-8, 4840.83e-8, 6257.74e-8, 7560.48e-8, 8701.29e-8, 9749.32e-8]  # u, g, r, i, z, y
        # passband_name = ['L_u', 'L_g', 'L_r', 'L_i', 'L_z', 'L_y']
    lambs = [5400e-8, 3694.25e-8, 4840.83e-8, 6257.74e-8, 7560.48e-8, 8701.29e-8, 9749.32e-8]  # V, u, g, r, i, z, y
    passband_name = ['L_V', 'L_u', 'L_g', 'L_r', 'L_i', 'L_z', 'L_y']

    for lam_passband, name_passband in zip(lambs, passband_name):
        flux_derived = []
        for Md in MD_al:
            flux_derived.append(
                find_flux(Mdot=Md, nu=c / lam_passband, rin=0.0025 * R_sun, rout=10 ** 11, i=0, d=1, M=M_sun))

        flux_derived = np.array(flux_derived)
        luminosity_derived = 2 * np.pi * flux_derived
        table[f'{name_passband}'] = luminosity_derived

    return table

def save_Mdot_lum_data(direc, table, obj_name):
    # direc = 'analysis_Mdot_test'
    os.makedirs(direc, exist_ok=True)
    # outbursts = data['outburst_index']
    # for table_out, obj_name, outburst_i in zip(Mdot_lum_tables, names, outbursts):
    table.write(os.path.join(direc, f'{obj_name}.csv'), format = 'ascii.csv', overwrite=True)




def main():
    # selected_obj = ['OGLE BLG-DN-0036', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0826']
    # selected_obj = ['OGLE BLG-DN-0174']
    # selected_obj = ['OGLE BLG-DN-0001_3', 'OGLE BLG-DN-0001_4', 'OGLE BLG-DN-0002_0', 'OGLE BLG-DN-0002_1',
    #          'OGLE BLG-DN-0002_2','OGLE BLG-DN-0036_0','OGLE BLG-DN-0087_0','OGLE BLG-DN-0168_0',
    #          'OGLE BLG-DN-0174_0','OGLE BLG-DN-0233_1','OGLE BLG-DN-0286_0','OGLE BLG-DN-0305_2',
    #          'OGLE BLG-DN-0373_0','OGLE BLG-DN-0376_1','OGLE BLG-DN-0421_2','OGLE BLG-DN-0444_0',
    #          'OGLE BLG-DN-0531_0','OGLE BLG-DN-0588_0','OGLE BLG-DN-0595_0','OGLE BLG-DN-0690_2',
    #          'OGLE BLG-DN-0783_0','OGLE BLG-DN-0826_0','OGLE BLG-DN-0899_0']
    # selected_obj = ['OGLE BLG-DN-0001', 'OGLE BLG-DN-0001', 'OGLE BLG-DN-0002', 'OGLE BLG-DN-0002',
    #          'OGLE BLG-DN-0002','OGLE BLG-DN-0036','OGLE BLG-DN-0087','OGLE BLG-DN-0168',
    #          'OGLE BLG-DN-0174','OGLE BLG-DN-0233','OGLE BLG-DN-0286','OGLE BLG-DN-0305',
    #          'OGLE BLG-DN-0373','OGLE BLG-DN-0376','OGLE BLG-DN-0421','OGLE BLG-DN-0444',
    #          'OGLE BLG-DN-0531','OGLE BLG-DN-0588','OGLE BLG-DN-0595','OGLE BLG-DN-0690',
    #          'OGLE BLG-DN-0783','OGLE BLG-DN-0826','OGLE BLG-DN-0899']
    selected_obj = ['OGLE BLG-DN-0001', 'OGLE BLG-DN-0001', 'OGLE BLG-DN-0002', 'OGLE BLG-DN-0002',
                    'OGLE BLG-DN-0002', 'OGLE BLG-DN-0036', 'OGLE BLG-DN-0087', 'OGLE BLG-DN-0168',
                    'OGLE BLG-DN-0174', 'OGLE BLG-DN-0233', 'OGLE BLG-DN-0286', 'OGLE BLG-DN-0305',
                    'OGLE BLG-DN-0373', 'OGLE BLG-DN-0421', 'OGLE BLG-DN-0444',
                    'OGLE BLG-DN-0531', 'OGLE BLG-DN-0588', 'OGLE BLG-DN-0595', 'OGLE BLG-DN-0690',
                    'OGLE BLG-DN-0783', 'OGLE BLG-DN-0826', 'OGLE BLG-DN-0899']
    for name in selected_obj:
        cand_t, cand_mag, cand_mag_err = load_candidates_with_id(name, 'phot')
        # print(f'obj:{name}\n cand_t: \n{cand_t} \n cand_m: \n{cand_m} \n cand_m_err:\n{cand_m_err}')
        cand_flux, cand_flux_err = mag_to_flux(cand_mag, True, cand_mag_err)
        # print(f'{cand_flux}, {cand_flux_err}')
        t_model, flux_model = fit_light_curve(cand_t, cand_flux, cand_flux_err)

        mag_model = flux_to_mag(flux_model, True)

        plot_light_curve_model(cand_t, cand_mag, t_model, mag_model, name)

        table = find_Mdot_lum_table(name, t_model, flux_model,mag_model, if_norm=False, conversion_coef=True, MD_pre=1e13)
        save_Mdot_lum_data(direc='Mdot_test_Jan28', table=table, obj_name=name)

if __name__ == '__main__':
    main()