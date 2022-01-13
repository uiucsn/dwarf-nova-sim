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

"""
def main():
    MD_pre = 1e18
    m_input = 18.6
    flux_input = 3.63e-20 * 10 ** (-0.4 * m_input)

    Mdot = find_Mdot(fv=flux_input, lamb=5400e-8, rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=5e3 * parsec_to_cm,
                                   M=M_sun, Mdot_previous=MD_pre)
    flux_output = find_flux(Mdot=Mdot['x'][0], nu=c/5400e-8, rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=5e3 * parsec_to_cm,
                    M=M_sun)

    m_output = np.log10(flux_output/ 3.63e-20) / (-0.4)
    luminosity_derived = 2 * np.pi * flux_output
    print(f'input m: {m_input}')
    print(f'input flux: {flux_input}')
    print(f"derived Mdot: {Mdot['x'][0]}")
    print(f'derived flux: {flux_output}')
    print(f'derived m: {m_output}')



    # print(4e-8 * M_sun / yr_to_sec)

"""


#### find Mdot for these three
def main():
    data = ascii.read("analysis/fitting_info.csv")
    dis_extin = ascii.read("analysis/extinction_1arcsec.csv")


    fit_num = data['method_num']
    fit = BazinFit('mcmc-lmsder', mcmc_niter=100_000)
    # fit = BazinFit('mcmc')
    selected_obj = ['OGLE BLG-DN-0036', 'OGLE BLG-DN-0174', 'OGLE BLG-DN-0826']
    # c_i = data['outburst_index']

    # t_st_al = data['t_start']
    # t_ed_al = data['t_end']

    distances_parsec = dict(zip(dis_extin['name'], dis_extin['distance']))
    extinction_baye = dict(zip(dis_extin['name'], dis_extin['extinction_baye']))
    extinction_sfd = dict(zip(dis_extin['name'], dis_extin['extinction_sfd']))
    t_st_dic = dict(zip(data['name'], data['t_start']))
    t_ed_dic = dict(zip(data['name'], data['t_end']))
    outburst_idx = dict(zip(data['name'], data['outburst_index']))

    Mdot_tables = []

    for name in selected_obj:
        k = name[-4:]
        print(f'object number: {k}')
        file_name = "phot/OGLE-BLG-DN-" + k + ".dat"
        obj_name = "OGLE BLG-DN-" + k
        t, m = data_loading(file_name)
        t_st = t_st_dic[f'{obj_name}']
        t_ed = t_ed_dic[f'{obj_name}']
        print(f'starting time:{t_st}')
        print(f'ending time:{t_ed}')
        idx = (t >= t_st) & (t <= t_ed)
        cand_t = t[idx]
        cand_m = m[idx]
        print(cand_t)
        print(cand_m)
        cand_flux = 3.63e-20 * 10 ** (-0.4 * cand_m)
        try:
            parameters = fit(cand_t, cand_flux)
            delt_t = 0.1  # days
            t_start = cand_t.min()
            t_end = cand_t.max() + 2 * delt_t
            t_model = np.arange(t_start, t_end, delt_t)
            flux_model = fit.model(np.arange(t_start, t_end, delt_t), parameters)
            m_model = np.log10(flux_model / 3.63e-20) / (-0.4)
    #plot modeled light_curve
            #         fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize=(6, 14))
            #         fig.suptitle(f'{obj_name}_{c_i}')
            #         ax1.plot(cand_t, cand_flux, 'x')
            #         ax1.plot(t_model, flux_model)
            #         ax1.set_ylabel('observed flux with model')

    #finding Mdot under two different situations
            Mdot_peak = 4e-8 * M_sun / (yr_to_sec)

            if obj_name in distances_parsec:  # ones that have distances
                print('have distance', obj_name)
                distance = distances_parsec[obj_name] * parsec_to_cm
                mass = M_sun
                MD_pre = 1e18
                MD_al = []

                if np.isnan(extinction_baye[obj_name]):
                    extinction = extinction_sfd[obj_name]
                else:
                    extinction = extinction_baye[obj_name]
                m_model_extins = m_model - extinction
                fv_model_extins = 3.63e-20 * 10 ** (-0.4 * m_model_extins)

                for fv_modle_extin in fv_model_extins:
                    #                 fv = 3.63e-20 * 10 ** (-0.4 * m_model_extin)
                    MD = find_Mdot(fv=fv_modle_extin, lamb=540e-7, rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=distance,
                                   M=mass, Mdot_previous=MD_pre)
                    MD_al.append(MD['x'][0])
                    MD_pre = MD['x'][0]

                if np.max(MD_al) > Mdot_peak:
                    print('max MD_al larger than Mdot_peak, start to normalize...')
                    flux_peak = find_flux(Mdot_peak, nu=c / 540e-7, rin=0.0025 * R_sun, rout=10 ** 11, i=np.pi / 4,
                                          d=distance, M=M_sun)
                    flux_norm = fv_model_extins
                    # flux_norm = norm_flux(F_old=fv_model_extins, F_peak=flux_peak)
                    print(f'Mdot_old = {MD_al}')

                    MD_al = []
                    MD_pre = 1e18
                    for f_norm in flux_norm:
                        #                 fv = 3.63e-20 * 10 ** (-0.4 * m_model_extin)
                        MD = find_Mdot(fv=f_norm, lamb=540e-7, rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=distance,
                                       M=mass, Mdot_previous=MD_pre)
                        MD_al.append(MD['x'][0])
                        MD_pre = MD['x'][0]
                    print(f'Normalize finished')
                    print(f'Mdot_norm = {MD_al}')
            else:  # the ones that don't have distances/
                print('method2_1kpc', obj_name)
                distance = 1e3 * parsec_to_cm
                mass = M_sun
                MD_pre = 1e18
                MD_al = []
                m_model_extins = m_model
                flux_peak = find_flux(Mdot_peak, nu=c / 540e-7, rin=0.0025 * R_sun, rout=10 ** 11, i=np.pi / 4,
                                      d=1e3 * parsec_to_cm, M=M_sun)
                flux_norm = norm_flux(F_old=flux_model, F_peak=flux_peak)

                for f_norm in flux_norm:
                    #                 fv = 3.63e-20 * 10 ** (-0.4 * m_model_extin)
                    MD = find_Mdot(fv=f_norm, lamb=540e-7, rin=0.0025 * R_sun, rout=1e11, i=np.pi / 4, d=distance, M=mass,
                                   Mdot_previous=MD_pre)
                    MD_al.append(MD['x'][0])
                    MD_pre = MD['x'][0]

            t = Table()
            t['t'] = t_model
            t['Mdot'] = MD_al

    # find luminosity at other passbands

            lambs = [5400e-8, 3694.25e-8, 4840.83e-8, 6257.74e-8, 7560.48e-8, 8701.29e-8, 9749.32e-8]  # u, g, r, i, z, y
            passband_name = ['L_V', 'L_u', 'L_g', 'L_r', 'L_i', 'L_z', 'L_y']

            for lam_passband, name_passband in zip(lambs, passband_name):
                flux_derived = []
                for Md in MD_al:
                    flux_derived.append(
                        find_flux(Mdot=Md, nu=c / lam_passband, rin=0.0025 * R_sun, rout=10 ** 11, i=0, d=1, M=M_sun))

                flux_derived = np.array(flux_derived)
                luminosity_derived = 2 * np.pi * flux_derived
                t[f'{name_passband}'] = luminosity_derived

        except RuntimeError:
            continue
        Mdot_tables.append(t)


    direc = 'analysis_Check_Model'
    os.makedirs(direc, exist_ok=True)
    names = selected_obj
    # outbursts = data['outburst_index']
    for table_out, obj_name in zip(Mdot_tables, names):
        table_out.write(os.path.join(direc, f'{obj_name}_{outburst_idx[obj_name]}.csv'), format = 'ascii.csv', overwrite=True)






"""

######for all objects, and LSST passbands
    #find derived Luminosity in other passband
    l_tables = []
    for name, index in zip(data['name'], data['outburst_index']):
        data_Mdot = ascii.read(f'analysis_Mdot/{name}_{index}.csv')
        MD_al = data_Mdot['Mdot']
        time = data_Mdot['t']

        lambs = [3694.25e-8, 4840.83e-8, 6257.74e-8, 7560.48e-8, 8701.29e-8, 9749.32e-8]  # u, g, r, i, z, y
        passband_name = ['u', 'g', 'r', 'i', 'z', 'y']

        l_table = Table()
        for lam_passband, name_passband in zip(lambs, passband_name):
            flux_derived = []
            for Md in MD_al:
                flux_derived.append(
                    find_flux(Mdot=Md, nu=c / lam_passband, rin=0.0025 * R_sun, rout=10 ** 11, i=0, d=1, M=M_sun))

            flux_derived = np.array(flux_derived)
            luminosity_derived = 2 * np.pi * flux_derived
            m_abs = np.log10(flux_derived / 3.63e-20) / (-0.4)

            l_table[f'L_{name_passband}'] = luminosity_derived
            l_table[f'm_abs_{name_passband}'] = m_abs

        l_tables.append(l_table)

    direc = 'analysis_luminosity'
    os.makedirs(direc, exist_ok=True)
    names = data['name']
    outbursts = data['outburst_index']
    for l_table, obj_name, outburst_i in zip(l_tables, names, outbursts):
        l_table.write(os.path.join(direc, f'{obj_name}_{outburst_i}_luminosity.csv'), format = 'ascii.csv', overwrite=True)
"""


if __name__ == '__main__':
    main()