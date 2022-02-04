#!/usr/bin/env python3

import pandas as pd
from scipy import optimize
from functools import partial
import matplotlib.pyplot as plt
import numpy as np

def loaddata(name, oid, filter, mjd0, mjd1 ):
    def str_to_np(s, dtype=np.float64):
        return np.fromstring(s[1:-1], sep=',', dtype=dtype)
    df = pd.read_csv(name, converters={'filter': partial(str_to_np, dtype=np.uint8), 'mjd': str_to_np, 'mag': str_to_np, 'magerr': str_to_np}, index_col='vsx_oid')

    row = df.loc[oid]
    idx = (row['filter'] == filter)  & (row['mjd'] >= mjd0) & (row['mjd'] <= mjd1)
    t = row['mjd'][idx]
    m = row['mag'][idx]
    sig = row['magerr'][idx]
    return(t, m, sig)

def piecewise(t,t0,dt0,m0,dm0,k0,k1):
    return np.piecewise(t , [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)] ,
                        [lambda t:k0*(t-t0) + m0,
                         lambda t:(t-t0)*(dm0)/dt0+m0,
                        lambda t:k1*(t-t0-dt0) + m0+dm0])
def jaco(t,t0,dt0,m0,dm0,k0,k1):
    dfdt0 = np.piecewise(t , [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)] ,
                          [lambda t: -k0,
                           lambda t:(-(dm0)*(dt0)+(t-t0)*(dm0))/((dt0)**2),
                           lambda t: -k1])
    dfdt1 = np.piecewise(t, [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)],
                          [lambda t: 0,
                           lambda t:-(t-t0)*(dm0)/((dt0)**2),
                           lambda t: 0])
    dfdm0 = np.piecewise(t, [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)],
                          [lambda t: 1,
                           lambda t:-(t-t0)/(dt0)+1,
                           lambda t: 0])
    dfdm1 = np.piecewise(t, [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)],
                          [lambda t: 0,
                           lambda t:(t-t0)/(dt0),
                           lambda t: 1])
    dfdk0 = np.piecewise(t, [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)],
                          [lambda t: (t - t0),
                           lambda t: 0,
                           lambda t: 0])
    dfdk1 = np.piecewise(t, [t <= t0, np.logical_and(t0<t, t<= (t0+dt0)),t>(t0+dt0)],
                          [lambda t: 0,
                           lambda t: 0,
                           lambda t: (t - (t0+dt0))])
    array = np.stack((dfdt0, dfdt1, dfdm0, dfdm1, dfdk0, dfdk1), axis= 0)
    return np.transpose(array)

def fit_piecewise(t, m):
    t0_l = 58322
    t0_u = 58323
    dt0_l = 0
    dt0_u = np.inf
    m0_l = 0
    m0_u = 14.9
    dm0_l = 0
    dm0_u = np.inf
    k0_l = -np.inf
    k0_u = 0
    k1_l = 0
    k1_u = np.inf
    bounds_lower = [t0_l, dt0_l, m0_l, dm0_l, k0_l, k1_l]
    bounds_upper = [t0_u, dt0_u, m0_u, dm0_u, k0_u, k1_u]
    p, e = optimize.curve_fit(piecewise, t, m, jac=jaco,
                              bounds=(bounds_lower , bounds_upper))
    return p

def plot_piecewise_fit(t, m, p):
    fig, axes = plt.subplots()
    plt.plot(t, piecewise(t, *p))
    plt.scatter(t, m, s=30, c='b')
    axes.invert_yaxis()

    plt.show()

def main():
    oid = 13913
    name = 'Dwarf nova.csv'
    filter = 2
    mjd0 = (58000 + 318.5)
    mjd1 = (58000+333)
    t, m, sig = loaddata(name, oid, filter, mjd0, mjd1 )

    p = fit_piecewise(t, m)
    plot_piecewise_fit(t, m, p)


if __name__ == '__main__':
    main()