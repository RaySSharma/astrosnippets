#!/usr/bin/env python

import matplotlib
matplotlib.use('tkAgg')

import matplotlib.pyplot as plt
import pynbody

# What do:
# 1 = Zoom out
# 2 = Zoom in
# 0 = Exit

res = 0.175

fname = '/data/REPOSITORY/cosmo/romulus25/cosmo25p.768sg1bwK1BHe75.'

timestep = str(input('Timestep: ')).zfill(6)
halo_num = int(input('Halo: '))

sim = pynbody.load(fname+timestep)

h = sim.halos(dosort=True)

rmax = h[halo_num].properties['Rvir']
wid = .25*rmax

print('N_gas = ', h[halo_num].properties['n_gas'])
print('N_star = ', h[halo_num].properties['n_star'])
print('M_g/M_b = ', h[halo_num].properties['M_gas']/(h[halo_num].properties['M_star']+h[halo_num].properties['M_gas']))

halo = h.load_copy(halo_num)
halo.physical_units()


f, ax = plt.subplots(1,2)
plt.ion()
whatdo = 0

while whatdo == 0:
    pynbody.analysis.angmom.faceon(halo)
    pynbody.plot.stars.render(halo.s, axes=ax[0], starsize=res, clear=False, width=wid)

    pynbody.analysis.angmom.sideon(halo)
    pynbody.plot.stars.render(halo.s, axes=ax[1], starsize=res, clear=False, width=wid)
    plt.show()

    z = int(input('What should I do (1: Zoom-out; 2: Zoom-in)? '))
    if z == 1:
        wid *= 1.5
    elif z == 2:
        wid *= 0.5
    else:
        whatdo = 1

'''
Created on Sep 6, 2017

@author: anna
'''
