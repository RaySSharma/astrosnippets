#!/usr/bin/env python
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import pynbody

'''
	Script to interactively display a faceon and sideon view of the stellar component of a halo, and allow for changing of field of view on the fly
	Vars:
		res: Resolution of stars in mock image, see pynbody.plot.stars.render() for more info.
		fname: Location of simulation and name of simulation up to the timestep of interest
			ex: ~/romulus25/cosmo25p.768sg1bwK1BHe75.
		timestep: Timestep of interest in simulation
		halo_num: Halo number compatible with pynbody halo indexing
	Usage:
		1 = Zoom out
		2 = Zoom in
		Anything else = Exit

Created on Sep 6, 2017
@author: Anna Wright
@author: Ray Sharma
'''

### Parameters to change ###
nrow = 5
ncol = 4 # Select an even number of columns
res = 0.175
fname = '/data/REPOSITORY/cosmo/romulus25/cosmo25p.768sg1bwK1BHe75.'
timestep = 7779
halo_num = 323
############################


timestep = str(timestep).zfill(6)
sim = pynbody.load(fname+timestep)
h = sim.halos(dosort=True)

rmax = h[halo_num].properties['Rvir']
wid = .25*rmax

print('Properties of halo ', halo_num)
print('N_gas = ', h[halo_num].properties['n_gas'])
print('N_star = ', h[halo_num].properties['n_star'])
print('M_g/M_b = ', h[halo_num].properties['M_gas']/(h[halo_num].properties['M_star']+h[halo_num].properties['M_gas']))
print('R_vir = ', rmax)

halo = h.load_copy(halo_num)
halo.physical_units()

f, ax = plt.subplots(1,2)
plt.ion()

while True:
    pynbody.analysis.angmom.faceon(halo)
    pynbody.plot.stars.render(halo.s, axes=ax[0], starsize=res, clear=False, width=wid)

    pynbody.analysis.angmom.sideon(halo)
    pynbody.plot.stars.render(halo.s, axes=ax[1], starsize=res, clear=False, width=wid)
    plt.show()

    z = int(input('What should I do (1: Zoom-out; 2: Zoom-in)? '))
    if z == 1: wid *= 1.5
    elif z == 2: wid *= 0.5
    else: break

'''
Created on Sep 6, 2017

@author: anna
'''
