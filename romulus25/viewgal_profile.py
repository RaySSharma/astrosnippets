#!/usr/bin/env python
import matplotlib
matplotlib.use('tkAgg')

import matplotlib.pyplot as plt
import pynbody

'''
	Script to display a density profile of the stellar component of the halo of interest.
	Vars:
		res: Resolution of stars in mock image, see pynbody.plot.stars.render() for more info.
		fname: Location of simulation and name of simulation up to the timestep of interest
			ex: ~/romulus25/cosmo25p.768sg1bwK1BHe75.
		timestep: Timestep of interest in simulation
		halo_num: Halo number compatible with pynbody halo indexing

Created 06/28/18
- Ray Sharma
'''

### Parameters to change ###
nrow = 5
ncol = 4 # Select an even number of columns
res = 0.175
fname = '/data/REPOSITORY/cosmo/romulus25/cosmo25p.768sg1bwK1BHe75.'
timestep = 7779
halo_num = 323
outfile = 'profile.png'
############################

sim = pynbody.load(fname+timestep)
h = sim.halos(dosort=True)

halo = h.load_copy(halo_num)
halo.physical_units()

pynbody.analysis.angmom.faceon(h[1])
p = profile.Profile(halo.s,min='.01 kpc', max='50 kpc')

plt.plot(p['rbins'].in_units('kpc'),p['density'],'k')
plt.savefig(outfile)