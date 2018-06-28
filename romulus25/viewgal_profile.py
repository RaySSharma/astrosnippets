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

halo = h.load_copy(halo_num)
halo.physical_units()

pynbody.analysis.angmom.faceon(h[1])
p = profile.Profile(halo.s,min='.01 kpc', max='50 kpc')

plt.plot(p['rbins'].in_units('kpc'),p['density'],'k')
plt.show()