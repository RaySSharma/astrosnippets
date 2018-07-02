#!/usr/bin/env python
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pynbody
import scalebar
'''
	Script to display a grid of halo and sideon mock images of galaxies. Two images are produced for each halo.
	Vars:
		nrow: Number of rows of the grid of images
		ncol: Number of columns of the grid of images
			nrow x ncol should be equal to twice the number of halos. ncol should be even for correct formatting.
		res: Resolution of stars in mock image, see pynbody.plot.stars.render() for more info.
		fname: Location of simulation and name of simulation up to the timestep of interest
			ex: ~/romulus25/cosmo25p.768sg1bwK1BHe75.
		timestep: Timestep of interest in simulation
		halos: List of halo numbers compatible with pynbody halo indexing
		outfile: Output image file
		del_panels: Number of panels to delete starting from bottom right
Created 06/28/18
- Ray Sharma
'''
### Parameters to change ###
nrow = 6
ncol = 6 # Select an even number of columns
res = 0.175
fname = '/data/REPOSITORY/cosmo/romulus25/cosmo25p.768sg1bwK1BHe75.'
timestep = 7779
halos = [99, 156, 163, 168, 182, 187, 208, 210, 218, 222, 232, 242, 248,
       259, 323, 333]
outfile = 'dwarfs.png'
del_panels = 4
############################


timestep = str(timestep).zfill(6)
sim = pynbody.load(fname+timestep)
h = sim.halos(dosort=True)

f = plt.figure(figsize=(ncol+1, nrow+1)) 

gs = gridspec.GridSpec(nrow, ncol,
         wspace=0.0, hspace=0.0, 
         top=1.-0.5/(nrow+1), bottom=0.5/(nrow+1), 
         left=0.5/(ncol+1), right=1-0.5/(ncol+1)) 

halo_count = 0

for i in range(nrow):
	j = 0
	while True:
		if j > ncol - 1:
			break

		print('Processing halo: '+str(halo_count))
		try: halo = h.load_copy(halos[halo_count])
		except: break
		halo.physical_units()
		pynbody.analysis.angmom.sideon(halo)

		ax1 = plt.subplot(gs[i,j])
		ax2 = plt.subplot(gs[i,j+1])

		pynbody.plot.sph.image(halo.d, cmap=plt.get_cmap('hot'), vmin=1e5, vmax=1e9, width=100, subplot=ax1, show_cbar=False)
		pynbody.plot.stars.render(halo.s, axes=ax2, starsize=res, clear=False, width=20)

		scalebar.add_scalebar(ax=ax1, matchx=True, matchy=False, hidex=True, hidey=True, sizey=0, labelx=None)
		scalebar.add_scalebar(ax=ax2, matchx=True, matchy=False, hidex=True, hidey=True, sizey=0, labelx=None)

		ax1.text(0.08,0.8, 'h'+str(halos[halo_count]), transform=ax1.transAxes, fontsize=14, color='w')
		ax1.text(0.55, 0.1, '50 kpc', transform=ax1.transAxes, fontsize=9, color='white')
		ax2.text(0.08,0.8, 'h'+str(halos[halo_count]), transform=ax2.transAxes, fontsize=14, color='w')
		ax2.text(0.6, 0.1, '10 kpc', transform=ax2.transAxes, fontsize=9, color='white')

		ax1.set_xticklabels([]); ax1.set_yticklabels([])
		ax2.set_xticklabels([]); ax2.set_yticklabels([])

		j += 2
		halo_count += 1

for i in range(1, del_panels+1):
	delaxis = plt.subplot(gs[nrow-1,ncol-i]); f.delaxes(delaxis)
f.savefig(outfile)
