import tangos as db
import scipy.spatial as spatial
import pynbody
from numpy import log10

'''
	Script to interactively find the nearest neighbors to a halo within some radius.
	Vars:
		fname: Location of simulation and name of simulation up to the timestep of interest
			ex: ~/romulus25/cosmo25p.768sg1bwK1BHe75.
		db_name: Tangos database name
			ex: cosmo25

Created on 06/29/18
- Ray Sharma
'''

### Parameters to change ###
fname = '/data/REPOSITORY/cosmo/romulus25/cosmo25p.768sg1bwK1BHe75.'
db_name = 'cosmo25'
############################

timestep = input('Timestep: ')
timestep = str(timestep).zfill(6)
sim = pynbody.load(fname+timestep)
ts = db.get_timestep(db_name+'/%'+timestep)

SSC_all = ts.calculate_all('SSC')[0]
tree = spatial.cKDTree(SSC_all)

while True:
	try: halo_num = int(input('Halo number: '))
	except: break
	Rvir = ts[halo_num]['Rvir']
	Mstar = log10(ts[halo_num]['Mstar'])
	SSC = ts[halo_num]['SSC']
	print('Virial radius: ', Rvir, ' kpc')
	print('Stellar mass: ', Mstar, 'log sol mass')

	while True:
		try: search_radius = float(input('Search radius: '))
		except: break 

		query = tree.query_ball_point(SSC, search_radius)
		print('Nearby objects: ', len(query))

		masses = [log10(ts[h]['Mstar']) for h in query]
		print('Stellar masses: ', masses, 'log sol mass')