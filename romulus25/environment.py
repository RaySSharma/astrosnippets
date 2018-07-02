'''
Created 07/02/18
- Ray Sharma
'''
def arrFilter(arr, condition):
    '''
        Filter numpy array in place by condition
        Usage: arrFilter(arr, condition)
    '''
    arr = arr[condition]
    return

def wrap(relpos, boxsize=25e3):
    '''
        Wrap coordinates relative to some point assuming periodic boundary conditions on all sides of the box and comoving coordinates. Scale by a(z) to use proper coordinates.
        Usage: wrap(relpos)
    @author: Michael Tremmel
    '''
    bphys = boxsize
    bad = np.where(np.abs(relpos) > bphys/2.)
    if type(bphys) == np.ndarray:
        relpos[bad] = -1.0 * (relpos[bad] / np.abs(relpos[bad])) * np.abs(bphys[bad] - np.abs(relpos[bad]))
    else:
        relpos[bad] = -1.0 * (relpos[bad] / np.abs(relpos[bad])) * np.abs(bphys - np.abs(relpos[bad]))
    return

def findCentrals(halo_sample, SSC_sample, Rvir_sample, Mstar_sample, halo_all, SSC_all, Mstar_all, boxsize=25e3):
    '''
    Identifies central halos in the current timestep using a scipy cKDTree and checking for halos within one virial radius of each halo.
    Args: 
        halo_sample: halo numbers of targets
        SSC_sample: [x,y,z] coordinates of targets
        Rvir_sample: virial radii of targets
        Mstar_sample: stellar mass of targets
        halo_all: halo numbers of all objects in timestep
        SSC_all: [x,y,z] coordinates of all objects in timestep
        Mstar_all: stellar mass of targets
        boxsize: size of simulation box in comoving kpc
    Output:
        central_flag: array to flag central objects with a 1.0
    '''
    wrap(SSC_all); wrap(SSC_sample)
    SSC_all += (boxsize / 2.) #Shift coordinates from [boxsize/2, boxsize/2] to [0, boxsize]
    SSC_sample += (boxsize / 2.)

    tree = cKDTree(SSC_all, compact_nodes=True, boxsize=boxsize*np.ones(len(SSC_all.T))) #Create periodic KDTree over all data, wrapped over [0, boxsize)
    central_flag = np.zeros(len(halo_sample)) #Array to flag central objects
    matches = [] #List to accumulate central objects
    
    for i in range(len(halo_sample)):
        query = np.array(tree.query_ball_point(SSC_sample[i], Rvir_sample[i])) #Query KDTree for objects within Rvir_sample[i] of ith object
        query.remove(arrFilter(query, halo_all[query] == halo_sample[i])) #Remove self from query

        if len(query) == 0:
            matches += [i]
        else:
            matches += [Mvir_all[query].argmax()]

    try: matches = np.unique(matches)
    except: matches = np.zeros(len(halo_sample))
    central_flag[matches] = 1.0

    return central_flag

def findIsolated(halo_sample, SSC_sample, Rvir_sample, Mstar_sample, halo_all, SSC_all, Mstar_all, boxsize=25e3):
    '''
    Identifies isolated halos in the given timestep with a scipy cKDTree using results from Geha (2012)
    Args: 
        halo_sample: halo numbers of targets
        SSC_sample: [x,y,z] coordinates of targets
        Rvir_sample: virial radii of targets
        Mstar_sample: stellar mass of targets
        halo_all: halo numbers of all objects in timestep
        SSC_all: [x,y,z] coordinates of all objects in timestep
        Mstar_all: stellar mass of targets
        boxsize: size of simulation box in comoving kpc
    Output:
        isolated_flag: array to flag isolated objects with a 1.0
    '''
    wrap(SSC_all); wrap(SSC_sample)
    SSC_all += (boxsize / 2.) #Shift coordinates from [boxsize/2, boxsize/2] to [0, boxsize]
    SSC_sample += (boxsize / 2.)

    tree = cKDTree(SSC_all, compact_nodes=True, boxsize=boxsize*np.ones(len(SSC_all.T))) #Create periodic KDTree over all data, wrapped over [0, boxsize)
    isolated_flag = np.zeros(len(halo_sample)) #Array to flag isolated objects
    matches = [] #List to accumulate isolated objects

    for i in range(len(halo_sample)):
        query = np.array(tree.query_ball_point(SSC_sample[i], Rvir_sample[i])) #Query KDTree for objects within Rvir_sample[i] of ith object
        query.remove(arrFilter(query, halo_all[query] == halo_sample[i])) #Remove self from query

        if len(query) == 0 & halo.Mstar > 10: #Check if object's stellar mass above 10^10 sol Mass
            matches += [i]
        elif len(query) == 0 & halo.Mstar <= 10: #Check if object's stellar mass below 10^10 sol Mass
            query = np.array(tree.query_ball_point(SSC_sample[i], 1500.0)) #Query KDTree for objects within 1.5 Mpc of ith object
            query.remove(arrFilter(query, halo_all[query] == halo_sample[i]))

            if all(Mstar[query] < np.log10(2.5e10)): #Check if any nearby objects bigger than 2.5*10^10 sol Mass
                matches += [i]

    try: matches = np.unique(matches)
    except: matches = np.zeros(len(halo_sample))
    isolated_flag[matches] = 1.0

    return isolated_flag