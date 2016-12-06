def findForce(system, forcetype, add=True):
    """ Finds a specific force in the system force list - added if not found."""
    for force in system.getForces():
        if isinstance(force, forcetype):
            return force
    if add==True:
        system.addForce(forcetype())
        return findForce(system, forcetype)
    return None

def setGlobalForceParameter(force, key, value):
    for i in range(force.getNumGlobalParameters()):
        if force.getGlobalParameterName(i)==key:
            print('setting force parameter', key, '=', value)
            force.setGlobalParameterDefaultValue(i, value);    

def atomIndexInResidue(residue):
    """ list of atom index in residue """
    index=[]
    for a in list(residue.atoms()):
        index.append(a.index)
    return index

def getResiduePositions(residue, positions):
    """ Returns array w. atomic positions of residue """
    ndx = atomIndexInResidue(residue)
    return np.array(positions)[ndx]

def uniquePairs(index):
    """ list of unique, internal pairs in list
    
    Parameters
    ----------
    index : list
        List of index
        
    Example
    -------
    >>> print( uniquePairs([0,1,2]) )
    [[0,1],[[1,2]]
    """
    return list(combinations( range(index[0],index[-1]+1),2 ) )

def addHarmonicConstraint(harmonicforce, pairlist, positions, threshold, k, verbose=False):
    """Add harmonic bonds between pairs if distance is smaller than threshold
    
    Parameters
    ----------
    harmoniceforce : OpenMM HarmonicForce object
        Bonds are added here
    pairlist : numpy array
        List of atom pair index
    positions : numpy array
        All positions in the system
    threshold : float
        Distance cutoff for adding constraint
    k : float
        Force constant of constraint
    verbose : bool
        Verbose output
    """
    if verbose:
        print('Constraint force constant =', k)
    for i,j in pairlist:
        distance = unit.norm( positions[i]-positions[j] )
        if distance<threshold:
            harmonicforce.addBond( i,j,
                                   distance.value_in_unit(unit.nanometer),
                                   k.value_in_unit( unit.kilojoule/unit.nanometer**2/unit.mole ))
            if verbose:
                print("added harmonic bond between", i, j, 'with distance',distance)

def addExclusions(nonbondedforce, pairlist):
    """ add nonbonded exclusions between pairs """
    for i,j in pairlist:
        nonbondedforce.addExclusion(i,j)

def rigidifyResidue(residue, harmonicforce, positions, nonbondedforce=None,
                    threshold=6.0*unit.angstrom, k=2500*unit.kilojoule/unit.nanometer**2/unit.mole):
    """Make residue rigid by adding constraints and nonbonded exclusions
    
    Parameters
    ----------
    residue : OpenMM residue
        Residue to rigidify
    harmonicforce : HarmonicBondForce
        Constraints are added to this
    positions : np.array
        All positions in the system
    nonbondedforce : NonbondedForce
        Exclusions are added to this
    """
    
    index    = atomIndexInResidue(residue)
    pairlist = uniquePairs(index)
    addHarmonicConstraint(harmonic, pairlist, pdb.positions, threshold, k)
    if nonbondedforce is not None:
        for i,j in pairlist:
            print('added nonbonded exclusion between', i, j)
            nonbonded.addExclusion(i,j)
            
def centerOfMass(index, pos, box):
    """Calculates the geometric center taking into account periodic boundaries
    
    Parameters
    ----------
    index : list
       Index of the particles to caculate COM for
    pos : numpy array
        All positions in the system
    box : numpy array
        Box side lengths for periodic boundaries
        
    Returns
    -------
    numpy array
        Mass center vector
    
    Notes
    -----
        More here: https://en.wikipedia.org/wiki/Center_of_mass#Systems_with_periodic_boundary_conditions
        Numpy style docstring: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html
    """
    theta=np.divide(positions[index], box).astype(np.float) * 2*np.pi
    x1=np.array( [np.cos(theta[:,0]).mean(), np.cos(theta[:,1]).mean(), np.cos(theta[:,2]).mean()] )
    x2=np.array( [np.sin(theta[:,0]).mean(), np.sin(theta[:,1]).mean(), np.sin(theta[:,2]).mean()] )
    return box * (np.arctan2(-x1,-x2)+np.pi) / (2*np.pi)
