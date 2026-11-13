@needs_openmm
def load_topology(topology, system=None, xyz=None, box=None, condense_atom_types=True):
    """
    Creates a :class:`parmed.structure.Structure` instance from an OpenMM
    Topology, optionally filling in parameters from a System

    Parameters
    ----------
    topology : :class:`openmm.app.Topology`
        The Topology instance with the list of atoms and bonds for this system
    system : :class:`openmm.System` or str, optional
        If provided, parameters from this System will be applied to the
        Structure. If a string is given, it will be interpreted as the file name
        of an XML-serialized System, and it will be deserialized into a System
        before used to supply parameters
    xyz : str or array of float
        Name of a file containing coordinate information or an array of
        coordinates. If file has unit cell information, it also uses that
        information unless ``box`` (below) is also specified
    box : array of 6 floats
        Unit cell dimensions
    condense_atom_types : bool, default=True
        If True, create unique atom types based on de-duplicating properties. If False,
        create one atom type for each atom in the system, even if its properties match
        an existing atom type.

    Returns
    -------
    struct : :class:`Structure <parmed.structure.Structure>`
        The structure from the provided topology

    Raises
    ------
    OpenMMWarning if parameters are found that cannot be interpreted or
    processed by ParmEd

    TypeError if there are any mismatches between the provided topology and
    system (e.g., they have different numbers of atoms)

    IOError if system is a string and it is not an existing file

    Notes
    -----
    Due to its flexibility with CustomForces, it is entirely possible that the
    functional form of the potential will be unknown to ParmEd. This function
    will try to use the energy expression to identify supported potential types
    that are implemented as CustomForce objects. In particular, quadratic
    improper torsions, when recognized, will be extracted.

    Other CustomForces, including the CustomNonbondedForce used to implement
    NBFIX (off-diagonal L-J modifications) and the 12-6-4 potential, will not be
    processed and will result in an unknown functional form.

    If an OpenMM Atom.id attribute is populated by a non-integer, it will be
    used to name the corresponding ParmEd AtomType object.
    """
    import openmm as mm
    from openmm import app

    struct = Structure()
    atommap = dict()
    for c in topology.chains():
        chain = c.id
        for r in c.residues():
            resname = r.name
            resnum = int(r.id) if r.id.lstrip('-').isdigit() else r.index
            inscode = r.insertionCode
            for a in r.atoms():
                if a.element is None:
                    atom = ExtraPoint(name=a.name)
                else:
                    try:
                        aid = int(a.id)
                    except ValueError:
                        aid = a.id
                    atype = aid if not isinstance(aid, int) else ''
                    atom = Atom(atomic_number=a.element.atomic_number,
                                name=a.name, mass=a.element.mass, type=atype)
                struct.add_atom(atom, resname, resnum, chain, inscode)
                atommap[a] = atom

    bond_orders_by_openmm_type = {
        app.Single: 1.0,
        app.Double: 2.0,
        app.Triple: 3.0,
        app.Amide: 1.25,
        app.Aromatic: 1.5,
    }
    for b in topology.bonds():
        if b.type in bond_orders_by_openmm_type:
            order = bond_orders_by_openmm_type[b.type]
        elif b.order is not None:
            order = float(b.order)
        else:
            order = 1.0
        struct.bonds.append(Bond(atommap[b.atom1], atommap[b.atom2], order=order))

    vectors = topology.getPeriodicBoxVectors()
    if vectors is not None:
        leng, ang = box_vectors_to_lengths_and_angles(*vectors)
        leng = leng.value_in_unit(u.angstroms)
        ang = ang.value_in_unit(u.degrees)
        struct.box = [leng[0], leng[1], leng[2], ang[0], ang[1], ang[2]]

    loaded_box = False

    if xyz is not None:
        if isinstance(xyz, str):
            xyz = load_file(xyz, skip_bonds=True)
            struct.coordinates = xyz.coordinates
            if struct.box is not None:
                if xyz.box is not None:
                    loaded_box = True
                    struct.box = xyz.box
        else:
            struct.coordinates = xyz

    if box is not None:
        loaded_box = True
        struct.box = box

    if struct.box is not None:
        struct.box = np.asarray(struct.box)

    if system is None:
        return struct

    if isinstance(system, str):
        system = load_file(system)

    if not isinstance(system, mm.System):
        raise TypeError('system must be an OpenMM System object or serialized XML of an OpenMM System object')

    # We have a system, try to extract parameters from it
    if len(struct.atoms) != system.getNumParticles():
        raise TypeError(f'Topology and System have different numbers of atoms ({len(struct.atoms)} vs. {system.getNumParticles()})')

    processed_forces = set()
    ignored_forces = (mm.CMMotionRemover, mm.AndersenThermostat, mm.MonteCarloBarostat,
                      mm.MonteCarloAnisotropicBarostat, mm.MonteCarloMembraneBarostat,
                      mm.CustomExternalForce, mm.GBSAOBCForce, mm.CustomGBForce)

    if system.usesPeriodicBoundaryConditions():
        if not loaded_box:
            vectors = system.getDefaultPeriodicBoxVectors()
            leng, ang = box_vectors_to_lengths_and_angles(*vectors)
            leng = leng.value_in_unit(u.angstroms)
            ang = ang.value_in_unit(u.degrees)
            struct.box = np.asarray([leng[0], leng[1], leng[2], ang[0], ang[1], ang[2]])
    else:
        struct.box = None

    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            if mm.HarmonicBondForce in processed_forces:
                # Try to process this HarmonicBondForce as a Urey-Bradley term
                _process_urey_bradley(struct, force)
            else:
                _process_bond(struct, force)
        elif isinstance(force, mm.HarmonicAngleForce):
            _process_angle(struct, force)
        elif isinstance(force, mm.PeriodicTorsionForce):
            _process_dihedral(struct, force)
        elif isinstance(force, mm.RBTorsionForce):
            _process_rbtorsion(struct, force)
        elif isinstance(force, mm.CustomTorsionForce):
            if not _process_improper(struct, force):
                struct.unknown_functional = True
                warnings.warn('Unknown functional form of CustomTorsionForce', OpenMMWarning)
        elif isinstance(force, mm.CMAPTorsionForce):
            _process_cmap(struct, force)
        elif isinstance(force, mm.NonbondedForce):
            _process_nonbonded(struct, force, condense_atom_types)
        elif isinstance(force, ignored_forces):
            continue
        else:
            struct.unknown_functional = True
            warnings.warn(f'Unsupported Force type {force.__class__.__name__}', OpenMMWarning)
        processed_forces.add(type(force))

    return struct
