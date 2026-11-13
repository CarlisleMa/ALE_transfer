    def addSolvent(self, forcefield, model='tip3p', boxSize=None, boxVectors=None, padding=None, numAdded=None, boxShape='cube', positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar, neutralize=True, residueTemplates=dict()):
        """Add solvent (both water and ions) to the model to fill a periodic box.

        The algorithm works as follows:

        1. Water molecules are added to fill the box.
        2. Water molecules are removed if their distance to any solute atom is less than the sum of their van der Waals radii.
        3. If the solute is charged and neutralize=True, enough positive or negative ions are added to neutralize it.  Each ion is added by
           randomly selecting a water molecule and replacing it with the ion.
        4. Ion pairs are added to give the requested total ionic strength.  Note that only monovalent ions are currently supported.

        The box size can be specified in any of several ways:

        1. You can explicitly give the vectors defining the periodic box to use.
        2. Alternatively, for a rectangular box you can simply give the dimensions of the unit cell.
        3. You can give a padding distance.  A bounding sphere containing the solute is determined, and the box size is
           set to (sphere diameter)+(padding).  This guarantees no atom in the solute will come closer than the padding
           distance to any atom of another periodic copy.  If the sphere diameter is less than the padding distance,
           the box size is set to 2*(padding) to ensure no atom is closer than the padding distance to two periodic
           copies of any other atom.
        4. You can specify the total number of molecules (both waters and ions) to add.  A box is then created whose size is
           just large enough hold the specified amount of solvent.
        5. Finally, if none of the above options is specified, the existing Topology's box vectors are used.

        When specifying either a padding distance or a number of molecules, you can specify a shape for the periodic box:
        cubic, rhombic dodecahedron, or truncated octahedron.  Using a non-rectangular box allows the same distance
        between periodic copies to be achieved with a smaller box.  The most compact option is a rhombic dodecahedron,
        for which the box volume is 70.7% the volume of a cubic box with the same amount of padding.

        There exist many different water models, many of which are very similar to each other. This method creates
        preequilibrated water boxes for a limited set of water models. In most cases, they work equally well for other models
        that involve the same number of particles. For example, to simulate a box of TIP3P-FB water, use this method to
        create a box of TIP3P water, construct a System using TIP3P-FB parameters, and perform a local energy minimization
        to correct for the small differences between the models. Likewise, a box of TIP4P-Ew water can be used for most
        four site water models.

        Parameters
        ----------
        forcefield : ForceField
            the ForceField to use for determining van der Waals radii and atomic charges
        model : str='tip3p'
            the water model to use.  Supported values are 'tip3p', 'spce', 'tip4pew', 'tip5p', and 'swm4ndp' (polarizable).
        boxSize : Vec3=None
            the size of the box to fill with water
        boxVectors : tuple of Vec3=None
            the vectors defining the periodic box to fill with water
        padding : distance=None
            the padding distance to use
        numAdded : int=None
            the total number of molecules (waters and ions) to add
        boxShape: str='cube'
            the box shape to use.  Allowed values are 'cube', 'dodecahedron', and 'octahedron'.  If padding and numAdded
            are both None, this is ignored.
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).
        """
        if len([x for x in (boxSize, boxVectors, padding, numAdded) if x is not None]) > 1:
            raise ValueError('At most one of the following arguments may be specified: boxSize, boxVectors, padding, numAdded')

        # Load the pre-equilibrated water box.

        vdwRadiusPerSigma = 0.5612310241546864907
        if model == 'tip3p':
            waterRadius = 0.31507524065751241*vdwRadiusPerSigma
        elif model == 'spce':
            waterRadius = 0.31657195050398818*vdwRadiusPerSigma
        elif model == 'tip4pew':
            waterRadius = 0.315365*vdwRadiusPerSigma
        elif model == 'tip5p':
            waterRadius = 0.312*vdwRadiusPerSigma
        elif model == 'swm4ndp':
            waterRadius = 0.318395*vdwRadiusPerSigma
        else:
            raise ValueError('Unknown water model: %s' % model)
        pdb = PDBFile(os.path.join(os.path.dirname(__file__), 'data', model+'.pdb'))
        pdbTopology = pdb.getTopology()
        pdbPositions = pdb.getPositions().value_in_unit(nanometer)
        pdbResidues = list(pdbTopology.residues())
        pdbBoxSize = pdbTopology.getUnitCellDimensions().value_in_unit(nanometer)

        # Pick a unit cell size.

        if numAdded is not None:
            # Select a padding distance which is guaranteed to give more than the specified number of molecules.

            padding = 2.2*(numAdded/((len(pdbResidues)/pdbBoxSize[0]**3)*8))**(1.0/3.0)
            if padding < 0.5:
                padding = 0.5 # Ensure we have enough when adding very small numbers of molecules
        if boxVectors is not None:
            if is_quantity(boxVectors[0]):
                boxVectors = (boxVectors[0].value_in_unit(nanometer), boxVectors[1].value_in_unit(nanometer), boxVectors[2].value_in_unit(nanometer))
            box = Vec3(boxVectors[0][0], boxVectors[1][1], boxVectors[2][2])
            vectors = boxVectors
        elif boxSize is not None:
            if is_quantity(boxSize):
                boxSize = boxSize.value_in_unit(nanometer)
            box = Vec3(boxSize[0], boxSize[1], boxSize[2])
            vectors = (Vec3(boxSize[0], 0, 0), Vec3(0, boxSize[1], 0), Vec3(0, 0, boxSize[2]))
        elif padding is not None:
            if is_quantity(padding):
                padding = padding.value_in_unit(nanometer)
            if len(self.positions) == 0:
                radius = 0
            else:
                positions = self.positions.value_in_unit(nanometer)
                minRange = Vec3(*(min((pos[i] for pos in positions)) for i in range(3)))
                maxRange = Vec3(*(max((pos[i] for pos in positions)) for i in range(3)))
                center = 0.5*(minRange+maxRange)
                radius = max(unit.norm(center-pos) for pos in positions)
            width = max(2*radius+padding, 2*padding)
            vectors = self._computeBoxVectors(width, boxShape)
            box = Vec3(vectors[0][0], vectors[1][1], vectors[2][2])
        else:
            box = self.topology.getUnitCellDimensions().value_in_unit(nanometer)
            vectors = self.topology.getPeriodicBoxVectors().value_in_unit(nanometer)
            if box is None:
                raise ValueError('Neither the box size, box vectors, nor padding was specified, and the Topology does not define unit cell dimensions')

        # Have the ForceField build a System for the solute from which we can determine van der Waals radii.

        system = forcefield.createSystem(self.topology, residueTemplates=residueTemplates)
        nonbonded = None
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), (NonbondedForce, AmoebaVdwForce)):
                nonbonded = system.getForce(i)
      
        if nonbonded is None:
            raise ValueError('The ForceField does not specify a NonbondedForce')
        cutoff = [waterRadius]*system.getNumParticles()

        for i in range(system.getNumParticles()):
            params = nonbonded.getParticleParameters(i)
            if params[2] != 0*kilojoules_per_mole:
                cutoff[i] += params[1].value_in_unit(nanometer)*vdwRadiusPerSigma

        waterCutoff = waterRadius
        if len(cutoff) == 0:
            maxCutoff = waterCutoff
        else:
            maxCutoff = max(waterCutoff, max(cutoff))

        # Copy the solute over.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(vectors*nanometer)
        newAtoms = {}
        newPositions = []*nanometer
        newResidueTemplates=dict()
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue in residueTemplates:
                    newResidueTemplates[newResidue] = residueTemplates[residue]
                for atom in residue.atoms():
                    newAtom = newTopology.addAtom(atom.name, atom.element, newResidue, atom.id, atom.formalCharge)
                    newAtoms[atom] = newAtom
                    newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        # Sort the solute atoms into cells for fast lookup.

        if len(self.positions) == 0:
            positions = []
        else:
            positions = deepcopy(self.positions.value_in_unit(nanometer))
        cells = _CellList(positions, maxCutoff, vectors, True)

        # Create a function to compute the distance between two points, taking periodic boundary conditions into account.

        periodicDistance = compiled.periodicDistance(vectors)

        # Find the list of water molecules to add.

        newChain = newTopology.addChain()
        if len(positions) == 0:
            center = Vec3(0, 0, 0)
        else:
            center = [(max((pos[i] for pos in positions))+min((pos[i] for pos in positions)))/2 for i in range(3)]
            center = Vec3(center[0], center[1], center[2])
        numBoxes = [int(ceil(box[i]/pdbBoxSize[i])) for i in range(3)]
        addedWaters = []
        for residue in pdbResidues:
            oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
            for boxx in range(numBoxes[0]):
                for boxy in range(numBoxes[1]):
                    for boxz in range(numBoxes[2]):
                        offset = Vec3(boxx*pdbBoxSize[0], boxy*pdbBoxSize[1], boxz*pdbBoxSize[2])
                        atomPos = pdbPositions[oxygen.index]+offset
                        if not any([atomPos[i] > box[i] for i in range(3)]):
                            # This molecule is inside the box, so see how close to it is to the solute.

                            atomPos += center-box/2
                            for i in cells.neighbors(atomPos):
                                if periodicDistance(atomPos, positions[i]) < cutoff[i]:
                                    break
                            else:
                                # Record this water molecule as one to add.

                                addedWaters.append((residue.index, atomPos))

        if numAdded is not None:
            # We added many more waters than we actually want.  Sort them based on distance to the nearest box edge and
            # only keep the ones in the middle.

            lowerBound = center-box/2
            upperBound = center+box/2
            distToEdge = (min(min(pos-lowerBound), min(upperBound-pos)) for index, pos in addedWaters)
            sortedIndex = [i[0] for i in sorted(enumerate(distToEdge), key=lambda x: -x[1])]
            addedWaters = [addedWaters[i] for i in sortedIndex[:numAdded]]

            # Compute a new periodic box size.

            maxSize = max(max((pos[i] for index, pos in addedWaters))-min((pos[i] for index, pos in addedWaters)) for i in range(3))
            maxSize += 0.1  # Add padding to reduce clashes at the edge.
            newTopology.setPeriodicBoxVectors(self._computeBoxVectors(maxSize, boxShape))
        else:
            # There could be clashes between water molecules at the box edges.  Find ones to remove.

            upperCutoff = center+box/2-Vec3(waterCutoff, waterCutoff, waterCutoff)
            lowerCutoff = center-box/2+Vec3(waterCutoff, waterCutoff, waterCutoff)
            lowerSkinPositions = [pos for index, pos in addedWaters if pos[0] < lowerCutoff[0] or pos[1] < lowerCutoff[1] or pos[2] < lowerCutoff[2]]
            filteredWaters = []
            cells.cells = {}
            for i in range(len(lowerSkinPositions)):
                cell = cells.cellForPosition(lowerSkinPositions[i])
                if cell in cells.cells:
                    cells.cells[cell].append(i)
                else:
                    cells.cells[cell] = [i]
            for entry in addedWaters:
                pos = entry[1]
                if pos[0] < upperCutoff[0] and pos[1] < upperCutoff[1] and pos[2] < upperCutoff[2]:
                    filteredWaters.append(entry)
                else:
                    if not any((periodicDistance(lowerSkinPositions[i], pos) < waterCutoff and norm(lowerSkinPositions[i]-pos) > waterCutoff for i in cells.neighbors(pos))):
                        filteredWaters.append(entry)
            addedWaters = filteredWaters

        # Add the water molecules.
        waterPos = {}
        for index, pos in addedWaters:
            newResidue = newTopology.addResidue(residue.name, newChain)
            residue = pdbResidues[index]
            oxygen = [atom for atom in residue.atoms() if atom.element == elem.oxygen][0]
            oPos = pdbPositions[oxygen.index]
            molAtoms = []
            for atom in residue.atoms():
                molAtoms.append(newTopology.addAtom(atom.name, atom.element, newResidue))
                newPositions.append((pos+pdbPositions[atom.index]-oPos)*nanometer)
                if atom.element == elem.oxygen:
                    waterPos[newResidue] = newPositions[-1]
            for atom1 in molAtoms:
                if atom1.element == elem.oxygen:
                    for atom2 in molAtoms:
                        if atom2.element == elem.hydrogen:
                            newTopology.addBond(atom1, atom2)

        self.topology = newTopology
        self.positions = newPositions

        # Total number of waters in the box
        numTotalWaters = len(waterPos)

        # Add ions to neutralize the system.
        self._addIons(forcefield, numTotalWaters, waterPos, positiveIon=positiveIon, negativeIon=negativeIon, ionicStrength=ionicStrength, neutralize=neutralize, residueTemplates=newResidueTemplates)
