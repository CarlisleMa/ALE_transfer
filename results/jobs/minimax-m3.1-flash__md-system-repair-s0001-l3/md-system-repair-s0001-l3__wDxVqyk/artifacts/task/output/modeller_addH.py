    def addHydrogens(self, forcefield=None, pH=7.0, variants=None, platform=None, residueTemplates=dict()):
        """Add missing hydrogens to the model.

        Some residues can exist in multiple forms depending on the pH and properties of the local environment.  These
        variants differ in the presence or absence of particular hydrogens.  In particular, the following variants
        are supported:

        Aspartic acid:
            ASH: Neutral form with a hydrogen on one of the delta oxygens
            ASP: Negatively charged form without a hydrogen on either delta oxygen

        Cysteine:
            CYS: Neutral form with a hydrogen on the sulfur
            CYX: No hydrogen on the sulfur (either negatively charged, or part of a disulfide bond)

        Glutamic acid:
            GLH: Neutral form with a hydrogen on one of the epsilon oxygens
            GLU: Negatively charged form without a hydrogen on either epsilon oxygen

        Histidine:
            HID: Neutral form with a hydrogen on the ND1 atom
            HIE: Neutral form with a hydrogen on the NE2 atom
            HIP: Positively charged form with hydrogens on both ND1 and NE2
            HIN: Negatively charged form without a hydrogen on either ND1 or NE2

        Lysine:
            LYN: Neutral form with two hydrogens on the zeta nitrogen
            LYS: Positively charged form with three hydrogens on the zeta nitrogen

        The variant to use for each residue is determined by the following rules:

        1. The most common variant at the specified pH is selected.
        2. Any Cysteine that participates in a disulfide bond uses the CYX variant regardless of pH.
        3. For a neutral Histidine residue, the HID or HIE variant is selected based on which one forms a better hydrogen bond.

        You can override these rules by explicitly specifying a variant for any residue.  To do that, provide a list for the
        'variants' parameter, and set the corresponding element to the name of the variant to use.

        A special case is when the model already contains a hydrogen that should not be present in the desired variant.
        If you explicitly specify a variant using the 'variants' parameter, the residue will be modified to match the
        desired variant, removing hydrogens if necessary.  On the other hand, for residues whose variant is selected
        automatically, this function will only add hydrogens.  It will never remove ones that are already present in the
        model, regardless of the specified pH.

        In all cases, the positions of existing atoms (including existing hydrogens) are not modified.

        Definitions for standard amino acids and nucleotides are built in.  You can call loadHydrogenDefinitions() to load
        additional definitions for other residue types.

        Parameters
        ----------
        forcefield : ForceField=None
            the ForceField to use for determining the positions of hydrogens.
            If this is None, positions will be picked which are generally
            reasonable but not optimized for any particular ForceField.
        pH : float=7.0
            the pH based on which to select variants
        variants : list=None
            an optional list of variants to use.  If this is specified, its
            length must equal the number of residues in the model.  variants[i]
            is the name of the variant to use for residue i (indexed starting at
            0). If an element is None, the standard rules will be followed to
            select a variant for that residue.  Alternatively, an element may specify
            exactly which hydrogens to add.  In that case, variants[i] should be
            a list of tuples [(name1, parent1), (name2, parent2), ...].  Each
            tuple specifies the name of a hydrogen and the name of the parent atom
            it should be bonded to.
        platform : Platform=None
            the Platform to use when computing the hydrogen atom positions.  If
            this is None, the default Platform will be used.
        residueTemplates : dict=dict()
            specifies which template the ForceField should use for particular residues.  The keys
            should be Residue objects from the Topology, and the values should be the names of the
            templates to use for them.  This is useful when a ForceField contains multiple templates
            that can match the same residue (e.g Fe2+ and Fe3+ templates in the ForceField for a
            monoatomic iron ion in the Topology).

        Returns
        -------
        list
             a list of what variant was actually selected for each residue,
             in the same format as the variants parameter
        """
        # Check the list of variants.

        residues = list(self.topology.residues())
        if variants is not None:
            if len(variants) != len(residues):
                raise ValueError("The length of the variants list must equal the number of residues")
        else:
            variants = [None]*len(residues)
        actualVariants = [None]*len(residues)

        # Load the residue specifications.

        Modeller._loadStandardHydrogenDefinitions()

        # Make a list of atoms bonded to each atom.

        bonded = {}
        for atom in self.topology.atoms():
            bonded[atom] = []
        for atom1, atom2 in self.topology.bonds():
            bonded[atom1].append(atom2)
            bonded[atom2].append(atom1)

        # Define a function that decides whether a set of atoms form a hydrogen bond, using fairly tolerant criteria.

        def isHbond(d, h, a):
            if norm(d-a) > 0.35*nanometer:
                return False
            deltaDH = h-d
            deltaHA = a-h
            deltaDH /= norm(deltaDH)
            deltaHA /= norm(deltaHA)
            return acos(dot(deltaDH, deltaHA)) < 50*degree

        # Loop over residues.

        newTopology = Topology()
        newTopology.setPeriodicBoxVectors(self.topology.getPeriodicBoxVectors())
        newAtoms = {}
        newPositions = []*nanometer
        newResidueTemplates = {}
        newIndices = []
        acceptors = [atom for atom in self.topology.atoms() if atom.element in (elem.oxygen, elem.nitrogen)]
        positions = self.positions.value_in_unit(nanometer)
        acceptorPositions = [positions[a.index] for a in acceptors]
        cells = _CellList(acceptorPositions, 0.35, None, False)
        for chain in self.topology.chains():
            newChain = newTopology.addChain(chain.id)
            for residue in chain.residues():
                newResidue = newTopology.addResidue(residue.name, newChain, residue.id, residue.insertionCode)
                if residue in residueTemplates:
                    newResidueTemplates[newResidue] = residueTemplates[residue]
                isNTerminal = (residue == chain._residues[0])
                isCTerminal = (residue == chain._residues[-1])
                if residue.name in Modeller._residueHydrogens or isinstance(variants[residue.index], list):
                    # Add hydrogens.  First select which variant to use.

                    variant = variants[residue.index]
                    if variant is None:
                        if residue.name == 'CYS':
                            # If this is part of a disulfide, use CYX.

                            sulfur = [atom for atom in residue.atoms() if atom.element == elem.sulfur]
                            if len(sulfur) == 1 and any((atom.residue != residue for atom in bonded[sulfur[0]])):
                                variant = 'CYX'
                        if residue.name == 'HIS' and pH > 6.5:
                            # See if either nitrogen already has a hydrogen attached.

                            nd1 = [atom for atom in residue.atoms() if atom.name == 'ND1']
                            ne2 = [atom for atom in residue.atoms() if atom.name == 'NE2']
                            if len(nd1) != 1 or len(ne2) != 1:
                                raise ValueError('HIS residue (%d) has the wrong set of atoms' % residue.index)
                            nd1 = nd1[0]
                            ne2 = ne2[0]
                            nd1HasHydrogen = any((atom.element == elem.hydrogen for atom in bonded[nd1]))
                            ne2HasHydrogen = any((atom.element == elem.hydrogen for atom in bonded[ne2]))
                            if nd1HasHydrogen and ne2HasHydrogen:
                                variant = 'HIP'
                            elif nd1HasHydrogen:
                                variant = 'HID'
                            elif ne2HasHydrogen:
                                variant = 'HIE'
                            else:
                                # Estimate the hydrogen positions.

                                nd1Pos = self.positions[nd1.index]
                                ne2Pos = self.positions[ne2.index]
                                hd1Delta = Vec3(0, 0, 0)*nanometer
                                for other in bonded[nd1]:
                                    hd1Delta += nd1Pos-self.positions[other.index]
                                hd1Delta *= 0.1*nanometer/norm(hd1Delta)
                                hd1Pos = nd1Pos+hd1Delta
                                he2Delta = Vec3(0, 0, 0)*nanometer
                                for other in bonded[ne2]:
                                    he2Delta += ne2Pos-self.positions[other.index]
                                he2Delta *= 0.1*nanometer/norm(he2Delta)
                                he2Pos = ne2Pos+he2Delta

                                # See whether either hydrogen would form a hydrogen bond.

                                nd1IsBonded = False
                                ne2IsBonded = False
                                for acceptorIndex in cells.neighbors(nd1Pos.value_in_unit(nanometer)):
                                    acceptor = acceptors[acceptorIndex]
                                    if acceptor.residue != residue:
                                        acceptorPos = self.positions[acceptor.index]
                                        if isHbond(nd1Pos, hd1Pos, acceptorPos):
                                            nd1IsBonded = True
                                            break
                                if not nd1IsBonded:
                                    for acceptorIndex in cells.neighbors(ne2Pos.value_in_unit(nanometer)):
                                        acceptor = acceptors[acceptorIndex]
                                        if acceptor.residue != residue:
                                            acceptorPos = self.positions[acceptor.index]
                                            if isHbond(ne2Pos, he2Pos, acceptorPos):
                                                ne2IsBonded = True
                                                break
                                if ne2IsBonded and not nd1IsBonded:
                                    variant = 'HIE'
                                else:
                                    variant = 'HID'
                        elif residue.name == 'HIS':
                            variant = 'HIP'
                    if isinstance(variant, list):
                        spec = Modeller._ResidueData(residue.name)
                        infinity = float('Inf')
                        spec.hydrogens = [Modeller._Hydrogen(name, parent, infinity, None, None) for name, parent in variant]
                    else:
                        spec = Modeller._residueHydrogens[residue.name]
                        if variant is not None and variant not in spec.variants:
                            raise ValueError('Illegal variant for %s residue: %s' % (residue.name, variant))
                    actualVariants[residue.index] = variant
                    removeExtraHydrogens = (variants[residue.index] is not None)

                    # Make a list of hydrogens that should be present in the residue.

                    parents = [atom for atom in residue.atoms() if atom.element != elem.hydrogen]
                    parentNames = [atom.name for atom in parents]
                    hydrogens = [h for h in spec.hydrogens if (variant is None and pH <= h.maxph) or (h.variants is None and pH <= h.maxph) or (h.variants is not None and variant in h.variants)]
                    hydrogens = [h for h in hydrogens if h.terminal is None or (isNTerminal and 'N' in h.terminal) or (isCTerminal and 'C' in h.terminal) or
                                 ((not isNTerminal) and (not isCTerminal) and '-' in h.terminal)]
                    hydrogens = [h for h in hydrogens if h.parent in parentNames]

                    # Loop over atoms in the residue, adding them to the new topology along with required hydrogens.

                    for parent in residue.atoms():
                        # Check whether this is a hydrogen that should be removed.

                        if removeExtraHydrogens and parent.element == elem.hydrogen and not any(parent.name == h.name for h in hydrogens):
                            continue

                        # Add the atom.

                        newAtom = newTopology.addAtom(parent.name, parent.element, newResidue)
                        newAtoms[parent] = newAtom
                        newPositions.append(deepcopy(self.positions[parent.index]))
                        if parent in parents:
                            # Match expected hydrogens with existing ones and find which ones need to be added.

                            existing = [atom for atom in bonded[parent] if atom.element == elem.hydrogen]
                            expected = [h for h in hydrogens if h.parent == parent.name]
                            if len(existing) < len(expected):
                                # Try to match up existing hydrogens to expected ones.

                                matches = []
                                for e in existing:
                                    match = [h for h in expected if h.name == e.name]
                                    if len(match) > 0:
                                        matches.append(match[0])
                                        expected.remove(match[0])
                                    else:
                                        matches.append(None)

                                # If any hydrogens couldn't be matched by name, just match them arbitrarily.

                                for i in range(len(matches)):
                                    if matches[i] is None:
                                        matches[i] = expected[-1]
                                        expected.remove(expected[-1])

                                # Add the missing hydrogens.

                                for h in expected:
                                    newH = newTopology.addAtom(h.name, elem.hydrogen, newResidue)
                                    newIndices.append(newH.index)
                                    delta = Vec3(0, 0, 0)*nanometer
                                    if len(bonded[parent]) > 0:
                                        for other in bonded[parent]:
                                            delta += self.positions[parent.index]-self.positions[other.index]
                                    else:
                                        delta = Vec3(random.random(), random.random(), random.random())*nanometer
                                    delta *= 0.1*nanometer/norm(delta)
                                    delta += 0.05*Vec3(random.random(), random.random(), random.random())*nanometer
                                    delta *= 0.1*nanometer/norm(delta)
                                    newPositions.append(self.positions[parent.index]+delta)
                                    newTopology.addBond(newAtom, newH)
                else:
                    # Just copy over the residue.

                    for atom in residue.atoms():
                        newAtom = newTopology.addAtom(atom.name, atom.element, newResidue)
                        newAtoms[atom] = newAtom
                        newPositions.append(deepcopy(self.positions[atom.index]))
        for bond in self.topology.bonds():
            if bond[0] in newAtoms and bond[1] in newAtoms:
                newTopology.addBond(newAtoms[bond[0]], newAtoms[bond[1]], bond.type, bond.order)

        # The hydrogens were added at random positions.  Now perform an energy minimization to fix them up.

        addedH = set(newIndices)  # keep track of Hs added

        if forcefield is not None:
            # Use the ForceField the user specified.

            system = forcefield.createSystem(newTopology, rigidWater=False, nonbondedMethod=CutoffNonPeriodic, residueTemplates=newResidueTemplates)
            for i in range(system.getNumParticles()):
                if i not in addedH:
                    # Existing atom, make it immobile.
                    system.setParticleMass(i, 0)
        else:
            # Create a System that restrains the distance of each hydrogen from its parent atom
            # and causes hydrogens to spread out evenly.

            system = System()
            nonbonded = CustomNonbondedForce('100/(r/0.1)^4')
            nonbonded.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic);
            nonbonded.setCutoffDistance(1*nanometer)
            bonds = HarmonicBondForce()
            angles = HarmonicAngleForce()
            system.addForce(nonbonded)
            system.addForce(bonds)
            system.addForce(angles)
            bondedTo = []
            for atom in newTopology.atoms():
                nonbonded.addParticle([])
                if atom.index not in addedH:  # make immobile
                    system.addParticle(0.0)
                else:
                    system.addParticle(1.0)
                bondedTo.append([])
            for atom1, atom2 in newTopology.bonds():
                if atom1.element == elem.hydrogen or atom2.element == elem.hydrogen:
                    bonds.addBond(atom1.index, atom2.index, 0.1, 100000.0)
                bondedTo[atom1.index].append(atom2)
                bondedTo[atom2.index].append(atom1)
            for residue in newTopology.residues():
                if residue.name == 'HOH':
                    # Add an angle term to make the water geometry correct.

                    atoms = list(residue.atoms())
                    oindex = [i for i in range(len(atoms)) if atoms[i].element == elem.oxygen]
                    if len(atoms) == 3 and len(oindex) == 1:
                        hindex = list(set([0,1,2])-set(oindex))
                        angles.addAngle(atoms[hindex[0]].index, atoms[oindex[0]].index, atoms[hindex[1]].index, 1.824, 836.8)
                else:
                    # Add angle terms for any hydroxyls.

                    for atom in residue.atoms():
                        index = atom.index
                        if atom.element == elem.oxygen and len(bondedTo[index]) == 2 and elem.hydrogen in (a.element for a in bondedTo[index]):
                            angles.addAngle(bondedTo[index][0].index, index, bondedTo[index][1].index, 1.894, 460.24)

        if platform is None:
            context = Context(system, VerletIntegrator(0.0))
        else:
            context = Context(system, VerletIntegrator(0.0), platform)
        context.setPositions(newPositions)
        LocalEnergyMinimizer.minimize(context, 1.0, 50)
        self.topology = newTopology
        self.positions = context.getState(positions=True).getPositions()
        del context
        return actualVariants
