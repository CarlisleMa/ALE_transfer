    def addMissingAtoms(self, seed=None):
        """Add all missing heavy atoms, as specified by the missingAtoms, missingTerminals, and missingResidues fields.

        Parameters
        ----------
        seed : int
            Integer to set the random seed number of the integrator used in the minimization of the
            coordinates of the newly-added atoms.

        Notes
        -----
        You must already have called findMissingAtoms() to have identified atoms to be added.

        Examples
        --------

        Find missing heavy atoms in Abl kinase structure.

        >>> fixer = PDBFixer(pdbid='2F4J')
        >>> fixer.findMissingResidues()
        >>> fixer.findMissingAtoms()
        >>> fixer.addMissingAtoms()

        """

        # Create a Topology that 1) adds missing atoms, 2) removes all hydrogens, and 3) removes unknown molecules.

        (newTopology, newPositions, newAtoms, existingAtomMap) = self._addAtomsToTopology(True, True)
        if len(newAtoms) == 0:

            # No atoms were added, but new bonds might have been created.

            newBonds = set(newTopology.bonds())
            for atom1, atom2 in self.topology.bonds():
                if atom1 in existingAtomMap and atom2 in existingAtomMap:
                    a1 = existingAtomMap[atom1]
                    a2 = existingAtomMap[atom2]
                    if (a1, a2) in newBonds:
                        newBonds.remove((a1, a2))
                    elif (a2, a1) in newBonds:
                        newBonds.remove((a2, a1))

            # Add the new bonds to the original Topology.

            inverseAtomMap = dict((y,x) for (x,y) in existingAtomMap.items())
            for atom1, atom2 in newBonds:
                self.topology.addBond(inverseAtomMap[atom1], inverseAtomMap[atom2])
        else:

            # Create a System for energy minimizing it.

            forcefield = self._createForceField(newTopology, False)
            system = forcefield.createSystem(newTopology)

            # Set any previously existing atoms to be massless, they so won't move.

            for atom in existingAtomMap.values():
                system.setParticleMass(atom.index, 0.0)

            # If any heavy atoms were omitted, add them back to avoid steric clashes.

            nonbonded = [f for f in system.getForces() if isinstance(f, mm.CustomNonbondedForce)][0]
            for atom in self.topology.atoms():
                if atom.element not in (None, hydrogen) and atom not in existingAtomMap:
                    system.addParticle(0.0)
                    nonbonded.addParticle([])
                    newPositions.append(self.positions[atom.index])

            # For efficiency, only compute interactions that involve a new atom.

            nonbonded.addInteractionGroup([atom.index for atom in newAtoms], range(system.getNumParticles()))

            # Do an energy minimization.

            integrator = mm.LangevinIntegrator(300*unit.kelvin, 10/unit.picosecond, 5*unit.femtosecond)
            if seed is not None:
                integrator.setRandomNumberSeed(seed)
            if self.platform is None:
                context = mm.Context(system, integrator)
            else:
                context = mm.Context(system, integrator, self.platform)
            context.setPositions(newPositions)
            mm.LocalEnergyMinimizer.minimize(context)
            state = context.getState(getPositions=True)
            if newTopology.getNumResidues() > 1:
                # When looking for pairs of atoms that are too close to each other, exclude pairs that
                # are in the same residue or are directly bonded to each other.

                exclusions = dict((atom, {a.index for a in atom.residue.atoms()}) for atom in newAtoms)
                for a1, a2 in newTopology.bonds():
                    if a1 in exclusions:
                        exclusions[a1].add(a2.index)
                    if a2 in exclusions:
                        exclusions[a2].add(a1.index)
                cutoff = 0.13
                nearest = self._findNearestDistance(context, newAtoms, cutoff, exclusions)
                if nearest < cutoff:

                    # Some atoms are very close together.  Run some dynamics while slowly increasing the strength of the
                    # repulsive interaction to try to improve the result.

                    for i in range(10):
                        context.setParameter('C', 0.15*(i+1))
                        integrator.step(200)
                        d = self._findNearestDistance(context, newAtoms, cutoff, exclusions)
                        if d > nearest:
                            nearest = d
                            state = context.getState(getPositions=True)
                            if nearest >= cutoff:
                                break
                    context.setState(state)
                    context.setParameter('C', 1.0)
                    mm.LocalEnergyMinimizer.minimize(context)
                    state = context.getState(getPositions=True)

            # Now create a new Topology, including all atoms from the original one and adding the missing atoms.

            (newTopology2, newPositions2, newAtoms2, existingAtomMap2) = self._addAtomsToTopology(False, False)

            # Copy over the minimized positions for the new atoms.

            for a1, a2 in zip(newAtoms, newAtoms2):
                newPositions2[a2.index] = state.getPositions()[a1.index]
            self.topology = newTopology2
            self.positions = newPositions2
