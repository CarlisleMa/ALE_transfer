    def _addIons(self, forcefield, numWaters, replaceableMols, ionCutoff=0.5*nanometer, positiveIon='Na+', negativeIon='Cl-', ionicStrength=0*molar, neutralize=True, residueTemplates=dict()):
        """Adds ions to the system by replacing certain molecules.

        Parameters
        ----------
        forcefield : ForceField
            the ForceField to use to determine the total charge of the system.
        numWaters : int
            the total number of water molecules in the simulation box, used to
            calculate the number of ions / concentration to add.
        replaceableMols : dict
            the molecules to replace by ions, as a dictionary of residue:positions
        ionCutoff: distance=0.5*nanometer
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

        if not neutralize and ionicStrength == 0 *molar:
            return

        posIonElements = {'Cs+': elem.cesium, 'K+': elem.potassium,
                          'Li+': elem.lithium, 'Na+': elem.sodium,
                          'Rb+': elem.rubidium}
        negIonElements = {'Cl-': elem.chlorine, 'Br-': elem.bromine,
                          'F-': elem.fluorine, 'I-': elem.iodine}

        ionPositions = []

        numReplaceableMols = len(replaceableMols)

        # Fetch ion elements from user input
        if positiveIon not in posIonElements:
            raise ValueError('Illegal value for positive ion: {}'.format(positiveIon))
        if negativeIon not in negIonElements:
            raise ValueError('Illegal value for negative ion: {}'.format(negativeIon))
        positiveElement = posIonElements[positiveIon]
        negativeElement = negIonElements[negativeIon]

        # Determine the total charge of the system
        system = forcefield.createSystem(self.topology, residueTemplates=residueTemplates)
        for i in range(system.getNumForces()):
            if isinstance(system.getForce(i), (NonbondedForce, AmoebaMultipoleForce)):
                nonbonded = system.getForce(i)
                break
        else:
            raise ValueError('The ForceField does not specify a NonbondedForce')

        totalCharge = 0.0
        if isinstance(nonbonded, AmoebaMultipoleForce):
            for i in range(nonbonded.getNumMultipoles()):
                totalCharge += nonbonded.getMultipoleParameters(i)[0].value_in_unit(elementary_charge)
        else:
            for i in range(nonbonded.getNumParticles()):
                totalCharge += nonbonded.getParticleParameters(i)[0].value_in_unit(elementary_charge)

        # Round to nearest integer
        totalCharge = int(floor(0.5 + totalCharge))

        # Figure out how many ions to add based on requested params/concentration
        numPositive, numNegative = 0, 0
        if neutralize:
            if abs(totalCharge) > numReplaceableMols:
                raise Exception('Cannot neutralize the system because the charge is greater than the number of available positions for ions')
            if totalCharge > 0:
                numNegative += totalCharge
            else:
                numPositive -= totalCharge

        if ionicStrength > 0 * molar:
            numIons = (numWaters - numPositive - numNegative) * ionicStrength / (55.4 * molar)  # Pure water is about 55.4 molar (depending on temperature)
            numPairs = int(floor(numIons + 0.5))
            numPositive += numPairs
            numNegative += numPairs
        totalIons = numPositive + numNegative

        if totalIons > 0:
            # Randomly select a set of waters
            # while ensuring ions are not placed too close to each other.
            modeller = Modeller(self.topology, self.positions)

            replaceableList = list(replaceableMols.keys())
            numAddedIons = 0
            numTrials = 10  # Attempts to add ions N times before quitting
            toReplace = []  # list of molecules to be replaced
            while numAddedIons < totalIons:
                pickedMol = random.choice(replaceableList)
                replaceableList.remove(pickedMol)
                # Check distance to other ions
                for pos in ionPositions:
                    distance = norm(pos - replaceableMols[pickedMol])
                    if distance <= ionCutoff:
                        numTrials -= 1
                        break
                else:
                    toReplace.append(pickedMol)
                    ionPositions.append(replaceableMols[pickedMol])
                    numAddedIons += 1

                    n_trials = 10

                if n_trials == 0:
                    raise ValueError('Could not add more than {} ions to the system'.format(numAddedIons))

            # Replace waters/ions in the topology
            modeller.delete(toReplace)
            ionChain = modeller.topology.addChain()
            for i, water in enumerate(toReplace):
                element = (positiveElement if i < numPositive else negativeElement)
                newResidue = modeller.topology.addResidue(element.symbol.upper(), ionChain)
                modeller.topology.addAtom(element.symbol, element, newResidue)
                modeller.positions.append(replaceableMols[water])

            # Update topology/positions
            self.topology = modeller.topology
            self.positions = modeller.positions
