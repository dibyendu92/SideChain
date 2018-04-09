#-------------------------------------------------------------------------------
# . File      : ProteinModel.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import exceptions, math, numpy


_BACKBONE_ATOMS = ("N", "CA", "C", "O", "+N", "+CA", "-C", )
_SKIP_RESIDUES = ("HOH", "WAT", )

_PDB_FORMAT_ATOM   = "%-6s%5d %-4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%22s\n"


class ProteinContainer (object):
    def __init__ (self, **kwargs):
        for (key, value) in kwargs.items ():
            setattr (self, key, value)


class ProteinChain (ProteinContainer):
    def AddResidue (self, residue):
        if not hasattr (self, "residues"):
            self.residues = []
        self.residues.append (residue)


class ProteinResidue (ProteinContainer):
    def AddAtom (self, atom, reorder=False):
        if not hasattr (self, "atoms"):
            self.atoms = []
        if (atom.label not in self):
            if (reorder):
                atom.serial = self.atoms[-1].serial + 1
            self.atoms.append (atom)

    def Label (self):
        return "%s.%s.%d" % (self.parent.label, self.label, self.serial)

    def __contains__ (self, label):
        if hasattr (self, "atoms"):
            for atom in self.atoms:
                if (label == atom.label):
                    return True
        return False

    def __getitem__ (self, label):
        if hasattr (self, "atoms"):
            for atom in self.atoms:
                if (label == atom.label):
                    return atom
        raise exceptions.StandardError ("Atom %s not found in residue %s.%d.%s" % (label, self.parent.label, self.serial, self.label))


class ProteinAtom (ProteinContainer):
    @property
    def coordinates (self):
        return numpy.array ((self.x, self.y, self.z))


class ProteinModel (object):
    def __init__ (self, pdb, library, internal, parameters, mutations, verbose=True):
        """Constructor."""
        self.pdb = pdb
        self.library = library
        self.internal = internal
        self.parameters = parameters
        self.mutations = mutations
        self.verbose = verbose
        self.prompt = "# . ProteinModel> "

    @property
    def nchains (self):
        if hasattr (self, "chains"):
            return len (self.chains)
        return 0

    @property
    def nresidues (self):
        if hasattr (self, "chains"):
            nresidues = 0
            for chain in self.chains:
                nresidues += chain.nresidues
            return nresidues
        return 0

    @property
    def natoms (self):
        if hasattr (self, "chains"):
            natoms = 0
            for chain in self.chains:
                for residue in chain.residues:
                    natoms += residue.natoms
            return natoms
        return 0

    def _Write (self, message):
        if (self.verbose):
            if hasattr (self, "prompt"):
                print ("%s%s" % (self.prompt, message))


    def _FindMutation (self, label, serial):
        for (i, chain) in enumerate (self.chains):
            for (j, residue) in enumerate (chain.residues):
                if ((chain.label == label) and (residue.serial == serial)):
                    return (i, j)
        return None


    def _CheckIfMutated (self, label, serial):
        for (chainLabel, residueSerial, residueTargetLabel) in self.mutations:
            if ((chainLabel == label) and (residueSerial == serial)):
                return residueTargetLabel
        return None


    def Build (self, skipChains=None, skipResidues=_SKIP_RESIDUES):
        """Builds a protein model.

        Residues that are subject to mutation start with empty side-chains."""
        chainLabels = set ()
        for pdbResidue in self.pdb.residues:
            if (skipChains != None):
                if (pdbResidue.chain in skipChains):
                    continue
            chainLabels.add (pdbResidue.chain)

        self.chains = []
        self.pairs = []

        for (i, chainLabel) in enumerate (chainLabels):
            chain = ProteinChain (label=chainLabel)

            for (j, pdbResidue) in enumerate (self.pdb.residues):
                if (pdbResidue.label in skipResidues):
                    continue

                if (pdbResidue.chain == chainLabel):
                    isMutated = self._CheckIfMutated (chainLabel, pdbResidue.serial)
                    if (isMutated):
                        pair = (i, j)
                        self.pairs.append (pair)
                        targetLabel = isMutated
                        residue = ProteinResidue (label = targetLabel, 
                                                  serial = pdbResidue.serial, 
                                                  parent = chain )

                        for pdbAtom in pdbResidue.atoms:
                            if (pdbAtom.label in _BACKBONE_ATOMS):
                                atom = ProteinAtom (label = pdbAtom.label, 
                                                    serial = pdbAtom.serial, 
                                                    x = pdbAtom.x, 
                                                    y = pdbAtom.y, 
                                                    z = pdbAtom.z, 
                                                    parent = residue )
                                residue.AddAtom (atom)
                        chain.AddResidue (residue)
                    else:
                        residue = ProteinResidue (label = pdbResidue.label, 
                                                  serial = pdbResidue.serial, 
                                                  parent = chain )
                        for pdbAtom in pdbResidue.atoms:
                            atom = ProteinAtom (label = pdbAtom.label, 
                                                serial = pdbAtom.serial, 
                                                x = pdbAtom.x, 
                                                y = pdbAtom.y, 
                                                z = pdbAtom.z, 
                                                parent = residue )
                            residue.AddAtom (atom)
                        chain.AddResidue (residue)
            self.chains.append (chain)

        for chain in self.chains:
            self._Write ("Chain: %s (%d residues)" % (chain.label, len (chain.residues)))
        self._Write ("Created model chain:residue:protein")


    def AddHydrogens (self):
        """Adds missing hydrogens."""

        for chain in self.chains:
            for residue in chain.residues:
                if (residue.label not in self.library):
                    raise exceptions.StandardError ("Residue %s not found in the amino-acid library." % residue.label)
                component = self.library[residue.label]

                hydrogens = []
                for atom in component.atoms:
                    if (atom.atomLabel[0] == "H"):
                        if (atom.atomLabel not in residue):
                            hydrogens.append (atom.atomLabel)

                if (residue.label not in self.internal):
                    raise exceptions.StandardError ("Residue %s not found in the internal coordinate library." % residue.label)
                ictable = self.internal[residue.label]

                for hydrogen in hydrogens:
                    for (i, j, k, l, distanceLeft, angleLeft, torsion, angleRight, distanceRight, improper) in ictable:
                        if (hydrogen == l):
                            internal = (i, j, k, distanceRight, angleRight, torsion, improper)
                            self._AddAtom (residue, hydrogen, internal)


    def _AddAtom (self, residue, label, internal):
        (a, b, c, distance, angle, torsion, improper) = internal

        pa = residue[a].coordinates
        pb = residue[b].coordinates
        pc = residue[c].coordinates
        
        CalculatePosition = self._CalculatePositionImproper if (improper) else self._CalculatePositionNormal
        pd = CalculatePosition (pa, pb, pc, distance, angle, torsion)
        
        atom = ProteinAtom (label=label, serial=-1, x=pd[0], y=pd[1], z=pd[2], parent=residue)
        residue.AddAtom (atom, reorder=True)

        self._Write ("Added atom %s to residue %s" % (label, residue.Label ()))


    @staticmethod
    def _UnmetDependencies (sequence, labels):
        for dependent in labels:
            if not ((dependent in sequence) or (dependent in _BACKBONE_ATOMS)):
                return True
        return False

    def _BuildSequence_r (self, internalCoordinates, connectivityTable, label, sequence, unmet):
        connections = connectivityTable[label]

        for otherLabel in connections:
            if ((otherLabel in _BACKBONE_ATOMS) or (otherLabel in sequence)):
                continue

            internal = self._SearchInternal (internalCoordinates, otherLabel)
            (a, b, c, distance, angle, torsion, improper) = internal

            if (self._UnmetDependencies (sequence, (a, b, c))):
                if (otherLabel not in unmet):
                    unmet.append (otherLabel)
                continue

            sequence.append (otherLabel)
            if (otherLabel in unmet):
                unmet.remove (otherLabel)

            self._BuildSequence_r (internalCoordinates, connectivityTable, otherLabel, sequence, unmet)


    def Mutate (self):
        """Generates initial conformations of mutated residues."""
        for (chainLabel, residueSerial, residueTargetLabel) in self.mutations:
            found = self._FindMutation (chainLabel, residueSerial)
            if (not found):
                raise exceptions.StandardError ("Residue %s.%d not found in the protein." % (chainLabel, residueSerial))
            (i, j) = found
            chain = self.chains[i]
            residue = chain.residues[j]

            if (residue.label not in self.library):
                raise exceptions.StandardError ("Residue %s not found in the amino-library." % residue.label)

            component = self.library[residue.label]
            component.GenerateConnectivities ()

            if (residue.label not in self.internal):
                raise exceptions.StandardError ("Residue %s not found in the internal coordinate library." % residue.label)
            ictable = self.internal[residue.label]

            #self._AddAtomInternal_r (ictable, component.connectivity, residue, "CA")

            rootLabel = "CA"
            sequence = []
            unmet = []
            # self._BuildSequence_r (ictable, component.connectivity, rootLabel, sequence, unmet)
            # 
            # self._Write ("Adding atoms for %s in sequence: %s" % (residue.Label (), " ".join (sequence)))
            # 
            # for label in sequence:
            #     internal = self._SearchInternal (ictable, label)
            #     self._AddAtom (residue, label, internal)

            for (i, j, k, l, distanceLeft, angleLeft, torsion, angleRight, distanceRight, improper) in ictable:
                tests = []
                for label in (i, j, k, l):
                    tests.append ((label[0] == "-") or (label[0] == "+"))
                if any (tests):
                    continue
                if (l not in residue):
                    internal = (i, j, k, distanceRight, angleRight, torsion, improper)
                    self._AddAtom (residue, l, internal)


    @staticmethod
    def _SearchInternal (internalCoordinates, label):
        for (i, j, k, l, distanceLeft, angleLeft, torsion, angleRight, distanceRight, improper) in internalCoordinates:
            if (improper):
                if (l == label):
                    return (i, j, k, distanceRight, angleRight, torsion, improper)
                elif (i == label):
                    return (l, j, k, distanceLeft, angleLeft, torsion, improper)
            else:
                if (l == label):
                    return (i, j, k, distanceRight, angleRight, torsion, improper)
                elif (i == label):
                    return (l, k, j, distanceLeft, angleLeft, torsion, improper)
        raise exceptions.StandardError ("Atom %s not found in internal coordinates." % label)

    @staticmethod
    def _SearchFormatInternal (internalCoordinates, label):
        entries = []
        for (i, j, k, l, distanceLeft, angleLeft, torsion, angleRight, distanceRight, improper) in internalCoordinates:
            entry = None
            if (improper):
                if (l == label):
                    entry = (i, j, k, distanceRight, angleRight, torsion, improper)
                elif (i == label):
                    entry = (l, j, k, distanceLeft, angleLeft, torsion, improper)
            else:
                if (l == label):
                    entry = (i, j, k, distanceRight, angleRight, torsion, improper)
                elif (i == label):
                    entry = (l, k, j, distanceLeft, angleLeft, torsion, improper)
            if (entry != None):
                entries.append (entry)
        if (entries is []):
            raise exceptions.StandardError ("Atom %s not found in internal coordinates." % otherLabel)
        return entries

    @staticmethod
    def _CheckAtomsPresent (residue, entries):
        for i, (a, b, c, distance, angle, torsion, improper) in enumerate (entries):
            checks = (a in residue, b in residue, c in residue)
            if all (checks):
                return i
        return -1


    def _AddAtomInternal_r (self, internalCoordinates, connectivityTable, residue, label):
        """Adds an atom to a residue based on internal coordinates.

        The format of internal coordinates is described in:
        https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/intcor/"""

        connections = connectivityTable[label]

        for otherLabel in connections:
            if (otherLabel not in residue):
                entries = self._SearchFormatInternal (internalCoordinates, otherLabel)
    
                i = self._CheckAtomsPresent (residue, entries)
                if (i < 0):
                    self._Write ("Warning: skipping atom %s." % otherLabel)
                    continue
                self._AddAtom (residue, label, entries[i])
                self._Write ("Added atom %s to residue %s" % (otherLabel, residue.Label ()))
            
                self._AddAtomInternal_r (internalCoordinates, connectivityTable, residue, otherLabel)


    @staticmethod
    def _VectorNormalize (vector):
        return (1.0 / numpy.linalg.norm (vector)) * vector


    def _CalculatePositionImproper (self, a, b, c, Rcd, Tbcd, Pabcd):
    #                I (a)    L (d)
    #                 \     /
    #                  \   /
    #                   *K (c)
    #                   |
    #                   |
    #                   J (b)
    #        values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)
        ca = a - c
        j  = self._VectorNormalize (b - c)
        k  = self._VectorNormalize (numpy.cross (ca, j))
        i  = self._VectorNormalize (numpy.cross (j, k))
    
        psi = Pabcd * math.pi / 180.0
        t   = math.cos (psi) * i + math.sin (psi) * k
    
        chi = Tbcd * math.pi / 180.0
        q   = math.cos (chi) * j + math.sin (chi) * t
    
        return (c + Rcd * q)
    
    
    def _CalculatePositionNormal (self, a, b, c, Rcd, Tbcd, Pabcd):
    #            (a) I
    #                 \
    #                  \
    #               (b) J----K (c)
    #                         \
    #                          \
    #                           L (d)
    #        values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)
        ab  = a - b
        j   = self._VectorNormalize (c - b)
        k   = self._VectorNormalize (numpy.cross (ab, j))
        i   = self._VectorNormalize (numpy.cross (j, k))
    
        psi = Pabcd * math.pi / 180.0
        t   = math.cos (psi) * i + math.sin (psi) * k
   
        chi = Tbcd * math.pi / 180.0
        q   = -math.cos (chi) * j + math.sin (chi) * t
 
        return (c + Rcd * q)


    def BuildEnergyModel (self):
        """Builds an energy model."""
        # try:
        #     component = components[residue.label]
        # except exceptions.KeyError:
        #     component = self.library[residue.label]
        #     self._Write ("Picked %s." % residue.label)
        # 
        #     component.GenerateAngles ()
        #     component.GenerateTorsions ()
        #     components[residue.label] = component


        # (torsionTypes, foo, bar) = component._TorsionsToTypes ()
        # for (torsion, types) in zip (component.torsions, torsionTypes):
        #     (typea, typeb, typec, typed) = types
        #     torsionalParameters = self.parameters.GetTorsion (typeb, typec)
        #     if (not torsionalParameters):
        #         raise exceptions.StandardError ("Parameters for torsion X-%s-%s-X not found." % (typeb, typec))

        # missingLabels = []
        # for atom in component.atoms:
        #         if (atom.atomLabel not in residue):
        #             missingLabels.append (atom.atomLabel)
        # if (missingLabels):
        #     self._Write ("Residue %s.%s.%d is missing atoms: %s" % (residue.parent.label, residue.label, residue.serial, " ".join (missingLabels)))
        pass


    def Optimize (self):
        """Optimizes side-chain positions with Monte Carlo."""
        if hasattr (self, "energyModel"):
            pass


    def DumpAsXYZ (self, filename="mutant.xyz"):
        output = open (filename, "w")
        output.write ("%d\n...\n" % self.natoms)
        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    output.write ("%4s  %8.3f  %8.3f  %8.3f\n" % (atom.label, atom.x, atom.y, atom.z))
        output.close ()


    def DumpAsPDB (self, filename="mutant.pdb"):
        output = open (filename, "w")
        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    output.write (_PDB_FORMAT_ATOM % ("ATOM", atom.serial, atom.label, residue.label, chain.label, residue.serial, "", atom.x, atom.y, atom.z, "SEGM"))
        output.close ()


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass
