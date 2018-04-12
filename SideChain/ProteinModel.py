#-------------------------------------------------------------------------------
# . File      : ProteinModel.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import exceptions, math, numpy, random

import Utilities

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

    def GetIndex (self, label):
        if (label in self):
            for (i, atom) in enumerate (self.atoms):
                if (label == atom.label):
                    return i
        return -1

    def LabelsToIndices (self, labels):
        indices = []
        for (i, atom) in enumerate (self.atoms):
            if (atom.label in labels):
                indices.append (i)
        return indices

    @property
    def natoms (self):
        if hasattr (self, "atoms"):
            return len (self.atoms)
        return 0


class ProteinAtom (ProteinContainer):
    @property
    def coordinates (self):
        return numpy.array ((self.x, self.y, self.z))

    @coordinates.setter
    def coordinates (self, array):
        self.x = array[0]
        self.y = array[1]
        self.z = array[2]


class RotatableTorsion (ProteinContainer):
    def __init__ (self, atoms, parent):
        self.atoms = atoms
        self.parent = parent
#        a   = atoms[0].coordinates
#        b   = atoms[1].coordinates
#        c   = atoms[2].coordinates
#        ab  = a - b
#        j   = Utilities.VectorNormalize (c - b)
#        k   = Utilities.VectorNormalize (numpy.cross (ab, j))
#        i   = Utilities.VectorNormalize (numpy.cross (j, k))

    def Rotate (self, degree):
        pass


class ProteinModel (object):
    def __init__ (self, pdb, library, internal, rotatable, parameters, mutations, verbose=True):
        """Constructor."""
        self.pdb = pdb
        self.library = library
        self.internal = internal
        self.rotatable = rotatable
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


    def _CalculateResidueAtoms (self, pdbResidue, skipChains, skipResidues):
        natoms=0
        if (pdbResidue.chain not in skipChains):
            if (pdbResidue.label not in skipResidues):
                isMutated = self._CheckIfMutated (pdbResidue.chain, pdbResidue.serial)
                if (isMutated):
                    for pdbAtom in pdbResidue.atoms:
                        if (pdbAtom.label in _BACKBONE_ATOMS):
                            natoms+=1
                    targetLabel = isMutated
                    component = self.library[targetLabel]
                    for atom in component.atoms:
                        if (atom.atomLabel not in _BACKBONE_ATOMS):
                            natoms+=1
                else:
                    natoms+=len (pdbResidue.atoms)
        return natoms


    def Build (self, skipChains=None, skipResidues=_SKIP_RESIDUES):
        """Builds a protein model.

        Residues that are subject to mutation start with empty side-chains."""
        # natoms=0
        # for residue in self.pdb.residues:
        #     natoms+=self._CalculateResidueAtoms (residue, skipChains, skipResidues)
        # self.coordinates = numpy.empty ((natoms, 3))
        # 
        # self._Write ("Array allocated for %d atoms." % natoms)

        chainLabels = set ()
        for pdbResidue in self.pdb.residues:
            if (skipChains != None):
                if (pdbResidue.chain in skipChains):
                    continue
            chainLabels.add (pdbResidue.chain)

        self.chains = []
        self.pairs = []
        begin = 0

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
                                                  begin = begin, 
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
                                                  begin = begin, 
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
        """Adds missing atoms to mutated residues based on internal coordinates.

        The format of internal coordinates is described in:
        https://www.charmm.org/charmm/documentation/by-version/c40b1/params/doc/intcor/"""

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

            # rootLabel = "CA"
            # sequence = []
            # unmet = []
            # self._BuildSequence_r (ictable, component.connectivity, rootLabel, sequence, unmet)
            # 
            # self._Write ("Adding atoms for %s in sequence: %s" % (residue.Label (), " ".join (sequence)))
            # 
            # for label in sequence:
            #     internal = self._SearchInternal (ictable, label)
            #     self._AddAtom (residue, label, internal)

            # . Internal coordinate tables from CHARMM seem to already contain atoms in "safe" sequences.

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
    def _CheckAtomsPresent (residue, entries):
        for i, (a, b, c, distance, angle, torsion, improper) in enumerate (entries):
            checks = (a in residue, b in residue, c in residue)
            if all (checks):
                return i
        return -1

    @staticmethod
    def _CalculatePositionImproper (a, b, c, Rcd, Tbcd, Pabcd):
    #                I (a)    L (d)
    #                 \     /
    #                  \   /
    #                   *K (c)
    #                   |
    #                   |
    #                   J (b)
    #        values (Rik),(Tikj),(Pijkl),T(jkl),(Rkl)
        ca = a - c
        j  = Utilities.VectorNormalize (b - c)
        k  = Utilities.VectorNormalize (numpy.cross (ca, j))
        i  = Utilities.VectorNormalize (numpy.cross (j, k))
    
        psi = Pabcd * math.pi / 180.0
        t   = math.cos (psi) * i + math.sin (psi) * k
    
        chi = Tbcd * math.pi / 180.0
        q   = math.cos (chi) * j + math.sin (chi) * t
    
        return (c + Rcd * q)
    
    @staticmethod
    def _CalculatePositionNormal (a, b, c, Rcd, Tbcd, Pabcd):
    #            (a) I
    #                 \
    #                  \
    #               (b) J----K (c)
    #                         \
    #                          \
    #                           L (d)
    #        values (Rij),(Tijk),(Pijkl),(Tjkl),(Rkl)
        ab  = a - b
        j   = Utilities.VectorNormalize (c - b)
        k   = Utilities.VectorNormalize (numpy.cross (ab, j))
        i   = Utilities.VectorNormalize (numpy.cross (j, k))
    
        psi = Pabcd * math.pi / 180.0
        t   = math.cos (psi) * i + math.sin (psi) * k
   
        chi = Tbcd * math.pi / 180.0
        q   = -math.cos (chi) * j + math.sin (chi) * t
 
        return (c + Rcd * q)


    def _BuildEnergyModel (self):
        """Builds an energy model."""
        pass

    def _IdentifyRotatableTorsions (self):
        for (i, j) in self.pairs:
            chain = self.chains[i]
            residue = chain.residues[j]
            self._Write ("Generating torsions for residue %s..." % residue.Label ())

            component = self.library[residue.label]
            component.GenerateAngles ()
            component.GenerateTorsions ()
            component.GenerateConnectivities ()

            rotatableTorsions = self.rotatable[residue.label]
            self.torsions = []

            for (labelb, labelc) in rotatableTorsions:
                roots = []
                rotatableLabels = []
                self._GetRotatableAtoms_r (labelb, labelc, component.connectivity, roots, rotatableLabels)
                #rotatable = residue.LabelsToIndices (rotatableLabels)

                self._Write ("Torsion X-%s-%s-X rotatable atoms: %s" % (labelb, labelc, " ".join (rotatableLabels)))
                self._Write ("Collecting energy terms:")

                terms = []
                for (a, b, c, d) in component.torsions:
                    if ((b == labelb) and (c == labelc)) or ((b == labelc) and (c == labelb)):
                        terms.append ((a, b, c, d))
                        self._Write ("  %s-%s-%s-%s" % (a, b, c, d))

                (atomb, atomc) = (component[labelb], component[labelc])
                (tb, tc) = (atomb.atomType, atomc.atomType)
                torsionalParameters = self.parameters.GetTorsion (tb, tc)
                #Torsion(typeb='C4', typec='C4', k=1.5, periodicity=3, phase=0.0)

                if (torsionalParameters is None):
                    raise exceptions.StandardError ("Parameters for torsion X-%s-%s-X not found." % (tb, tc))

                connections = component.connectivity[labelb]
                for connection in connections:
                    if (connection != labelc):
                        break
                labela = connection

                atoms = [residue[labela], residue[labelb], residue[labelc], ]
                for label in rotatableLabels:
                    atoms.append (residue[label])

                torsion = RotatableTorsion (atoms=atoms, parent=residue)
                self.torsions.append (torsion) 


    def _GetRotatableAtoms_r (self, precedingLabel, rootLabel, connectivityTable, roots, atoms):
        connections = connectivityTable[rootLabel]
        roots.append (rootLabel)

        for connection in connections:
            if (connection != precedingLabel):
                atoms.append (connection)
                if (connection not in roots):
                    self._GetRotatableAtoms_r (rootLabel, connection, connectivityTable, roots, atoms)


    def Optimize (self):
        """Optimizes side-chain positions with Monte Carlo."""
        #if not hasattr (self, "energyModel"):
        #    self._BuildEnergyModel ()

        self._IdentifyRotatableTorsions ()


    def SavePDB (self, filename="mutant.pdb"):
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
