#!/usr/bin/python
#-------------------------------------------------------------------------------
# . File      : Driver.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import sys, os

sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools.Parser import PDBFile
from MolarisTools.Library import ParametersLibrary, AminoLibrary

from SideChain import ProteinModel, InternalLibrary, RotatableBonds


library = AminoLibrary (filename=os.path.join ("toppar", "amino98_custom.lib"))
internalLibrary = InternalLibrary ()
rotatableBonds = RotatableBonds ()

pdb = PDBFile ("4pti.pdb")
parameters = ParametersLibrary (filename=os.path.join ("toppar", "parm.lib"))

mutations = (("A", 50, "GLU"), )    # ("A", 33, "MET"), 
protein = ProteinModel (pdb, library, internalLibrary, rotatableBonds, parameters, mutations)

protein.Build ()
protein.Optimize ()
protein.SavePDB ()

values = (0.0, 5.0, 10.0, 15.0, 20.0, )
torsion = protein.torsions[0]
torsion.Print ()

for (i, value) in enumerate (values):
    torsion.Rotate (degree=value)
    protein.SavePDB (filename="test%02d.pdb" % i)
