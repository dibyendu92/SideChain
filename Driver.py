#!/usr/bin/python

import sys, os

sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools.Parser import PDBFile
from MolarisTools.Library import ParametersLibrary, AminoLibrary

from SideChain import ProteinModel, InternalLibrary, RotatableBonds


#-------------------------------------------------------------------------------
library = AminoLibrary (filename=os.path.join ("toppar", "amino98_custom.lib"))
internalLibrary = InternalLibrary ()
rotatableBonds = RotatableBonds ()

pdb = PDBFile ("4pti.pdb")
parameters = ParametersLibrary (filename=os.path.join ("toppar", "parm.lib"))

mutations = (("A", 50, "TRP"), )    # ("A", 33, "MET"), 
protein = ProteinModel (pdb, library, internalLibrary, parameters, mutations)

protein.Build ()
protein.Mutate ()
protein.Optimize ()
protein.DumpAsPDB ()
