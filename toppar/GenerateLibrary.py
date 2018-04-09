#!/usr/bin/python

import sys, os
sys.path.append (os.path.join (os.environ["HOME"], "devel", "MolarisTools"))

from MolarisTools.Library import AminoLibrary


library = AminoLibrary (filename="top_all27_prot_na.inp", topologyFormat="CHARMM")

aminoacids = ("ALA",  "ARG",  "ASN",  "ASP",  "CYS",  "GLN",  "GLU",  "GLY",  "HSD",  
              "HSE",  "HIS",  "ILE",  "LEU",  "LYS",  "MET",  "PHE",  "PRO",  "SER",  
              "THR",  "TRP",  "TYR",  "VAL", )
for label in aminoacids:
    filename = "%s.ic" % label
    if os.path.exists (filename):
        continue

    component = library[label]
    component.WriteInternalCoordinates ()
