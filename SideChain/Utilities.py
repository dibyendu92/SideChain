#-------------------------------------------------------------------------------
# . File      : Utilities.py
# . Program   : SideChain
# . Copyright : USC, Mikolaj Feliks (2018)
# . License   : GNU GPL v3.0       (http://www.gnu.org/licenses/gpl-3.0.en.html)
#-------------------------------------------------------------------------------
import numpy

def VectorNormalize (vector):
    return (1.0 / numpy.linalg.norm (vector)) * vector


#===============================================================================
# . Main program
#===============================================================================
if (__name__ == "__main__"): pass