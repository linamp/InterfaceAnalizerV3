import os
import re
import types
import numpy as np
from pymol import cmd
from math import cos, sin, radians

def createfrontbackRing(nSym, radius, rotation, namestring, cmd=cmd, saveflag=False):
    """Function to create a ring with monomer termini on the
    front and back surface of the ring
    nSym = circular symmetry definition (integer)
    radius = radius in Angstroms (float)
    rotation = array specifying the x,y,z rotation angles (in that order)
    namestring = Name of the model to be used as the monommer
    A PDB object with the namestring must exist in Pymol
    saveflag = Save a PDB containing two dimers (boolean)
    """
    
    # input checking
    if checkParams(nSym, radius, rotation, namestring, saveflag):
        print("There was an error with a parameter.  Please see")
        print("the above error message for how to fix it.")
        return None
        
    object_set = set(cmd.get_object_list())
    
    if namestring in object_set:
        cmd.create('dummy', namestring)
    else:
        raise NameError('Monomer object does not exist in Pymol')

    cmd.do('rotate x, 90, dummy, camera=0')
    cmd.do('rotate x, %s, dummy, camera=0'% rotation[0])
    cmd.do('rotate y, %s, dummy, camera=0'% rotation[1])
    cmd.do('rotate z, %s, dummy, camera=0'% rotation[2])

    # cmd.rotate('x', str(90), 'dummy', camera=0)
    # cmd.rotate('x', str(rotation[0]), 'dummy', camera=0) 
    # cmd.rotate('y', str(rotation[1]), 'dummy', camera=0)
    # cmd.rotate('z', str(rotation[2]), 'dummy', camera=0)

    adeg = 360/nSym
    arad = radians(adeg)

    monrange = np.arange(0,nSym)
    for ind in monrange:
        objname = 'su_'+ str(ind)
        cmd.create(objname, 'dummy')
        # cmd.rotate('z', str(adeg*ind), objname, camera=0)
        cmd.do('rotate z, %s, %s, camera=0'%(adeg*ind, objname))
        cmd.do('translate [%s, %s, 0], %s, camera=0'%(radius*cos(arad*ind), radius*sin(arad*ind), objname))
        # cmd.translate([radius*cos(arad*ind), radius*sin(arad*ind), 0], objname, camera=0)

    cmd.alter('%s and chain A'% 'su_0', 'chain="C"')
    cmd.alter('%s and chain B'% 'su_0', 'chain="D"')
    cmd.delete('dummy')
    
    if saveflag == True:
        savename = namestring + '_' + str(nSym) + '_' + str(radius)+'_'+str(rotation[0])+'_'+str(rotation[1])+'_'+str(rotation[2])
        cmd.save('%s.pdb'% savename, 'su_0 or su_1')
        cmd.delete('su_0 or su_1')
        return savename
    
    cmd.extend("createfrontbackRing", createfrontbackRing)

def checkParams(nSym, radius, rotation, namestring, cmd=cmd, saveflag=False):
    """
    Check input format validity
    """ 
    if not type(nSym) is int:
        raise TypeError("Symmetry number must be an integer")
         
    if not type(radius) is float:
        raise TypeError("Ring radius must be a float")
    
    assert isinstance(namestring, str), "Namestring must be a string with the model name"
        
    if len(namestring) == 0:
        print("Error: namestring must have non-zero length")
        return False
        
    if not type(saveflag) is bool:
        raise TypeError("saveflag must be True or False")
        
    if not len(rotation) == 3:
        raise ValueError("The rotation variable must specify three angles")
       

    
