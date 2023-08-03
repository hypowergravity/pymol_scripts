from pymol import cmd, stored, math
import numpy as np

def rmsf2B (pdb,source="rmsf.dat",max_b=None,min_b=None):
    """ This script replaces the rmsf value of the data generated from cpptraj and replaces in B-factor column of the pdb object.
  usage: rmsf2B pdb,source="rmsf.dat"
  pdb = any object selection with name 
  source = source file generated from cpptraj like what generated from cpptraj (atomicfluct byres out rmsf.dat :1-480&!@H*= bfac),
  this will only work on average pdb generated from cpptraj (average average.pdb pdb), i.e with only single chain.
  max_b,min_b can be set the same scale if multiple rmsf files to be used
  example: rmsf2B amylase, rmsf.dat"""
    obj=cmd.get_object_list(pdb)[0]
    data = np.loadtxt(source)
    residue_number = data.T[0].astype(int)
    b_factor = data.T[1].astype(float)
    max_b_factor = max_b if max_b is not None else max(b_factor)
    min_b_factor = min_b if min_b is not None else min(b_factor)
    res_2_set = dict(zip(residue_number, b_factor))
    for resi,b_fac in res_2_set.items():
        cmd.alter("resi %s and n. CA"%(resi), "b=%s"%b_fac)

    cmd.show_as("cartoon",pdb)
    cmd.show("sticks","org")
    cmd.hide("sticks","hydrogen")
    cmd.cartoon("putty", pdb)
    cmd.set("cartoon_putty_scale_min", min_b_factor,obj)
    cmd.set("cartoon_putty_scale_max", max_b_factor,obj)
    cmd.set("cartoon_putty_transform", 0,pdb)
    cmd.set("cartoon_putty_radius", 0.1,pdb)
    cmd.spectrum("b","rainbow", "%s and n. CA " %obj)
    cmd.ramp_new("count", obj, [min_b_factor, max_b_factor], "rainbow")
    cmd.recolor()

cmd.extend("rmsf2B", rmsf2B);
