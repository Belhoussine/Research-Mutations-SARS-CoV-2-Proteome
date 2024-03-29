package require psfgen

resetpsf

topology minimization/scripts/top_all27_prot_lipid_na.inp

pdbalias residue HIS HSD
pdbalias residue CYX CYS
pdbalias atom ILE CD1 CD

segment A {

first NTER
last CTER
pdb chainA.pdb

}
coordpdb chainA.pdb A
guesscoord

writepsf 6M0J.E492C_autopsf.psf

writepdb 6M0J.E492C_autopsf.pdb

