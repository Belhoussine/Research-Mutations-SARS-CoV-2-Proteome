#!/bin/bash

# File : runEnergyMinimization.sh
# Author : Nurit Haspel (major)
# Author : Filip Jagodzinski (minor)
# Author : Katie Hursh (minor)
# Last Modify Date : 4 May 2018
#
# Description : script for energy minization of a PDB structure file
#
# Sample invocation:
#
# ./runEnergyMinimization.sh 1HHP.A34G
#
# which assumes that the processed.pdb.knr file exists
#
# param 1 : PDB ID.{target chain}{target residue}{target mutation}

if [ $# != 1 ] ; then
    echo "Usage: runEnergyMinimization.sh <pdbID>.<target chain><target res><target mut>"
    exit
fi

# location where script "files" are located
scriptlib="minimization/scripts"

# create variable for pdbFile, and fetch pdb file
pdbFileBase=$1
pdbFile=${pdbFileBase}.pdb

cp ../../${pdbFile} ./${pdbFile}

# check to make sure that file exists
if [ ! -e $pdbFile ] ; then
    echo "the file $pdbFile does not exist"
    echo "quiting"
    exit
fi

if [ ! -e $pdbFileBase ] ; then
    mkdir $pdbFileBase
fi

if [ ! -d $pdbFileBase ] ; then
    echo "couldn't created folder $pdbFileBase"
    echo "quiting"
    exit
fi

# process chain (this is internal to namd,
# and NOT the target chain)
chain=" "
echo "Chain:" $chain

# To do
# DETECT DISULFIDES
#

# extract chain and all HOH molecules
\grep " $chain " $pdbFile > chainAFirst.pdb
\grep -v "^HETATM" chainAFirst.pdb > chainA.pdb

# Check if first res. is GLY or PRO
nter="NTER"
firstRes=`grep ^ATOM $pdbFile | head -1 | awk '{print $4}'`
echo $firstRes
if [ "$firstRes" == "PRO" ] ;  then
    nter="PROP"
elif [ "$firstRes" == "GLY" ] ; then
    nter="GLYP"
fi
echo "nter:" $nter

# write PSF file
cat > newpsf.tcl << EOF
package require psfgen

resetpsf

topology $scriptlib/top_all27_prot_lipid_na.inp

pdbalias residue HIS HSD
pdbalias residue CYX CYS
pdbalias atom ILE CD1 CD

segment A {

first $nter
last CTER
pdb chainA.pdb

}
coordpdb chainA.pdb A
EOF

# disulfide and processed PDB files
#disuFile="M_C/data/${pdbcode}.${retained}${chainTarget}${resNumTarget}${mutantTarget}.mut.cur.out/${pdbcode}.${retained}${chainTarget}${resNumTarget}${mutantTarget}.DisulfideBonds.bnd.knr"
disuFile="testing.txt"
processedFile=$pdbFile
#"M_C/data/${pdbcode}.${retained}${chainTarget}${resNumTarget}${mutantTarget}.mut.cur.out/${pdbcode}.${retained}${chainTarget}${resNumTarget}${mutantTarget}.processed.pdb.knr"

echo "disuFile:" $disuFile
echo "processedFile:" $processedFile

if [ -e $disuFile ] ; then
    for lines in `cat $disuFile | awk '{print $2"and"$4}'` ; do
        atom1=`echo $lines | sed s/and/" "/ | awk '{print $1}'`
        atom2=`echo $lines | sed s/and/" "/ | awk '{print $2}'`
        res1=`grep ^ATOM $processedFile | awk -v lines=$atom1 '{if ($2 == lines) print $6}'`
        res2=`grep ^ATOM $processedFile | awk -v lines=$atom2 '{if ($2 == lines) print $6}'`
        echo "patch DISU $chain:$res1 $chain:$res2" >> newpsf.tcl
    done
fi
cat > newpsf1.tcl << EOF
guesscoord

writepsf ${pdbFileBase}_autopsf.psf

writepdb ${pdbFileBase}_autopsf.pdb

EOF

# concatenate
cat newpsf1.tcl >> newpsf.tcl ; rm -f newpsf1.tcl

# generate PSF
${scriptlib}/psfgen newpsf.tcl

# remove the intermediate chainA.pdb file that was created
rm -f chainA.pdb

# write namd param file
cat > namd-temp.namd << EOF
coordinates      ${pdbFileBase}_autopsf.pdb
structure        ${pdbFileBase}_autopsf.psf
paraTypeCharmm   on
parameters       $scriptlib/par_all27_prot_lipid_na.inp
exclude          scaled1-4

temperature		 0
cutoff			 12
dielectric		 1.0
switchdist		 10
outputEnergies   500
set ts 0
set outp 0
binaryoutput     no
binaryrestart off
restartname decoy_tmp
restartfreq 500
outputname out
#   minimization on
minimize 500
file copy -force decoy_tmp.coor ${pdbFileBase}_min.pdb
EOF

# invoke charmrun
${scriptlib}/charmrun ++local ${scriptlib}/namd2 namd-temp.namd > ${pdbFileBase}_min.log

# make the necessary output directories; copy files
mv ${pdbFileBase}_min.* ${pdbFileBase}/.
mv ${pdbFileBase}.pdb ${pdbFileBase}/${pdbFileBase}_preEM.pdb
mv ${pdbFileBase}/${pdbFileBase}_min.pdb ${pdbFileBase}/${pdbFileBase}.pdb

# rename all protonated histidines into HIS
sed -i 's/HSD/HIS/g' ${pdbFileBase}/${pdbFileBase}.pdb
sed -i 's/HSE/HIS/g' ${pdbFileBase}/${pdbFileBase}.pdb
sed -i 's/HSP/HIS/g' ${pdbFileBase}/${pdbFileBase}.pdb

# rename all CD of ILE to CD1
sed -i 's/CD  ILE/CD1 ILE/g' ${pdbFileBase}/${pdbFileBase}.pdb

# remove all of the H atoms, as they are added during curation
echo "" > ${pdbFileBase}/hRemoved.pdb
OLD_IFS="$IFS"
IFS=
while read line
do
    atom=${line: -2:-1}
    if [ ${atom} != "H" ]; then
	echo $line >> ${pdbFileBase}/hRemoved.pdb
    fi
done < ${pdbFileBase}/${pdbFileBase}.pdb
IFS="$OLD_IFS"

# copy the w/o H pdb file into the output file
cp ${pdbFileBase}/hRemoved.pdb ${pdbFileBase}/${pdbFileBase}.pdb

mv ${pdbFileBase}_autopsf* ${pdbFileBase}/.
mv decoy_* ${pdbFileBase}/.
mv out.* ${pdbFileBase}/.
mv namd-temp.namd ${pdbFileBase}/.
mv newpsf.tcl ${pdbFileBase}/.

cp ${pdbFileBase}/${pdbFileBase}.pdb ../../${pdbFileBase}_em.pdb
