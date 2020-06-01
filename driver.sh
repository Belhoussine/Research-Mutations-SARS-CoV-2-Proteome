#!/bin/bash

# After making promute, copy this script into the build directory
# Change the "FILES" variable to the script (or scripts) you would like to execute
# to profile different mutation scripts, you can run this script with the "time" command 

# generate list of mutations
./proMuteBatch 6M0J E:E 493:493 X mutList

declare -A residues
residues+=( 
  ["PHE"]="F" 
  ["LEU"]="L" 
  ["ILE"]="I"
  ["MET"]="M"
  ["VAL"]="V"
  ["SER"]="S"
  ["PRO"]="P"
  ["THR"]="T"
  ["ALA"]="A"
  ["TYR"]="Y"
  ["HIS"]="H"
  ["GLN"]="Q"
  ["ASN"]="N"
  ["LYS"]="K"
  ["ASP"]="D"
  ["GLU"]="E"
  ["CYS"]="C"
  ["TRP"]="W"
  ["ARG"]="R"
  ["GLY"]="G"
)

# append the correct flags
sed -i 's/$/ em/' mutList

# Customize with your scripts
FILES=("mutList")

# Execute the lines in the script
for FILE in $FILES; do


  # Create a file to store SDM script
  SDMFILE="${FILE}_sdm"
  touch $SDMFILE

  # Read each mutation into an array
  readarray -t LINES < "$FILE"

  # Routine to execute for every mutation
  for LINE in "${LINES[@]}"; do

    # Make a directory to store outputs
    # Structure is PDBID.Chain.Location.Res2.output
    OUTPUTDIRNAME=$(echo $LINE | awk '{ print $2"."$3"."$4"."$5".output" }')
    echo [INFO] Creating Directory $OUTPUTDIRNAME
    mkdir -p $OUTPUTDIRNAME

    # Directory to look for EM files that were generated
    # Structure is PDBID.ChainLocationRes2
    EMDIRNAME=$(echo $LINE | awk '{ print $2"."$3$4$5 }')

    # Actually call proMute
    $LINE

    # translate the promute command to an SDM script line
    PDBID=$(echo $LINE | awk '{print $2}')
    RESNUM=$(echo $LINE | awk '{print $4}')
    CHAIN=$(echo $LINE | awk '{print $3}')
    FINALRES=$(echo $LINE | awk '{print $5}')
    ORIGRESCODE=$(cat ${PDBID}.pdb | grep ATOM | awk -v rn="$RESNUM" -v chn=$CHAIN '$6 == rn && $5 == chn {print $4}' | sed 's/.*\(...\)/\1/' | head -n 1)
    ORIGRES=${residues[$ORIGRESCODE]}
    # put this line in the script
    echo $(echo $CHAIN $ORIGRES$RESNUM$FINALRES) >> $SDMFILE


    # Move all output files: PDB, Fasta, and EM data
    mv *.txt $OUTPUTDIRNAME
    mv *.pdb $OUTPUTDIRNAME
    mv external/em/$EMDIRNAME $OUTPUTDIRNAME
    echo [INFO] Moving files to $OUTPUTDIRNAME

  done
  # cleanup output
  rm -rf "${FILE}_output"
  mkdir -p "${FILE}_output"
  echo [INFO] Cleaning outputs from last run
  echo [INFO] Moving new files to "${FILE}_output"
  mv $SDMFILE "${FILE}_output/"
  mv *.output "${FILE}_output/"

done
