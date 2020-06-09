#!/bin/bash

# After making promute, copy this script into the build directory
# Change the "FILES" variable to the script (or scripts) you would like to execute
# to profile different mutation scripts, you can run this script with the "time" command 

# generate list of mutations
# ACE2 Mutations
# ./proMuteBatch 6M0J A:A 31:31 X K31 # exhaustively mutate K at residue 31, chain A
# ./proMuteBatch 6M0J A:A 34:34 X H34 # exhaustively mutate H at residue 34, chain A
# ./proMuteBatch 6M0J A:A 35:35 X E35 # exhaustively mutate E at residue 35, chain A
# RBD Domain mutations
# ./proMuteBatch 6M0J E:E 453:453 X Y453 # exhaustively mutate Y at residue 453, chain E
# ./proMuteBatch 6M0J E:E 492:492 X L492 # exhaustively mutate L at residue 492, chain E
# ./proMuteBatch 6M0J E:E 494:494 X S494 # exhaustively mutate S at residue 494, chain E

# Mutations on second PDB 6Y2E
./proMuteBatch 6Y2E A:A 25:25 X T25 # exhaustively mutate T at residue 25, chain A
./proMuteBatch 6Y2E A:A 26:26 X T26 # exhaustively mutate T at residue 26, chain A
./proMuteBatch 6Y2E A:A 27:27 X L27 # exhaustively mutate L at residue 27, chain A
./proMuteBatch 6Y2E A:A 28:28 X N28 # exhaustively mutate N at residue 28, chain A
./proMuteBatch 6Y2E A:A 29:29 X G29 # exhaustively mutate G at residue 29, chain A
./proMuteBatch 6Y2E A:A 30:30 X L30 # exhaustively mutate L at residue 30, chain A
./proMuteBatch 6Y2E A:A 31:31 X W31 # exhaustively mutate W at residue 31, chain A
./proMuteBatch 6Y2E A:A 32:32 X L32 # exhaustively mutate L at residue 32, chain A


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
# sed -i 's/$/ em/' K31
# sed -i 's/$/ em/' H43
# sed -i 's/$/ em/' E35
# sed -i 's/$/ em/' Y453
# sed -i 's/$/ em/' L492
# sed -i 's/$/ em/' S494

sed -i 's/$/ em/' T25
sed -i 's/$/ em/' T26
sed -i 's/$/ em/' L27
sed -i 's/$/ em/' N28
sed -i 's/$/ em/' G29
sed -i 's/$/ em/' L30
sed -i 's/$/ em/' W31
sed -i 's/$/ em/' L32


# Customize with your scripts
# FILES=("K31" "H34" "E35" "Y453" "L492" "S494")
FILES=("T25" "T26" "L27" "N28" "G29" "L30" "W31" "L32")

# Execute the lines in the script
for FILE in ${FILES[@]}; do


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
    if [ "$ORIGRES" != "$FINALRES" ]
    then
      # put this line in the script
      echo $(echo $CHAIN $ORIGRES$RESNUM$FINALRES) >> $SDMFILE
    fi

    # Move all output files: PDB, Fasta, and EM data
    mv *.txt $OUTPUTDIRNAME
    mv *.pdb $OUTPUTDIRNAME
    mv external/em/$EMDIRNAME $OUTPUTDIRNAME
    mv $FILE $OUTPUTDIRNAME
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
