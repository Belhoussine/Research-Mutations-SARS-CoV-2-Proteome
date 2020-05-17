#!/bin/bash

# After making promute, copy this script into the build directory
# Change the "FILES" variable to the script (or scripts) you would like to execute
# to profile different mutation scripts, you can run this script with the "time" command 

# generate two promute scripts
./makeScript 1hhp A:A 1:30 X emScript
./makeScript 1hhp A:A 1:30 X emsrScript

# append the correct flags
sed -i 's/$/ emsr/' emsrScript
sed -i 's/$/ em/' emScript

FILES=("emScript")

# Execute the lines in the script
for FILE in $FILES; do
  readarray -t LINES < "$FILE"
  for LINE in "${LINES[@]}"; do
    $LINE
  done
  # cleanup output
  mkdir -p "${FILE}_output"
  mv *.txt "${FILE}_output"
  mv *.pdb "${FILE}_output"
done

