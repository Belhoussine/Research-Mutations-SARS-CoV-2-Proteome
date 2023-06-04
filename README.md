# Research: Impactful Mutations in Mpro of the SARS-CoV-2 Proteome

### Goal: Single amino acid substitution in COVID-19 spike protein

# Notes

Promute builds fine on `archlinux 5.6.11-arch1-1`

Surface Racer needs [libstdc++5](https://www.archlinux.org/packages/extra/x86_64/libstdc++5/)

### `driver.sh`

This script is meant to compare the time and performance between proMute with Surface Racer enabled and disabled.

Customize the script to run one or both of the generated mutant scripts from `proMuteBatch`. Then *copy* the script into the proMute build directory. Now I know to check the makefile
before writing scripts in the `build` directory...

If you have `time` installed on the system (not the built in shell command, run with `/usr/bin/time`) you can profile the performance of the desired script with

```bash
/usr/bin/time --verbose ./driver.sh
```

### dataParsing.cpp

This program should be run in the em data folder for *every* mutation. Compile with `g++ -o dataParser dataParsing.cpp`

It will generate a `results.csv` file containing the results.
