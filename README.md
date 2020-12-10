# Experiments for MetaGeneMarkS

Georgia Institute of Technology, Atlanta, Georgia, USA

Reference: PAPER LINK

## Overview
This repository contains the data and source code needed to reproduce all results found in the MetaGeneMarkS paper.

## Program Versions
MetaGeneMarKS is a standalone tool, but building the initial set of models relies on GeneMarkS-2 predictions. Similarly, results are compared to multiple external tools, whose versions are shown here:

- GeneMarkS-2:
- (Meta)Prodigal:
- MetaGeneAnnotator:
- FragGeneScan:
- MetaGeneMark:

MetaGeneMarkS is a C++ program. That said, experiments and results are all executed and analyzed in `python`. To get all the packages (for reproducibility), it is recommended that the user creates a `conda` environment from the file in `install/conda_mgm2.yml` through the following command:

    conda env create -f install/conda_mgm2.yml --name mgms

This can then be activate via

    conda activate mgms

See `info/reproduce.[html|pdf]` for more information.

## Reproducing Results

We provide a document detailing how to reproduce all results. This can be found at `info/reproduce.[html|pdf]`

## Folder structure

The following directory structure should be maintained while running experiments

    .
    ├── bin                                   # Executables constructed from python/bash drivers (via install.sh)
    ├── bin_external                          # External tools
    ├── config                                # Configuration files, e.g. MetaGeneMarkS learning parameters
    ├── config.sh                             # Load bash variables for paths to different directories
    ├── install                               # Conda environment file for easy installation
    ├── lists                                 # Lists of genomes (main input method to scripts)
    ├── info                                  # Information about reproducing results
    ├── metadata                              # Non-genomic data, including taxonomy information
    ├── data                                  # Data Location: where all raw data will be stored during runs
    │   ├── GCFID 1                           # ID of genome 1
    │   │   ├── ncbi.gff                      # RefSeq annotation
    │   │   ├── sequence.fasta                # Genomic sequence file
    │   ├── GCFID 2                           # ID of genome 2
    │   │   ├── ncbi.gff                      # RefSeq annotation
    │   │   ├── sequence.fasta                # Genomic sequence file
    │   │   ...
    ├── code                                  # Source code
    │   ├── python                            # Python code
    │   │   ├── driver                        # Drivers that can be executed
    │   │   ├── lib                           # Library files
    |   |── mgms                              # MetaGeneMarkS (C++) source code and Makefile
    │   ├── bash                              # Bash scripts
    │   │   ├── driver                        # Drivers that can be executed
    │   │   ├── lib                           # Library files
    ├── runs                                  # Data Location: where all raw data will be stored during runs
    │   ├── GCFID 1                           # ID of genome 1
    │   │   ├── startlink                     # StartLink runs
    │   │   ├── mgms                          # MGMS runs
    |   |   ├── others...                     # Other tools
    │   ├── GCFID 2                           # ID of genome 2
    │   │   ├── startlink                     # StartLink runs
    │   │   ├── mgms                          # MGMS runs
    |   |   ├── others...                     # Other tools
    │   │   ...

