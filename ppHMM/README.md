# ppHMM to Identify Colicin-like Bacteriocins

[![python2](https://img.shields.io/badge/python-2.7-blue.svg)](https://biopython.org/wiki/Download)

A simple python script to identify two linked pair-profile Hidden Markov Models (ppHMM) from the output of a HMMer search. This pipeline is intended to identify the cytotoxic domain and corresponding immunity protein of large colicin like nuclease bacteriocins. However, this pipeline can be adapted to identify any two pHMMs which are associated with each other.

## Dependencies

ppHMM requires a few python packages which might not be standard and some fairly common bioinformatics software to run

- python2.7
- [Biopython](https://biopython.org/wiki/Download)
- [pygal](http://www.pygal.org/en/stable/)
- [HMMer3](http://hmmer.org)
- [cd-hit](https://github.com/weizhongli/cdhit)
- [Emboss](http://emboss.sourceforge.net/download/)
- [Prodigal](https://github.com/hyattpd/Prodigal/wiki/Introduction)

Additionally, if you want to remove sequences or check what domains are found in your proteins, a local copy of Pfam-A.hmm is required

- [Pfam-A.hmmm](https://pfam.xfam.org)

## Introduction

Nuclease bacteriocins (NBs) are potent antimicrobial proteins which can degrade DNA, tRNA and rRNA. To prevent death of the producing cell, NBs are encoded as a toxin gene followed immediately by a gene encoding an immunity protein. The cytotoxic domains and immunity proteins  contain domains identifiable by pHMMs. The immunity proteins of NBs are always found adjacent to the cytotoxic domain to ensure that an abundance of immunity protein is always in the cell.

## Input
###### Required
- ```-fin --file_in``` The output file of HMMscan
- ```-fout --file_out``` Turns the HMMscan into a '.csv' file for easy reading
- ```-g --genomes``` Directory which stores nucleotide genomes. Each genome should be in its own folder.
- ```-o --orf``` Directory containing the open reading frame predictions for contigs containing the ppHMM hits
- ```-p --profile``` "Text file with names of profiles to be identified formatted as: NAME,1st-profile,2nd-profile,inter-profile-distance"

###### Optional
- ```-m --meta``` Metadata file. Allows ppHMM to find genome files if filenames differ. Also contains species and other information
- ```-r``` Pfam HMM file with binaries. Gives the user the option to remove sequences based on the Pfam domains they contain
- ```-rc``` Pfam HMM file with binaries. Clusters the sequences at 98% sequence identity using cd-hit and removes entire clusters based on Pfam domains

```-fin```

ppHMM Does not perform the HMMscan itself. To identify the two pHMMs, hmmscan should be used on a whole genome translation in all 6 frames ```transeq -sequence genome.fa -outseq 6frame.fa -clean -table 11 -frame 6```. The output of hmmscan should be in the ```--domtblout``` form. If large databases are being scanned hmmscan outputs should be concatenated with a ```#``` separating results.

```-p```

A comma separated file containing the 'NAME' of the profile, the N-pHMM (e.g. NB cytotoxic domain), the C-pHMM (e.g. NB immunity protein), the maximum intergenic region allowed (bp) e.g.

```
HNH,Colicin-DNase,Colicin_Pyocin,60
ColD,Colicin_D,Colicin_immun,60
```

```-r```

If the user is only interested in a certain family of proteins that contain the ppHMM, they have the option to remove sequences which contain different pHMMs. If this option is selected all identified sequences are scanned against Pfam-A. A GUI will appear asking the user to identify Pfam domains that they **DO NOT** want included in the final database.

```-rc```

Open reading frames often make mistakes at the N-terminus. As this region can contain important information about toxin secretion pathways, sequences are clustered at 98% sequence identity and the profiles of the longest sequence are used to include to remove the entire cluster. Gives increased accuracy but can remove some allowed sequences.

## Outputs

The ppHMM will output a summary html containing some basic statistics of the run and a comma separated file of all the hits identified. The html will incorporate a couple of handy plots of the data which are also available as '.svg' files. All of these files will be placed into a directory created when you run the script. The directory will be named after the time of the run. If multiple runs are being performed at the same time then they will be created with an additional suffix do _N until a unique filename is achieved. Output will also include four fasta files of sequences which have passed and failed the profile inspection, one for each profile in the ppHMM pair (if your profiles are in the same ORF it will return the C-terminal ORF).

The output will also include a pickled python variable which contain all of the information for those which passed (```profilePass = True```) or failed the profile inspection (```profilePass = False```). This can be opened using 'cPickle' and allows the user to further analyse the results.
