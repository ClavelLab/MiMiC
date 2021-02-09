# MiMiC
Minimal Microbial Consortia creation tools


# What is MiMiC?

MiMiC proposes minimal microbial consortia from the functional potential of a given metagenomic sample. This is done by producing a binary vector (presence,absence) of all Pfams within the sample and comparing this again our pre-calculated genome database. Through an iterative process, the genomes that best match the metagenomic samples functiona repertoire are selected. Once all functions have been accounted for, or the additiona of further genome accounts for no more Pfams, the process is halted and the user provided with the list and statistics for selection.

We hope that this tool will lead to the further adoption and use of minimal microbial consortia for next generation probiotics as well as basic research into microbial communities.


# Installation

Before installating Protologger, all the dependancies below must be installed;
- [PRODIGAL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- [HMMSCAN](https://academic.oup.com/nar/article/39/suppl_2/W29/2506513)

Additional R modules used by MiMiC that must be installed are;
- dplyr
- plyr
- inflection
- tidyverse
- ggplot2
- taRifx
- readr

Once the denpendancies have been installed, the user may test the installation by running the following commands;

TBA


# Running MiMiC

MiMiC is currently designed as a series of four scripts which the user must apply to a dataset subsequentially (step1...4).

Some of these files requires the user to alter the code to include their files (all files will accept command line arguments in future editions).

<b>Step_1;</b> This script conducts metagenomics assembly, protein prediction and runs hmmscan to identify Pfams. Knowing that the majority of users will have their own assembly pipelines, the only step required is;
`hmmscan --cpu 20 --acc --noali --cut_ga --tblout $RESULTS.table_results Pfam-A_2019.hmm $ASSEMBLY.faa`
In this command, the user must change $RESULTS to their wanted prefix to their results, the Pfam-A_2019.hmm file should link to their prepared version of Pfam that will be used throughout the analysis and the $ASSEMBLY command will link to the predicted proteome of the metagenomic assembly (FAA format).

<b>Step_2;</b> This file parses the output of hmmscan into the format used by MiMiC. Users must modify line 12 and 13 with their directory with their files and their output folders, respectively.

<b>Step_3;</b> This script generates the Pfam vector of the input metegenome using the command; `Rscript script_2_pfam_vector.R Pfam.A.clans $PATH`, where the user must provide the link to the Pfam.A.clans file they are using and the $PATH to their folder with the output from Step_2. 

<b>Step_4;</b> This script conducts the final comparison of the metagenomes Pfam vector file to the genome database. Both such files need to be provided within the R script on line 19 and 122 (genome database file) and line 22 (metagenome Pfam file). The name of the output file can be stated on line 110.


# Example dataset

The 'example_data' folder contains the initial run of MiMiC on the PiBAC collection to generate the published minimal consortia (https://www.nature.com/articles/s41467-020-19929-w).

These files include;
- Pfam_data_base_Pfam-A.clans_30_8_18.tsv, this file is the pfam list used for the analysis.
- PiBC_Genome_Binary_Vector_111_update_3oct_2019.txt, this file has the pfam binary vector of 111 bacterial species (PiBAC).  
- PiBC_Vector_284_updated_3rdOct.txt, pfam vetor for the Metagenome used in PiBAC analysis 


# Citation
We ask that anyone who uses MiMiC cites not only our publication but the list of publications below which provide tools and databases which are integral for MiMiCs working;


### Tools
- [PRODIGAL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- [HMMSCAN](https://academic.oup.com/nar/article/39/suppl_2/W29/2506513)

### Databases 
- [Pfam](https://academic.oup.com/nar/article/47/D1/D427/5144153)
