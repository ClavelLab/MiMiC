![logo](/logo.png)

# MiMiC
Minimal Microbial Consortia creation tools


# What is MiMiC?

MiMiC proposes minimal microbial consortia from the functional potential of a given metagenomic sample. This is done by producing a binary vector (presence,absence) of all Pfams within the sample and comparing this again our pre-calculated genome database. Through an iterative process, the genomes that best match the metagenomic samples functiona repertoire are selected. Once all functions have been accounted for, or the addition of further genome accounts for no more Pfams, the process is halted and the user provided with the list and statistics for selection.

We hope that this tool will lead to the further adoption and use of minimal microbial consortia for next generation probiotics as well as basic research into microbial communities.


# Installation

Before running MiMiC, all the dependancies below must be installed;
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

Once the dependancies have been installed, the user may test the installation by running the following commands;

TBA


# Running MiMiC

MiMiC is currently designed as a series of four scripts which the user must apply to a dataset subsequentially (step1...4).

Some of these files requires the user to alter the code to include their files (all files will accept command line arguments in future editions).

<b>Step_1;</b>  step_1_contig_to_Pfam.sh

Description: Identifies the proteins within your genomes/metagenome and annotate them against the Pfam database.

Usage: `bash step_1_contig_to_Pfam.sh -i PathToInputFolder -d PathToPfamDatabase -t optionforProdigal(meta or single)`

Options:
-i path to input folder where contigs/assemblies are kept 
-d path to hmmpressed pfam database
-t if metagenomic assembly then `meta`, otherwise `single`

Output folder:
output_faa : all faa files are generaed in this folder 
output_pfam : all hmmscan output files are generated in this folder 
output_fna : all fna (prodigal) files are generated in this folder 


<b>Step_2;</b> step_2_hmmscan_parsed.pl

Description: Parses the Hmmscan output into a format MiMiC can utilise in Step 3.

Usage: `perl step_2_hmmscan_parsed.pl path_to_input_folder`

Notes:
The input folder should contain only faa files inside.


<b>Step_3;</b> step_3_pfam_vector.R

Description: This script makes a binary vector for genome/metagenome profile of Pfams.

Usage: `Rscript --vanilla script_3_pfam_vector.R -i pfam.A.clans -p Step2Output -o Step3Output`

Options:
-i Pfam.A.clans file from pfam database 
-p path for input folder where all pfam parsed files are available
-o name of output file name (default is PfamVector.txt)

<b>Step_4;</b> step_4_mimic.R

Description: This script calculates the minimal microbial consortia for a metagenome based the minimum number of species covering maximum number of metagenomic Pfams.

Usage: `Rscript --vanilla step_4_mimic.txt -m metagenomefileName -g genomeVectorFileName -i iterationNumber -o mimicOutputName -k kneepointbasedOutputName`

Options:
-m input metagenome Pfam vector
-g input reference genome Pfam vector
-i iteration number to be ran for each metagenome
-o mimic outputfile name (optional: default MiMiC.txt)
-k mimic output file name with knee point (optional : MiMiC_KneePoint.txt)

# Example dataset

The 'example_data' folder contains the initial run of MiMiC on the PiBAC collection to generate the published minimal consortia (https://www.nature.com/articles/s41467-020-19929-w).

These files include;
- Pfam_data_base_Pfam-A.clans_30_8_18.tsv, this file is the pfam list used for the analysis.
- PiBC_Genome_Binary_Vector_111_update_3oct_2019.txt, this file has the pfam binary vector of 111 bacterial species (PiBAC).  
- PiBC_Vector_284_updated_3rdOct.txt, pfam vetor for the Metagenome used in PiBAC analysis 



# Reference data sets

The reference database file (`MiMiC_reference_Database.tar.gz`) contains four folders inside;

- 1.Host specific, this folder has three host specific species pfam binary vector 
  - human 
  - pig

- 2.Metagenome used in study
  - MetagenomeVector_937.txt : pfam vecor of metagenomes used in the MiMiC paper study
  - All_meta_937.txt : meta data information with sample name and the source of the metagenome

- 3.Mock community 
  - Pfam vector of the mBarc mock community metagenome 

- 4.NCBI RefSeq
  - Pfam vector of the ncbi refseq genomes used in the MiMiC study




# Citation
We ask that anyone who uses MiMiC cites not only our publication but the list of publications below which provide tools and databases which are integral for MiMiCs working;


### Tools
- [PRODIGAL](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119)
- [HMMSCAN](https://academic.oup.com/nar/article/39/suppl_2/W29/2506513)

### Databases 
- [Pfam](https://academic.oup.com/nar/article/47/D1/D427/5144153)
