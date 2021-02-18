#Metagenome pipeline.
#Neeraj Kumar
#last update Feb 2021
#use : bash metagenome_pipeline.sh
#
#*********************************************#
#tools requirement :

#3.Prodigal
#4.Hmmer
#5.Pfam database

#******* Usage ******************************************

# bash step_1_contig_to_Pfam.sh -h
# bash step_1_contig_to_Pfam.sh -i PathToInputFolder -d PathToPfamDatabase -t optionforProdigal(meta or single)  

#*************** making option for the user to pass the arguments from command line **********

	if [ "$1" == "-h" ]; then
		echo -e "\n\n"
		echo -e "********  Usage ***********************\n\n"
		echo -e "Usage: bash step_1_contig_to_Pfam.sh -i PathToInputFolder -d PathToPfamDatabase -t optionforProdigal(meta or single)\n\n"
		echo -e "[options]\n"		
		echo -e "	-i : this is the path for genome or metagenome assembly folder to be annotated"
		echo -e "	-d : this is the folder for hmmpress pfam database"
		echo -e "	-t : this provides an option to prodigal to be in metagenomic mode or single genome mode mode\n\n"
		exit 0
	fi

	while getopts i:d:t: flag
		do
		    case "${flag}" in
		        i) input_fna=${OPTARG};;
		        d) PfamPath=${OPTARG};;
		        t) Type=${OPTARG};;
		    esac
		done


#****** creating directory for output ***************************

	mkdir output_faa
	mkdir output_pfam
	mkdir output_fna

#******************  path for tools *****************************************
#Prodigal
	ProdigalPath=prodigal
#hmmscan
	hmmscanPath=hmmscan

# running script for each genome or metagenome

	for i in `cat $input_fna/filenames.txt`; # each sample generate two samples R1 and R2. sample names should be the same. only sufix in the file should be different. 
		do
			echo $i


#prodigal --> gene prediction
			$ProdigalPath -i $input_fna$i -p $Type -a output_faa/$i.faa -c -m -q -f sco -d output_fna/$i.fna	
			echo "$i	prodigal finished">>output_faa/Prodigallogfile.txt
			ProteinCount=$(grep -c "^" output_faa/'$i'.faa) #-------------- Counting proteins
		
			echo -e "'$i'\t$ProteinCount" >> output_faa/protein_Count.txt 
			echo "$i	hmmscan started">>output_pfam/hmmscan_logfile.txt

#hmmscan ----> pfam prediction
	
			$hmmscanPath --cpu 20 --acc --noali --cut_ga --tblout  output_pfam/$i.faa.table_results $PfamPath/Pfam-A_2019.hmm   output_faa/$i.faa
			echo "$i	hmmscan finished" >> output_pfam/hmmscan_logfile.txt	

		done

echo -e "processing is finished \n\n" 
#Reference
#Prodigal
#2.Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal:
#prokaryotic gene recognition and translation initiation site identification. BMC 
#Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PubMed PMID:
#20211023; PubMed Central PMCID: PMC2848648.

#Hmmscan

#http://hmmer.org/

#Pfam data base




