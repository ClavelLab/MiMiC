#Metagenome pipeline.
#Neeraj Kumar
#last update 20 March 2020
#use : bash metagenome_pipeline.sh
#
#*********************************************#
#tools requirement :

#1.Spade 
#2.Seqkit
#3.Prodigal
#4.Hmmer
#5.Pfam database

#****** path for tools ***************************

	input="" #----- path to files for metagenome
	path_to_processed_file="" # path to processed files


#******************  path for tools *****************************************
#fastp
#fastpPath=/fastp

#SpadesPath
#	metaspade=spades.py

#seqkit
#	seqkitPath=seqkit
#Prodigal
#	ProdigalPath=prodigal
#hmmscan
#	hmmscanPath=hmmscan
#Pfam database
 PfamPath="pfam_2019"

# files to be annotated  

for i in `cat test_sample.txt`; # each sample generate two samples R1 and R2. sample names should be the same. only sufix in the file should be different. 
do
r1=$i."rmhost.pair.1.fq.gz" # change here as per your file suffix
r2=$i."rmhost.pair.2.fq.gz"


#********* running tools step by step *************************************************

#fastp -> quality control

	echo "$i	fastp running">>fastplogfile.txt
	$fastpPath -i $input/$r1 -I $input/$r2 -o $path_to_processed_file/$i.1.fastq_QC.fastq.gz -O $path_to_processed_file/$i.2.fastq_QC.fastq.gz -3 20 -l 35 -q 15 -h $i.html -j $i.json -w 16 
	echo "$i	fastp completed" >>fastplogfile.txt

#assembly -> raw reads to contigs

	echo "$i	assembly started">>$path_to_processed_file/ASMlogfile.txt
	#spades.py -1 $input/$r1 -2 $input/$r2 --meta -o $i.spade 
	$metaspade -1 $path_to_processed_file/$i.1.fastq_QC.fastq.gz -2 $path_to_processed_file/$i.2.fastq_QC.fastq.gz --meta -o $path_to_processed_file/$i.spade
	echo "$i	assembly finished" >> $path_to_processed_file/ASMlogfile.txt
	ContigCountAssembly=$(grep -c "^>" $path_to_processed_file/$i.spade/contigs.fasta) # --------> counting contings
	#rm $input/$r1
	#rm $input/$r2
#seqkit --> contig length cutoff

	$seqkitPath seq -m 500 $path_to_processed_file/$i.spade/contigs.fasta >> $path_to_processed_file/$i.seqkit_contig.fasta
	echo "$i	prodigal statred">>$path_to_processed_file/Prodigallogfile.txt


	ContigCountSeqKit_500=$(grep -c "^>" $path_to_processed_file/'$i'.seqkit_contig.fasta) # ----> contigCount after cutoff

#prodigal --> gene prediction

	$ProdigalPath -i $path_to_processed_file/$i.seqkit_contig.fasta -p meta -a $path_to_processed_file/$i.seqkit_contigs_prod.faa -c -m -q -f sco -d $path_to_processed_file/$i.fna	
	echo "$i	prodigal finished">>$path_to_processed_file/Prodigallogfile.txt
	ProteinCount=$(grep -c "^>" $path_to_processed_file/'$i'.seqkit_contigs_prod.faa) #-------------- Counting proteins
	
	echo -e "'$i'\t$ContigCountAssembly\t$ContigCountSeqKit_500\t$ProteinCount" >> $path_to_processed_file/Contig_protein_Count.txt 
	echo "$i	hmmscan started">>$path_to_processed_file/hmmscan_logfile.txt

#hmmscan ----> pfam prediction

	$hmmscanPath --cpu 20 --acc --noali --cut_ga --tblout  $path_to_processed_file/$i.faa.table_results $PfamPath/Pfam-A_2019.hmm   $path_to_processed_file/$i.seqkit_contigs_prod.faa
	echo "$i	hmmscan finished" >> $path_to_processed_file/hmmscan_logfile.txt	

done



#Reference
#1: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ
#preprocessor. Bioinformatics. 2018 Sep 1;34(17):i884-i890. doi:
#10.1093/bioinformatics/bty560. PubMed PMID: 30423086; PubMed Central PMCID:
#PMC6129281.


#2.Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal:
#prokaryotic gene recognition and translation initiation site identification. BMC 
#Bioinformatics. 2010 Mar 8;11:119. doi: 10.1186/1471-2105-11-119. PubMed PMID:
#20211023; PubMed Central PMCID: PMC2848648.

#Bankevich A, Nurk S, Antipov D, Gurevich AA, Dvorkin M, Kulikov AS, Lesin VM, 
#Nikolenko SI, Pham S, Prjibelski AD, Pyshkin AV, Sirotkin AV, Vyahhi N, Tesler G,
#Alekseyev MA, Pevzner PA. SPAdes: a new genome assembly algorithm and its
#applications to single-cell sequencing. J Comput Biol. 2012 May;19(5):455-77.
#doi: 10.1089/cmb.2012.0021. Epub 2012 Apr 16. PubMed PMID: 22506599; PubMed
#Central PMCID: PMC3342519.


#http://hmmer.org/




