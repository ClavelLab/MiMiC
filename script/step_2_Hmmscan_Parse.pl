#Author Neeraj Kumar
#Last updated 2021

#Description
# script parsed the hmmscan out put for step 3 

#Usage
# perl step_2_hmmscan_Parsed_v2.pl path_to_input_folder*
# the input folder should contain only faa files inside.


#!usr/bin/perl 
use strict;
use warnings;
use File::Copy;

my $srcdir = $ARGV[0]; # path for source directory where hmmscan output tables are present.

sub main {
    my $directory = "hmmscan_parsed";
    
    unless(mkdir $directory) {
        die "Unable to create $directory\n";
    }
}

main();

my $dest ="hmmscan_parsed/"; #target directory

#my $srcdir = "/media/sparekh/External_storage/mbarc_analysis_2019/mbarc_2019/table"; #source directory
#my $dest ="/media/sparekh/External_storage/mbarc_analysis_2019/mbarc_2019/table_Parsed"; #target directory

opendir (DIR, $srcdir) or die "can not open $srcdir,$!";
#my $file_to_write = 'report.txt';
#open(my $fh_t, '>', $file_to_write) or die "Could not open file '$file_to_write' $!";

my @files; #programe does not work if you do not defined variable a local before using it.(very important)
my @line_array;
my @column1;
my ($target_name, $accession_FS, $queryname, $accession,$E_value_FS, $score_FS,$bias_FS,$E_value_BD, $score_BD, $bias_BD, $exp_DNE, $reg_DNE, $clu_DNE, $ov_DNE,
$env_DNE,$dom_DNE,$rep_DNE,$inc_DNE,$description_of_target)=0;



@files = grep {!/^\.+$/ } readdir(DIR); #grep all filenames in array

foreach my $file (@files)
	{
		my $old = "$srcdir/$file";
		my $new = "$dest/$file";
		open(my $fh,'<:encoding(UTF-8)',$old) or die "could not open file '$old'$! ";
		open(my $fh_t, '>', $new) or die "Could not open file '$new' $!";
		while (my $line = <$fh>) 
		{ 
			unless ($line =~ /#/) 
			{ 
				#print $line;
				#@line_array = split('\s+',$line);#''\s+' split the line based on multiple array.
				#@column1 = ("$line_array[0]\t",$line_array[1]);
				($target_name ,$accession_FS, $queryname,  $accession,  $E_value_FS, $score_FS, $bias_FS,  $E_value_BD, $score_BD, $bias_BD,  $exp_DNE,$reg_DNE, $clu_DNE, $ov_DNE, $env_DNE, $dom_DNE, $rep_DNE, $inc_DNE, $description_of_target)=split('\s+',$line);
				unless ( $E_value_FS >= 0.01 )
					{
						
						print $fh_t "$accession_FS\n";
						#print $fh_t "$queryname\t$accession_FS\t$E_value_FS\n";
						#print $fh_t "$E_value_FS\n";
						print "$accession_FS\t";
						print "$E_value_FS\n";
					}
				#print "@column1\n";
				 
			} 
		}
#		move($old, $dest) or die "Move $old -> $dest failed: $!"; #to move files from one folder to another.
	}  
