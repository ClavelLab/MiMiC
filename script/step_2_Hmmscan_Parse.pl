#!usr/bin/perl 
#Neeraj kumar
#script to parse hmmscan output
#last update 20 mrach 2020

#use : perl step_1_Hmmscan.pl
#need to mention the input and output directory in the variable $srcdir (input directory), $dest(output directory)
use strict;
use warnings;
use File::Copy;

my $srcdir = "/path to input directory";  #source directory
my $dest ="path to output directory";  #target directory

opendir (DIR, $srcdir) or die "can not open $srcdir,$!";

my @files; 
my @line_array;
my @column1;
my ($target_name, $accession_FS, $queryname, $accession,$E_value_FS, $score_FS,$bias_FS,$E_value_BD, $score_BD, $bias_BD, $exp_DNE, $reg_DNE, $clu_DNE, $ov_DNE,
$env_DNE,$dom_DNE,$rep_DNE,$inc_DNE,$description_of_target)=0;

@files = grep {!/^\.+$/ } readdir(DIR);

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
				
				($target_name ,$accession_FS, $queryname,  $accession,  $E_value_FS, $score_FS, $bias_FS,  $E_value_BD, $score_BD, $bias_BD,  $exp_DNE,$reg_DNE, $clu_DNE, $ov_DNE, $env_DNE, $dom_DNE, $rep_DNE, $inc_DNE, $description_of_target)=split('\s+',$line);
				unless ( $E_value_FS >= 0.01 )
					{
						
						print $fh_t "$accession_FS\n";
						print "$accession_FS\t";
						print "$E_value_FS\n";
					}
				#print "@column1\n";
				 
			} 
		}
#		move($old, $dest) or die "Move $old -> $dest failed: $!"; #to move files from one folder to another.
	}  


#script ends here



