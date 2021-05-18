#!/usr/bin/perl

#####################################
#                                   #
#     Written by Kai Battenberg     #
#     Plant Sciences UC Davis       #
#                                   #
#####################################

use strict;
use warnings;

use Getopt::Long;
use Cwd;
use lib cwd;
use File::Path;

#####SCRIPT DESCRIPTION#####
#Script "step-05-06.pl" finds orthologs in the DATABASE file for each sequence in the QUERY file.
#All output files will be in the working directory.
##########



#####Options#####
#setting the default values.
my $help = "";

my $procedures = "01 02 03 04 05 06 07 08 09"; #e.g. 01 02 03...09.

my $logprint = "YES"; #"YES" to output a progress log or "NO" to not keep it.
my $silent = "NO"; #"NO" to output progress into stdout or "YES" to keep the stdout silent.
my $wrap = "NO"; #"NO" when this script is ran as is or "YES" when it is parallelized.

my $spp_list = "";
my $og = ""; #Full names of outgroup species.

my $query = "";
my $q_seq_type = ""; #"DNA" or "AA".
my $database = ""; #path to the directory with the database files included
my $db_seq_type = ""; #"AA" or "DNA and AA".

my $blast_type = "SWIPE"; #type of BLAST. "AB", "NCBI" or "SWIPE".
my $w_threshold = 3; #minimum word length required for a hit(W in AB-BLAST).
my $eval_threshold = "1e-3"; #minimum e-value.
my $identity_threshold = 0; #minimum %identity.
my $length_threshold = 0; #minimum alignment length.
my $sander_schneider = "YES"; #"YES" or "NO", enables threshold described in Sander and Schneider (1991).
my $loci_threshold = 5; #maximum number of loci to keep per species.

my $i = 1.2; #inflation rate in MCL (between 1.2-6), the higher the I the more things are removed.

my $rooting = "MI"; #Option for rooting the tree. "MI" for maximum ingroup method (RT method) described in Yang and Smith (2014) or "MD" for most distant outgroup.

my $vraxml = "AVX"; #the type of RAxML compatible with the environment. "AVX" or "SSE3".
my $branch_threshold = 2; #the longest branch length allowed within the last ML tree.

my $overlap_threshold = 0; #the number of overlapping species tolerated by orthologID to be called a speciation event.

my $threads = 2; #number of threads to use for the run.

my $conversion_list = ""; #This variable is optional.

#making the options into external arguments.
GetOptions (
	'help' => \$help,
	
	'procedures=s' => \$procedures,
	
	'log_print=s' => \$logprint,
	'silent=s' => \$silent,
	'wrap=s' => \$wrap,
	
	'spp_list=s' => \$spp_list,
	'og=s@' => \$og,
	
	'query=s' => \$query,
	'q_seq_type=s' => \$q_seq_type,
	'database=s' => \$database,
	'db_type=s' => \$db_seq_type,
	
	'blast_type=s' => \$blast_type,
	'w_threshold=i' => \$w_threshold,
	'eval_threshold=s' => \$eval_threshold,
	'identity_threshold=i' => \$identity_threshold,
	'length_threshold=i' => \$length_threshold,
	'sander_schneider=s' => \$sander_schneider,
	'loci_threshold=i' => \$loci_threshold,
	
	'i=s' => \$i,
	
	'rooting=s' => \$rooting,
	
	'vraxml=s' => \$vraxml,
	'branch_threshold=s' => \$branch_threshold,
	
	'overlap_threshold=i' => \$overlap_threshold,
	
	'threads=i' => \$threads,
	
	'conversion_list=s' => \$conversion_list
);

#printing help
if ($help) {
	print "\n";
	print "\tstep-05-06.pl, for each sequence in the query, searches the database for orthologs.\n\n";
	print "\t--help\tPrint this massage and die.\n\n";
	print "\t--query\tPath to a single file with all the formatted query sequences.\n\n";
	print "\t--q_seq_type\tSequence type of the query. Values 'DNA' or 'AA' are accepted.\n\n";
	print "\t--database\tFile path to the directory with the database files.\n\n";
	print "\t--db_type\tSequence type of database. 'AA' when DNA sequences are not available, or 'DNA and AA' when DNA sequences are avialable along with AA sequences.\n\n";
	print "\t--blast_type\tType of BLAST that will be run. 'AB' for running AB-BLAST, 'NCBI' for running NCBI-BLAST, 'SWIPE' for running SWIPE. (DEFAULT: SWIPE)\n\n";
	print "\t--vraxml\tType of RAxML installed. 'AVX' or 'SSE3'. (DEFAULT: AVX)\n\n";
	print "\t--spp_list\tPath to the file with the lsit of species. See spp_list_EXAMPLE.txt for the proper format.\n\n";
	print "\t--og\tA list of one or more full species name within the species list selected as the outgroup.\n\n";
	print "\t--conversion_list\tPath to a single file with the gene name converison. (optional)\n\n";
	print "\t--w_threshold\tMinimum word count for AB-BLAST to establish a seed. Only meaningful when AB-BLAST is used. (DEFAULT: 3)\n\n";
	print "\t--eval_threshold\tMinimum e-value for a hit to be kept. (DEFAULT: 1e-3)\n\n";
	print "\t--identity_threshold\tMinimum \% identity for a hit to be kept. (DEFAULT: 0)\n\n";
	print "\t--length_threshold\tMinimum alignment length for a hit to be kept. (DEFAULT: 0)\n\n";
	print "\t--sander_schneider\tUse the function provided in Sander \& Schneider (1991) to screen hits. 'YES' to implement, or 'NO' to not implement. (DEFAULT: YES)\n\n";
	print "\t--loci_threshold\tNumber of loci kept from hits per species. (DEFAULT: 5)\n\n";
	print "\t--i\tInflation rate for MCL. Range from 1.2 (most relaxed) to 6.0 (most stringent). (DEFAULT: 1.2)\n\n";
	print "\t--rooting\tOpriont for rooting the tree. 'MI' for outgroup that is most inclusive, 'MD' for outgroup that is most distant. (DEFAULT: MI)\n\n";
	print "\t--branch_threshold\tMaximum branch length allowed to be kept in the final tree for orthology prediction. (DEFAULT: 2)\n\n";
	print "\t--overlap_threshold\tMaximum number of overlapping species across child clades tolerated by orthologID for a node be called a speciation event.\n\n";
	print "\t--threads\tNumber of threads used. (DEFAULT: 2)\n\n";
	print "\t--procedures\tProcedures (subsections) that are run. '01 02 03 04 05 06 07 08 09' to run the complete process. (DEFAULT: 01 02 03 04 05 06 07 08 09)\n\n";
	print "\t--log_print\tPrint the progress into a log file for the procedures conducted. 'YES' to print out or 'NO' to not print out. (DEFAULT: YES)\n\n";
	print "\t--silent\tPrint the progress into the standard output. 'YES' to print out or 'NO' to not print out. (DEFAULT: YES)\n\n";
	print "\t--wrap\tProgram run using a SLURM wrapper. Use 'YES' when it is run by a wrapper, and use 'NO' when it is not run by a wrapper. (DEFAULT: NO).\n";
	die "\n";
}

#checking for required options.
if (!$query) {
	die "USAGE: option --query is required.\n";
}
if (!$q_seq_type) {
	die "USAGE: option --q_seq_type is required.\n";
}
if (!$database) {
	die "USAGE: option --database is required.\n";
}
if (!$db_seq_type) {
	die "USAGE: option --db_type is required.\n";
}
if (!$spp_list) {
	die "USAGE: option --spp_list is required.\n";
}
if (!$og) {
	die "USAGE: option --og is required.\n";
}
##########



#####(Procedure-00)Setting up the environment.#####
###Identifying the working directory.
my $pwd = cwd();
my $pwd_ortho;
my @ls;
my $ls;
my @ls_ortho;

opendir my $dh, $pwd or die "Couldn't open dir '$pwd': $!";
@ls = readdir $dh;
closedir $dh;
$ls = join ("\t", @ls);
###

###Checking all the required files###
#"OrthoReD_vXXXXXXXX" Directory with all scripts needed.
my $test_st = 1;
foreach my $file (@ls) {
	if ($file =~ m/OrthoReD/) {
		$test_st = 0;
		
		$pwd_ortho = $pwd."/".$file;
		opendir my $dh, $pwd_ortho or die "Couldn't open dir '$pwd_ortho': $!";
		@ls_ortho = readdir $dh;
		closedir $dh;
		
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-05-06.pl> requires folder <OrthoReD_vXXXXXXXX> to be in the working directory.\n";
}
else {
	$test_st = 1;
}

#"step-05-06.pl" script for this process.
foreach my $file (@ls_ortho) {
	if ($file =~ m/step-05-06.pl/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-05-06.pl> needs to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}

#"OrthoReD_library.pm" Perl module with all the functions.
foreach my $file (@ls_ortho) {
	if ($file =~ m/OrthoReD_library.pm/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-05-06.pl> requires file <OrthoReD_library.pm> to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}

#"header_cleaner.pm" Perl module with header_cleaner function.
foreach my $file (@ls_ortho) {
	if ($file =~ m/header_cleaner.pm/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-05-06.pl> requires file <OrthoReD_library.pm> to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}

#BLAST
if ($blast_type =~ m/AB/) {
	my $blast_check = `which xdformat`;
	if (!$blast_check) {
		die "USAGE: <step-05-06.pl> requires AB-BLAST to be running unless --blast_type is set to NCBI or SWIPE.\n";
	}
}
elsif ($blast_type =~ m/NCBI/) {
	my $blast_check = `which segmasker`;
	if (!$blast_check) {
		die "USAGE: <step-05-06.pl> requires NCBI-BLAST to be running unless --blast_type is set to AB.\n";
	}
}
elsif ($blast_type =~ m/SWIPE/) {
	my $blast_check = `which swipe`;
	if (!$blast_check) {
		die "USAGE: <step-05-06.pl> requires SWIPE to be running unless --blast_type is set to AB or NCBI.\n";
	}
	$blast_check = `which blastdbcmd`;
	if (!$blast_check) {
		die "USAGE: <step-05-06.pl> requires NCBI-BLAST to be running unless --blast_type is set to AB.\n";
	}
}
else {
	die "USAGE: <step-05-06.pl> requires --blast_type to be AB or NCBI.\n";
}

#MCL
my $mcl_check = `which mcl`;
if (!$mcl_check) {
	die "USAGE: <step-05-06.pl> requires MCL to be running.\n";
}

#MAFFT
my $mafft_check = `which mafft`;
if (!$mafft_check) {
	die "USAGE: <step-05-06.pl> requires MAFFT to be running.\n";
}

#RAxML
my $raxml_check;
if ($vraxml =~ m/^AVX$/) {
	$raxml_check = `which raxmlHPC-PTHREADS-AVX`;
	if (!$raxml_check) {
		die "USAGE: <step-05-06.pl> requires raxmlHPC-PTHREADS-AVX to be running.\n";
	}
}
elsif ($vraxml =~ m/^SSE3$/) {
	$raxml_check = `which raxmlHPC-PTHREADS-SSE3`;
	if (!$raxml_check) {
		die "USAGE: <step-05-06.pl> requires raxmlHPC-PTHREADS-SSE3 to be running.\n";
	}
}
else {
	die "USAGE: <step-05-06.pl> requires --vraxml to be AVX or SSE3.\n";
}

#Newick utilities
my $newick_check = `which nw_reroot`;
if (!$newick_check) {
	die "USAGE: <step-05-06.pl> requires Newick utilities to be running.\n"
}
###

###Setting up all the building blocks needed
#setting up to use the perl libraries
use OrthoReD_v20180830::OrthoReD_library;
use OrthoReD_v20180830::header_cleaner;

#generating a uniq id for the run
my $id = int(rand(1000000000000));

#Making the subject directory
my $input_dir = "Step-05_INPUT";
mkdir $input_dir unless -d $input_dir;

#Making the product directory
my $product_dir = "Step-06_ORTHOLOG";
mkdir $product_dir unless -d $product_dir;

#Identify all the species included in the the species list.
my $full_name = 0;
my $short_name = 1;
my $acronym_name = 2;

my %spp1 = %{ &OrthoReD_library::species_lister ($spp_list, $acronym_name, $full_name) };
my %spp2 = %{ &OrthoReD_library::species_lister ($spp_list, $short_name, $full_name) };

#Checking for wrap status
if ($wrap !~ m/^YES$/ && $wrap !~ m/^NO$/) {
	die "USAGE: only YES or NO is allowed for option --wrap.\n";
}

#Identify all the outgroup species in the species list.
my @ogs = @{ $og };
if ($wrap =~ m/YES/) {
	foreach my $og (@ogs) {
		$og =~ s/_/ /;
	}
}

#Identify the QUERY file and copy it into directory Step-05_INPUT.
my @file_name;
my $new_query;
@file_name = split (/\//, $query);
my $temp_query = $pwd."/".$input_dir."/".$file_name[-1];
if ($q_seq_type =~ m/AA/) {
	$new_query = $temp_query;
	system "cp -n $query $new_query 2>/dev/null";
}
if ($q_seq_type =~ m/DNA/) {
	$new_query = $temp_query;
	$new_query =~ s/.fas/_AA.fas/;
	&OrthoReD_library::translator2($query, $new_query);
}

#Identify all the queries
my %querys = %{ &OrthoReD_library::fasta2hash($new_query) };

#Identify the database
$database =~ s/\/$//;
my $new_database = $pwd."/".$input_dir."/DATABASE";
if ($wrap !~ m/YES/) {
	system "rsync -a $database/* $new_database";
}

#opendir $dh, $new_database or die "Couldn't open dir '$new_database': $!";
my @new_db;
#@new_db = readdir $dh;
#closedir $dh;
#shift @new_db;
#shift @new_db;
my $new_db = `ls $new_database`;
@new_db = split (/\n/, $new_db);

print "\n\nIs the database actually there in $new_database?\n";
foreach my $db (@new_db) {
	print "\t$db\n";
}
die "\n\n\n";

my $database_dna_file;
my $database_aa_file;
if ($db_seq_type =~ m/DNA/) {
	foreach my $file (@new_db) {
		if ($file =~ m/DNA/ and $file =~ m/fas$/) {
			$database_dna_file = $new_database."/".$file;
		}
	}
	if (!$database_dna_file) {
		die "USAGE: DNA database could not be identified.\n";
	}
}
foreach my $file (@new_db) {
	if ($file =~ m/AA/ and $file =~ m/fas$/) {
		$database_aa_file = $new_database."/".$file;
	}
}
if (!$database_aa_file) {
	die "USAGE: AA database could not be identified.\n";
}

my $database_dna;
my $database_aa;
my $database_temp;
if ($blast_type =~ m/AB/) {
	$database_temp = $new_database."/AB";
}
elsif ($blast_type =~ m/NCBI/ || $blast_type =~ m/SWIPE/) {
	$database_temp = $new_database."/NCBI";
}
if ($db_seq_type =~ m/DNA/) {
	$database_dna = `ls $database_temp | grep DNA`;
	$database_dna =~ s/\n$//;
	$database_dna =~ s/\r$//;
	if (!$database_dna) {
		die "USAGE: DNA database could not be identified.\n";
	}
	$database_dna = $database_temp."/".$database_dna;
}
$database_aa = `ls $database_temp | grep AA`;
$database_aa =~ s/\n$//;
$database_aa =~ s/\r$//;
if (!$database_aa) {
	die "USAGE: AA database could not be identified.\n";
}
$database_aa = $database_temp."/".$database_aa;

#making a hash for the entire database.
my %database_dna;
my %database_aa;
if ($db_seq_type =~ m/DNA/) {
	(%database_dna) = %{ &OrthoReD_library::fasta2hash($database_dna_file) };
}
(%database_aa) = %{ &OrthoReD_library::fasta2hash($database_aa_file) };

#adjust procedures and og for wrap
if ($wrap =~ m/YES/) {
	$procedures =~ s/_/ /g;
	foreach my $ogs (sort @ogs) {
		$ogs =~ s/_/ /g;
	}
}
###

###Generating a log.
if ($wrap =~ m/YES/) {
	$logprint = "NO";
}

if ($logprint !~ m/^YES$/ && $logprint !~ m/^NO$/) {
	die "USAGE: only YES or NO is allowed for option --log_print.\n";
}

my $log = "step-05-06_LOG.txt";
if ($logprint =~ m/^YES$/) {
	#Making a log file
	open (LOG, ">$log") or die "cannot open $log.\n";
	
	#printing title
	print LOG "RUNNING SCRIPT: $pwd_ortho step-05-06.pl\n";
	print LOG "(Procedure-00)Setting up the environment.\n";
	
	#printing the command called
	print LOG "options were called as bellow:\n";
	print LOG "\t--query $query\n";
	print LOG "\t--q_seq_type $q_seq_type\n";
	print LOG "\t--database $database\n";
	print LOG "\t--db_seq_type $db_seq_type\n";
	print LOG "\t--blast_type $blast_type\n";
	print LOG "\t--spp_list $spp_list\n";
	print LOG "\t--og @ogs\n";
	print LOG "\t--w_threshold $w_threshold\n";
	print LOG "\t--eval_threshold $eval_threshold\n";
	print LOG "\t--identity_threshold $identity_threshold\n";
	print LOG "\t--length_threshold $length_threshold\n";
	print LOG "\t--sander_schneider $sander_schneider\n";
	print LOG "\t--loci_threshold $loci_threshold\n";
	print LOG "\t--vraxml $vraxml\n";
	print LOG "\t--i $i\n";
	print LOG "\t--rooting $rooting\n";
	print LOG "\t--branch_threshold $branch_threshold\n";
	print LOG "\t--overlap_threshold $overlap_threshold\n";
	print LOG "\t--threads $threads\n";
	print LOG "\t--conversion_list $conversion_list\n";
	
	#printing species
	print LOG "Following species were included in the species list:\n";
	foreach my $sp (sort keys %spp1) {
		print LOG "\t$spp1{$sp}\n";
	}
	
	#printing the selected outgroups
	print LOG "Following species were used as outgroup:\n";
	foreach my $ogsp (@ogs) {
		print LOG "\t$ogsp\n";
	}
	
	#printing the QUERY file
	print LOG "Sequences in the following files were used as queries:\n";
	print LOG "\t$query\n";
	
	#printing the DATABASE file
	print LOG "Contents of the following folders are selected as the database:\n";
	print LOG "\t$database\n";
	
	#Identifying the called procedures
	print LOG "Following procedures were conducted:\n";
	print LOG "\t$procedures\n";
	
	#description for each called procedure
	if ($procedures =~ m/01/) {
		print LOG "(Procedure-01)For each query sequence, BLAST was run against the masked DATABASE to generate a tabulated hit file.\n";
		print LOG "BLASTP parameters were optimized according to Moreno-Hagelsieb and Latimer (2008).\n";
		print LOG "outputformat:\n";
		print LOG "\tquery_seq_id\tsubject_seq_id\te-value\talignment_score\talignment_length\t\%identity\n";
	}
	
	if ($procedures =~ m/02/) {
		print LOG "(Procedure-02)For each HIT file, database file with only unique top hits was generated.\n";
		print LOG "For subjects with multiple hits to the query, all the hits but the one wiht the best e-value were removed.\n";
		print LOG "Top $loci_threshold loci were kept.\n";
		print LOG "Only the isoform with the best e-value was kept.\n";
	}
	if ($procedures =~ m/03/) {
		print LOG "(Procedure-03)For each HITS_TOP database, All_v_All BLASTP hit file was generated.\n";
		print LOG "For subjects with multiple hits to the query, all the hits but the one wiht the best e-value were removed.\n";
	}
	if ($procedures =~ m/04/) {
		print LOG "(Procedure-04)For each HITS_TOP_BLAST file, a FASTA file with all the members of the MCL cluster with the query was generated.\n";
		print LOG "Hits that did not match the following thresholds were removed:\n";
		print LOG "\te-value <= $eval_threshold\n\t\%identity >= $identity_threshold\n\talignment length >= $length_threshold\n";
		if ($sander_schneider =~ m/YES/) {
			print LOG "\tDynamic threshold was implemented for \%identity and alignment length, according to Sander and Schneider (1991).\n";
		}
		print LOG "The MCL inflation rate was set to $i\n";
	}
	if ($procedures =~ m/05/) {
		print LOG "(Step-05)For each HITS_TOP_MCL file, an aligned FASTA file was generated.\n";
		print LOG "MAFFT was used to align the sequences.\n";
		print LOG "Stop codons\(\*\) were not included in the alignement.\n";
	}
	if ($procedures =~ m/06/) {
		print LOG "(Step-06)For each MCL_AA_ALN file, a rooted maximum likelihood tree was generated.\n";
		print LOG "RAxML was used to construct the trees.\n";
		if ($rooting =~ m/MD/) {
			print LOG "The most distant outgroup member from the gene of interest was selected as the root.\n";
		}
		elsif ($rooting =~ m/MI/) {
			print LOG "The most distant outgroup that would result in the largest number of ingroup species was selected as the root.\n";
		}
		print LOG "When no tips of outgroup is present in the tree, the tree was rooted at midpoint.\n";
	}
	if ($procedures =~ m/07/) {
		print LOG "(Step-07)For each MCL_AA_TRE file, a single-lined FASTA file with orthologs was generated.\n";
		print LOG "Before the selection of orthologs, all branches longer than $branch_threshold substitutions/site were cut from MCL_AA_TRE file.\n";
		print LOG "Any splcie variants present in the database were also included.\n";
	}
	if ($procedures =~ m/08/) {
		print LOG "(Procedure-08)Summarizing the output.\n";
	}
	if ($procedures =~ m/09/) {
		print LOG "(Procedure-09)Formatting directories.\n";
	}
}
###

###printing out procedures
if ($wrap =~ m/YES/) {
	$silent = "YES";
}

if ($silent !~ m/^YES$/ && $silent !~ m/^NO$/) {
	die "USAGE: only YES or NO is allowed for option --silent.\n";
}

if ($silent !~ m/^YES$/) {
	#print title
	print "RUNNING SCRIPT: $pwd_ortho step-05-06.pl\n";
	print "(Procedure-00)Setting up the environment.\n";
	
	#printing the command called
	print "options were called as bellow:\n";
	print "\t--query $query\n";
	print "\t--q_seq_type $q_seq_type\n";
	print "\t--database $database\n";
	print "\t--db_seq_type $db_seq_type\n";
	print "\t--blast_type $blast_type\n";
	print "\t--spp_list $spp_list\n";
	print "\t--og @ogs\n";
	print "\t--w_threshold $w_threshold\n";
	print "\t--eval_threshold $eval_threshold\n";
	print "\t--identity_threshold $identity_threshold\n";
	print "\t--length_threshold $length_threshold\n";
	print "\t--sander_schneider $sander_schneider\n";
	print "\t--loci_threshold $loci_threshold\n";
	print "\t--vraxml $vraxml\n";
	print "\t--i $i\n";
	print "\t--rooting $rooting\n";
	print "\t--branch_threshold $branch_threshold\n";
	print "\t--overlap_threshold $overlap_threshold\n";
	print "\t--threads $threads\n";
	print "\t--conversion_list $conversion_list\n";
	
	#printing species
	print "Following species are included in the species list:\n";
	foreach my $sp (sort keys %spp1) {
		print "\t$spp1{$sp}\n";
	}
	
	#printing the selected outgroups.
	print "Following species are used as outgroup:\n";
	foreach my $ogsp (@ogs) {
		print "\t$ogsp\n";
	}
	
	#printing the QUERY file
	print "Sequences in the following files were used as queries:\n";
	print "\t$query\n";
	
	#printing the DATABASE file
	print "Contents of the following folders are selected as the database:\n";
	print "\t$database\n";
	
	#printing the called procedures
	print "Following procedures were conducted:\n";
	print "\t$procedures\n";
	
	#description for each called procedure
	if ($procedures =~ m/01/) {
		print "(Procedure-01)For each query sequence, BLAST is run against the masked DATABASE to generate a tabulated hit file.\n";
		print "BLASTP parameters were optimized according to Moreno-Hagelsieb and Latimer (2008).\n";
		print "outputformat:\n";
		print "\tquery_seq_id\tsubject_seq_id\te-value\talignment_score\talignment_length\t\%identity\n";
	}
	if ($procedures =~ m/02/) {
		print "(Procedure-02)For each BLAST_HIT file, database file with only unique top hits is generated.\n";
		print "For subjects with multiple hits to the query, all the hits but the one with the best e-value are removed.\n";
		print "Top $loci_threshold loci are kept.\n";
		print "Only the isoform with the best e-value is kept.\n";
	}
	if ($procedures =~ m/03/) {
		print "(Procedure-03)For each TOP_BLAST database file, All_v_All BLASTP hit file is generated.\n";
		print "For subjects with multiple hits to the query, all the hits but the one with the best e-value are removed.\n";
	}
	if ($procedures =~ m/04/) {
		print "(Procedure-04)For each TOP_BLAST_HITS file, a FASTA file with all the members of the MCL cluster with the query is generated.\n";
		print "Hits that did not match the following thresholds were removed:\n";
		print "\te-value <= $eval_threshold\n\t\%identity >= $identity_threshold\n\talignment length >= $length_threshold\n";
		if ($sander_schneider =~ m/YES/) {
			print "\tDynamic threshold was implemented for \%identity and alignement length, according to Sander and Schneider (1991).\n";
		}
		print "The MCL inflation rate was set to $i\n";
	}
	if ($procedures =~ m/05/) {
		print "(Step-05)For each HITS_TOP_MCL file, an aligned FASTA file is generated.\n";
		print "MAFFT is used to align the sequences.\n";
		print "Stop codons\(\*\) were not included in the alignement.\n";
	}
	if ($procedures =~ m/06/) {
		print "(Step-06)For each MCL_AA_ALN file, a rooted maximum likelihood tree is generated.\n";
		print "RAxML was used to construct the trees.\n";
		if ($rooting =~ m/MD/) {
			print "The most distant outgroup member from the gene of interest was selected as the root.\n";
		}
		elsif ($rooting =~ m/MI/) {
			print "The most distant outgroup that would result in the largest number of ingroup species was selected as the root.\n";
		}
		print "When no tips of outgroup is present in the tree, the tree was rooted at midpoint.\n";
	}
	if ($procedures =~ m/07/) {
		print "(Step-07)For each MCL_AA_TRE file, a single-lined FASTA file with orthologs is generated.\n";
		print "Before the selection of orthologs, all branches longer than $branch_threshold substitutions/site are cut from MCL_AA_TRE file.\n";
		print "Any splcie variants present in the database were also included.\n";
	}
	if ($procedures =~ m/08/) {
		print "(Procedure-08)Summarizing the output.\n";
	}
	if ($procedures =~ m/09/) {
		print "(Procedure-09)Formatting directories.\n";
	}
}
###
##########



#####(Procedure-01)For each query sequence, file with all BLASTP hits against the DATABASE is generated.#####
#Setting file names
my $procedure_01 = "01";
my $suffix_01 = "_HITS";
my $procedure_01_dir = $procedure_01.$suffix_01;

if ($procedures =~ m/01/) {
	mkdir $procedure_01_dir unless -d $procedure_01_dir;
	my $id_01 = "Procedure01_".$id;
	
	#adjusting the file name for BLASTP.
	my @database_aa = split (/\//, $database_aa);
	my $database_in = $database_aa."/".$database_aa[-1];
	
	my $blastp_out;
	if ($blast_type =~ m/AB/) {
		$blastp_out = &OrthoReD_library::blaster_AB ($new_query, $database_in, $w_threshold, $id_01);
	}
	elsif ($blast_type =~ m/NCBI/) {
		$blastp_out = &OrthoReD_library::blaster_NCBI ($new_query, $database_in, $id_01, $threads);
	}
	elsif ($blast_type =~ m/SWIPE/) {
		$blastp_out = &OrthoReD_library::blaster_SWIPE ($new_query, $database_in, $id_01, $threads);
	}
	
	my $inquerys = `grep '>' $new_query`;
	$inquerys =~ s/\r$//;
	$inquerys =~ s/\n$//;
	$inquerys =~ s/>//g;
	my @inquerys = split (/\n/, $inquerys);
	
	my $procedure01_temp_dir = "Procedure01_".$id;
	mkdir ($procedure01_temp_dir);
	my (@outquerys) = @{ &OrthoReD_library::query_separator ($blastp_out, $procedure01_temp_dir) };
	
	my %outquerys;
	foreach my $outquery (@outquerys) {
		$outquerys{$outquery} = "PRESENT";
	}
	
	foreach my $inquery (@inquerys) {
		my $outfile = $procedure_01_dir."/".$inquery.$suffix_01.".txt";
		my $fromfile = $procedure01_temp_dir."/".$inquery;
		#handling files that did not have a hit
		if (!exists $outquerys{$inquery}) {
			system "touch $outfile";
		}
		else {
			system "mv $fromfile $outfile";
		}
	}
	system "rm -rf $procedure01_temp_dir";
	unlink ($blastp_out);
}
##########



#####(Procedure-02)For each HIT file, a database with only unique top hits is generated.#####
#Setting file names
my $procedure_02 = "02";
my $suffix_02 = "_HITS_TOP";
my $procedure_02_dir = $procedure_02.$suffix_02;

if ($procedures =~ m/02/) { if ($ls !~ m/$procedure_01_dir/ and $procedures !~ m/01/) { die "USAGE: Procedure-02 requires directory $procedure_01_dir.\n"};
	#making required folders
	mkdir $procedure_02_dir;
	
	#Selecting files
	my @hit_files;
	foreach my $query (sort keys %querys) {
		my $file_name = $query.$suffix_01.".txt";
		push (@hit_files, $file_name);
	}
	
	foreach my $element (@hit_files) {
		#query name
		my $query = $element;
		$query =~ s/.txt$//;
		$query =~ s/$suffix_01$//;
		
		#hit file
		my $hit_file = $procedure_01_dir."/".$element;
		
		#top folder
		my $top_folder = $element;
		$top_folder =~ s/$suffix_01/$suffix_02/;
		$top_folder = $procedure_02_dir."/".$top_folder;
		$top_folder =~ s/.txt$//;
		
		#top file
		my $top_file = $element;
		$top_file =~ s/$suffix_01/$suffix_02/;
		$top_file =~ s/.txt$//;
		$top_file = $top_file."_AA.fas";
		
		#database name
		my $top_db = $element;
		my @top_db = split (/\./, $top_db);
		$top_db = $top_db[0];
		$top_db =~ s/$suffix_01$/$suffix_02/;
		
		#Counting the number of hits.
		my $lines = `wc -l $hit_file`;
		$lines =~ s/^\s*//g;
		my @hit_count = split (/ /, $lines);
		my $hit_count = $hit_count[0];
		
		#Handling files with no hits.
		if ($hit_count < 1) {
			mkdir $top_db;
			system "touch $top_file";
		}
		else {
			#Screening for unique hits.
			my $temp_02_01 = $query."_Proceture02_01.txt";
			&OrthoReD_library::uniq_finder ($hit_file, $temp_02_01);
			
			#Screening for top hits.
			my $temp_02_02 = $query."_Procedure02_02.txt";
			&OrthoReD_library::top_finder ($temp_02_01, $temp_02_02, $loci_threshold);
			
			#Making a fasta file of the hits.
			my $headers = `cut -f 2 $temp_02_02`;
			$headers =~ s/\r$//;
			$headers =~ s/\n$//;
			my @headers = split (/\n/, $headers);
			my $querytest = "ABSENT";
			foreach my $header (@headers) {
				if ($header =~ m/$query/) {
					$querytest = "PRESENT";
				}
				else {
					next;
				}
			}
			if ($querytest =~ m/ABSENT/) {
				push (@headers, $query);
			}
			$headers = \@headers;
			my $in_database = \%database_aa;
			&OrthoReD_library::fasta_maker ($headers, $in_database, $top_file);
			
			#Making a top database for all-against-all blast
			if ($blast_type =~ m/AB/) {
				&OrthoReD_library::database_maker_AB ($top_file, "AA", $top_db);
			}
			elsif ($blast_type =~ m/NCBI/ || $blast_type =~ m/SWIPE/) {
				&OrthoReD_library::database_maker_NCBI ($top_file, "AA", $top_db);
			}
			unlink ($temp_02_01);
			unlink ($temp_02_02);
		}
		system "mv $top_file $top_db";
		system "mv $top_db $procedure_02_dir";
	}
}
##########



#####(Procedure-03)For each TOP database, all-against-all BLASTP results will be generated.#####
#Setting file names
my $procedure_03 = "03";
my $suffix_03 = "_TOP_BLAST";
my $procedure_03_dir = $procedure_03.$suffix_03;

if ($procedures =~ m/03/) { if ($ls !~ m/$procedure_02_dir/ and $procedures !~ m/02/) { die "USAGE: Procedure-03 requires directory $procedure_02_dir.\n"};
	#making required folders
	mkdir $procedure_03_dir;
	my $id_03 = "Procedure03_".$id;
	
	#Selecting files
	my @top_dbs;
	foreach my $query (sort keys %querys) {
		my $file_name = $query.$suffix_02;
		push (@top_dbs, $file_name);
	}
	
	foreach my $element (@top_dbs) {
		#query
		my $query = $element;
		$query =~ s/$suffix_02//;
		
		#db folder
		my $top_db = $procedure_02_dir."/".$element."/".$element;
		
		#query file
		my $query_file = $procedure_02_dir."/".$element."/".$element."_AA.fas";
		
		#blast file
		my $blast_file = $element;
		$blast_file =~ s/$suffix_02/$suffix_03/;
		$blast_file = $procedure_03_dir."/".$blast_file.".txt";
		
		#Counting the number of file.
		my $element_dir = $procedure_02_dir."/".$element;
		$element_dir = `ls $element_dir`;
		
		#Handliioiing files with no hits.
		if (!-s $query_file) {
			system "touch $blast_file";
			next;
		}
		#Files with hits.
		else {
			#Running all-against-all BLASTP.
			my $blast_out;
			if ($blast_type =~ m/AB/) {
				$blast_out = &OrthoReD_library::blaster_AB ($query_file, $top_db, $w_threshold, $id_03);
			}
			elsif ($blast_type =~ m/NCBI/) {
				$blast_out = &OrthoReD_library::blaster_NCBI ($query_file, $top_db, $id_03, $threads);
			}
			elsif ($blast_type =~ m/SWIPE/) {
				$blast_out = &OrthoReD_library::blaster_NCBI ($query_file, $top_db, $id_03, $threads);
			}
			my $temp_03_01 = $query."_Procedure03_01.txt";
			&OrthoReD_library::blast_organizer ($blast_out, $temp_03_01);
			&OrthoReD_library::uniq_pair_finder ($temp_03_01, $blast_file);
			unlink ($blast_out);
			unlink ($temp_03_01);
		}
	}
}
##########



#####(Procedure-04)For each HITS_TOP_BLAST file, distance matrix (reciprocal e-values) of screened TOP_HITS will be generated, then is used to collect all the MCL members.#####
#Setting file names
my $procedure_04 = "04";
my $suffix_04 = "_TOP_MCL";
my $procedure_04_dir = $procedure_04.$suffix_04;

if ($procedures =~ m/04/) { if ($ls !~ m/$procedure_03_dir/ and $procedures !~ m/03/) { die "USAGE: Procedure-04 requires directory $procedure_03_dir.\n"};
	#making required folders
	mkdir $procedure_04_dir;
	
	#Selecting files
	my @top_blast;
	foreach my $query (sort keys %querys) {
		my $file_name = $query.$suffix_03.".txt";
		push (@top_blast, $file_name);
	}
	
	foreach my $element (@top_blast) {
		#query
		my $query = $element;
		my @query = split (/_/, $query);
		$query = $query[0]."_".$query[1]."_".$query[2]."_".$query[3];
		
		#blast file
		my $blast_file = $procedure_03_dir."/".$element;
		
		#mcl aa file
		my $mcl_aa_file = $element;
		$mcl_aa_file =~ s/$suffix_03/$suffix_04/;
		$mcl_aa_file =~ s/\.txt/_AA.fas/;
		$mcl_aa_file = $procedure_04_dir."/".$mcl_aa_file;
		
		#selecting the database for fasta making
		my $in_database = \%database_aa;
		
		#Screening out low quality hits.
		my $temp_04_01 = $query."_Procedure04_01.txt";
		&OrthoReD_library::blast_parser ($blast_file, $temp_04_01, $eval_threshold, $identity_threshold, $length_threshold, $sander_schneider);
		
		#Handling files with no hits.
		my $lines = `wc -l $temp_04_01`;
		$lines =~ s/^\s*//g;
		my @hit_count = split (/ /, $lines);
		my $hit_count = $hit_count[0];
		if ($hit_count < 1) {
			my @cluster;
			push (@cluster, $query);
			my $cluster = \@cluster;
			
			&OrthoReD_library::fasta_maker ($cluster, $in_database, $mcl_aa_file);
		}
		else {
			#preparing the file for running mcl
			my $temp_04_02 = $query."_Procedure04_02.txt";
			&OrthoReD_library::mcl_prepper ($temp_04_01, $temp_04_02, $query);
			
			my $cluster = &OrthoReD_library::mcl_runner ($temp_04_02, $query, $i, $query, $threads);
			
			&OrthoReD_library::fasta_maker ($cluster, $in_database, $mcl_aa_file);
			
			unlink ($temp_04_02);
		}
		
		unlink ($temp_04_01);
	}
}
##########



#####(Step-05)For each AA FASTA file, an aligned AA FASTA file is generated.#####
#Setting file names
my $procedure_05 = "05";
my $suffix_05 = "_MCL_AA_ALN";
my $procedure_05_dir = $procedure_05.$suffix_05;

if ($procedures =~ m/05/) { if ($ls !~ m/$procedure_04_dir/ and $procedures !~ m/04/) { die "USAGE: Procedure-05 requires directory $procedure_04_dir.\n"};
	#making required folders
	mkdir $procedure_05_dir;
	
	#Selecting files
	my @aa_files;
	foreach my $query (sort keys %querys) {
		my $file_name = $query.$suffix_04."_AA.fas";
		push (@aa_files, $file_name);
	}
	
	foreach my $element (@aa_files) {
		#aa file
		my $aa_file = $procedure_04_dir."/".$element;
		
		#aa aln file
		my $aln_file = $element;
		$aln_file =~ s/_AA//;
		$aln_file =~ s/$suffix_04/$suffix_05/;
		$aln_file = $procedure_05_dir."/".$aln_file;
		
		#query
		my $query = $element;
		my @query = split (/_/, $query);
		$query = $query[0]."_".$query[1]."_".$query[2]."_".$query[3];
		
		#Counting the number of sequences.
		my $lines = `wc -l $aa_file`;
		$lines =~ s/^\s*//g;
		my @sequence_count = split (/ /, $lines);
		my $sequence_count = $sequence_count[0] / 2;
		
		#Handling files with only one sequence.
		if ($sequence_count == 1) {
			system "cp $aa_file $aln_file";
			next;
		}
		
		#Handling other files.
		my ($temp_log_maffter) = &OrthoReD_library::maffter ($aa_file, $aln_file, "AA", $threads, $query);
	}
}
##########



#####(Step-06)For each aligned AA FASTA file, a rooted maximum likelihood tree is generated.#####
#Setting file names
my $procedure_06 = "06";
my $suffix_06 = "_MCL_AA_TRE";
my $procedure_06_dir = $procedure_06.$suffix_06;

if ($procedures =~ m/06/) { if ($ls !~ m/$procedure_05_dir/ and $procedures !~ m/05/) { die "USAGE: Procedure-06 requires directory $procedure_05_dir.\n"};
	#making required folders
	mkdir $procedure_06_dir;
	
	#Selecting files
	my @aln_files;
	foreach my $query (sort keys %querys) {
		my $file_name = $query.$suffix_05.".fas";
		push (@aln_files, $file_name);
	}
	
	foreach my $element (@aln_files) {
		#aa aln file
		my $aln_file = $procedure_05_dir."/".$element;
		
		#tre file
		my $tre_file = $element;
		$tre_file =~ s/$suffix_05/$suffix_06/;
		$tre_file =~ s/.fas/.tre/;
		$tre_file = $procedure_06_dir."/".$tre_file;
		
		#gene of interest
		my $gene = $element;
		$gene =~ s/$suffix_05//;
		$gene =~ s/.fas//;
		
		#assigning file names
		my $conv = $gene."_conv.txt";
		my $aln_temp = $gene."_aln_reduced.fas";
		my $temp_dist = $gene."_dist.txt";
		my $tree_u = $gene."_tree_u.tre";
		my $tree_u2 = $gene."_tree_exu.tre"; 
		
		#reducing identical sequences
		&OrthoReD_library::seq_reducer ($aln_file, $aln_temp, $conv);
		
		#Counting the number of sequences.
		my $lines = `wc -l $aln_temp`;
		$lines =~ s/^\s*//g;
		my @sequence_count = split (/ /, $lines);
		my $sequence_count = $sequence_count[0] / 2;
		
		#Handling files with three sequences or less.
		my $root = "midpoint";
		if ($sequence_count < 4) {
			&OrthoReD_library::mini_tree_maker ($aln_temp, $tree_u);
			&OrthoReD_library::tree_expander ($tree_u, $conv, $tree_u2);
			system "touch $temp_dist";
			system "cp $tree_u2 $tre_file";
		}
		else {
			my $test = &OrthoReD_library::raxmler_speed ($aln_temp, $tree_u, $threads, $gene, $vraxml);
			#recovering the reduced tips
			&OrthoReD_library::tree_expander ($tree_u, $conv, $tree_u2);
			
			if ($test =~ m/GOOD/) {
				#picking a given tip as the root
				$root = &OrthoReD_library::root_finder ($gene, $tree_u2, $temp_dist, \@ogs, $spp_list, $rooting);
				&OrthoReD_library::rooter ($tree_u2, $root, $tre_file);
			}
			elsif ($test =~ m/BACKUP/) {
				system "touch $temp_dist";
				system "cp $tree_u2 $tre_file";
			}
		}
		unlink ($conv);
		unlink ($aln_temp);
		unlink ($temp_dist);
		unlink ($tree_u);
		unlink ($tree_u2);
		system "rm -rf $aln_temp.reduced";
	}
}
###########



#####(Step-07)For each TRE file, a single-lined FASTA file with orthologs is generated.#####
my $procedure_07 = "07";
my $suffix_07;
if ($db_seq_type =~ m/DNA/) {
	$suffix_07 = "_ORTHO_DNA";
}
if ($db_seq_type !~ m/DNA/) {
	$suffix_07 = "_ORTHO_AA";
}
my $procedure_07_dir = $procedure_07.$suffix_07;

if ($procedures =~ m/07/) { if ($ls !~ m/$procedure_06_dir/ and $procedures !~ m/06/) { die "USAGE: Procedure-07 requires directory $procedure_06_dir.\n"};
	#making required folders
	mkdir $procedure_07_dir;
	
	#Selecting files
	my @tre_files;
	foreach my $query (sort keys %querys) {
		my $file_name = $query.$suffix_06.".tre";
		push (@tre_files, $file_name);
	}
	
	#preparing species list
	my @prefixes;
	foreach my $key (sort keys %spp1) {
		push(@prefixes, $key);
	}
	my $prefixes = \@prefixes;
	
	foreach my $element (@tre_files) {
		#tre file
		my $tree_file = $procedure_06_dir."/".$element;
		
		#preparing gene name
		my $gene = $element;
		my @gene = split (/_/, $gene);
		$gene = $gene[0]."_".$gene[1]."_".$gene[2]."_".$gene[3];
		
		#ortho file
		my $ortho_file = $procedure_07_dir."/".$gene.$suffix_07.".fas";
		
		#preparing parenthetical tree
		my $cut_tre = $gene."_cut.tre";
		&OrthoReD_library::tree_cutter ($tree_file, $gene, $branch_threshold, $cut_tre);
		my $tree = &OrthoReD_library::branch_length_stripper ($cut_tre);
		
		
		#make a tree object
		my $tree_object = &OrthoReD_library::toTreeObj($tree, $prefixes, $overlap_threshold);
		
		#select a set of orthologs
		my %pre_orthos;
		my @orthogroups = &OrthoReD_library::orthologGroups($tree_object, $prefixes);

		foreach my $group (@orthogroups) {
			my @group = @$group;
			foreach my $tip (@group) {
				if ($tip =~ m/$gene/) {
					foreach my $ortho (@group) {
						my @key = split (/_/, $ortho);
						my $key = $key[0]."_".$key[1];
						$pre_orthos{$key} = "TRUE";
					}
				}
			}
		}
		
		#generate a single lined FASTA file from the selected orthologs with its splice variants.
		#add all the splice variants into the list of orthologs.
		my @orthologs;
		push (@orthologs, $gene);
		
		my $in_database;
		if ($db_seq_type =~ m/DNA/) {
			$in_database = \%database_dna;
		}
		if ($db_seq_type !~ m/DNA/) {
			$in_database = \%database_aa;
		}
		open (DATABASE, "<$database_aa_file");
		while (my $line = <DATABASE>) {
			if ($line =~ m/^>/) {
				$line =~ s/\r//sig;
				$line =~ s/\n//sig;
				my $header = $line;
				$header =~ s/^>//;
				my @header = split (/_/, $header);
				$header = $header[0]."_".$header[1];
				if (exists $pre_orthos{$header}) {
					$line =~ s/^>//;
					push(@orthologs, $line);
				}
			}
		}
		close (DATABASE);
		
		my @unique;
		my %seen;
		foreach my $value (@orthologs) {
			if (! $seen{$value}) {
				push @unique, $value;
				$seen{$value} = 1;
			}
		}
		@orthologs = sort @unique;
		
		my $orthologs = \@orthologs;
		
		&OrthoReD_library::fasta_maker ($orthologs, $in_database, $ortho_file);
		unlink ($cut_tre);
	}
}
##########



#####(Procedure-08)Summarizing the output.#####
my $summary_file = "08_SUMMARY.txt";
my $rank_file = "08_RANKING.txt";

if ($procedures =~ m/08/) { if ($ls !~ m/$procedure_07_dir/ and $procedures !~ m/07/) { die "USAGE: Procedure-08 requires directory $procedure_07_dir.\n"};
	#summary
	if (!$conversion_list) {
		&OrthoReD_library::summarizer ($spp_list, $procedure_07_dir, $summary_file);
	}
	else {
		&OrthoReD_library::summarizer ($spp_list, $procedure_07_dir, $summary_file, $conversion_list);
	}
	
	#ranking
	&OrthoReD_library::ranker ($procedure_07_dir, $procedure_01_dir, $rank_file);
}
##########



#####(Procedure-09)Formatting directories.#####
if ($procedures =~ m/09/) { 
	system "mv $procedure_01_dir $product_dir 2>/dev/null";
	system "mv $procedure_02_dir $product_dir 2>/dev/null";
	system "mv $procedure_03_dir $product_dir 2>/dev/null";
	system "mv $procedure_04_dir $product_dir 2>/dev/null";
	system "mv $procedure_05_dir $product_dir 2>/dev/null";
	system "mv $procedure_06_dir $product_dir 2>/dev/null";
	system "mv $procedure_07_dir $product_dir 2>/dev/null";
	system "mv $summary_file $product_dir 2>/dev/null";
	system "mv $rank_file $product_dir 2>/dev/null";
	
	if ($logprint =~ m/YES/) {
		print LOG "SCRIPT COMPLETE\n";
		close (LOG);
		system "mv $log $product_dir";
	}
}
##########
__END__
