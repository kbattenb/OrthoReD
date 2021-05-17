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
#Script "step-03-04.pl" runs a quality check on all subjects and generates a single-lined FASTA file with all subjects and a BLAST database file.
#All output files will be in the working directory.
##########



#####Options#####
#setting the default values.
my $help = "";
my $query = ""; #file path to the query file with all the queries.
my $q_seq_type = ""; #"DNA" or "AA".
my $subject = ""; #file path with the subject files.
my $s_seq_type = ""; #"DNA" or "AA".
my $s_clean = "YES"; #"YES" to clean the headers or "NO" to keep headers as is.
my $s_full_length = "NO"; ##"YES" to only keep full-length sequences or "NO" to not have this screen.
my $s_qc_handle = "REMOVE"; #"KEEP" to keep all subject sequences or "REMOVE" to automatically remove problematic sequences or "REPORT" to report and kill the process when such sequence is found.
my $database_name = ""; #name of the database.
my $blast_type = "BOTH"; #"AB" for AB-BLAST, "NCBI" for NCBI blast, and "BOTH" for both.
my $spp_list = ""; #File with the species list.

#making the options into external arguments.
GetOptions (
	'help' => \$help,
	'query=s' => \$query,
	'q_seq_type=s' => \$q_seq_type,
	'subject=s@' => \$subject,
	's_seq_type=s' => \$s_seq_type,
	's_clean=s' => \$s_clean,
	's_full_length=s' => \$s_full_length,
	's_qc_handle=s' => \$s_qc_handle,
	'database_name=s' => \$database_name,
	'blast_type=s' => \$blast_type,
	'spp_list=s' => \$spp_list
);

#printing help
if ($help) {
	print "\n";
	print "\tstep-03-04.pl generates BLAST database files with all the queries and subjects formatted for the downstream processes of OrthoReD based on all the query and subject sequence file(s) provided.\n\n";
	print "\t--help\tPrint this massage and die.\n\n";
	print "\t--query\tPath to a single file with all the formatted query sequences.\n\n";
	print "\t--q_seq_type\tSequence type of the query. Values 'DNA' or 'AA' are accepted.\n\n";
	print "\t--subject\tFile path to the directory with the subject files.\n\n";
	print "\t--s_seq_type\tSequence type of the subject. Values 'DNA' or 'AA' are accepted.\n\n";
	print "\t--s_clean\t'YES' to format the headers of the subject sequences using the module <header_cleaner.pm> or 'NO' to leave the headers unchanged. (DEFAULT: YES)\n\n";
	print "\t--s_qc_handl\t'KEEP' to keep all subject sequences or 'REMOVE' to automatically remove any sequence(s) with unclear reading frame. (DEFAULT: REMOVE)\n\n";
	print "\t--s_full_length\t'YES' to only keep full-length sequences or 'NO' to not have this screen. (DEFAULT: NO)\n\n";
	print "\t--database_name\tUser-defined name for the database being generated.\n\n";
	print "\t--blast_type\tType of BLAST database generated. 'AB' for only generating AB-BLAST database, 'NCBI' for only generating NCBI-BLAST database, and 'BOTH' for generating both databases. (DEFAULT: BOTH)\n\n";
	print "\t--spp_list\tPath to the file with the lsit of species. See spp_list_EXAMPLE.txt for the proper format.\n";
	die "\n";
}

#checking for required options.
if (!$query) {
	die "USAGE: option --query is required.\n";
}
if ($q_seq_type !~ m/^DNA$/ && $q_seq_type !~ m/^AA$/) {
	die "USAGE: option --q_seq_type must be set as DNA or AA.\n";
}
if (!$subject) {
	die "USAGE: option --subject is required.\n";
}
if ($s_seq_type !~ m/^DNA$/ && $s_seq_type !~ m/^AA$/) {
	die "USAGE: option --s_seq_type must be set as DNA or AA.\n";
}
if (!$database_name) {
	die "USAGE: option --database_name is required.\n";
}
if ($blast_type !~ m/^AB$/ && $blast_type !~ m/^NCBI$/ && $blast_type !~ m/^BOTH$/) {
	die "USAGE: option --blast_type must be set as AB, NCBI, or BOTH.\n";
}
if (!$spp_list) {
	die "USAGE: option --spp_list is required.\n";
}
##########



#####(Procedure-00)Setting up the environment.#####

###Identifying the working directory.
my $pwd = cwd();
my $pwd_ortho;
my @ls;
my @ls_ortho;

opendir my $dh, $pwd or die "Couldn't open dir '$pwd': $!";
@ls = readdir $dh;
closedir $dh;
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
	die "USAGE: <step-03-04.pl> requires folder <OrthoReD_vXXXXXXXX> to be in the working directory.\n";
}
else {
	$test_st = 1;
}

#"step-03-04.pl" script for this process.
foreach my $file (@ls_ortho) {
	if ($file =~ m/step-03-04.pl/) {
		$test_st = 0;
		last;
	}
}
if ($test_st ne 0) {
	die "USAGE: <step-03-04.pl> needs to be in folder <OrthoReD_vXXXXXXXX>.\n";
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
	die "USAGE: <step-03-04.pl> requires file <OrthoReD_library.pm> to be in folder <OrthoReD_vXXXXXXXX>.\n";
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
	die "USAGE: <step-03-04.pl> requires file <OrthoReD_library.pm> to be in folder <OrthoReD_vXXXXXXXX>.\n";
}
else {
	$test_st = 1;
}
###

###setting up to use the libraries
use OrthoReD_v20180830::OrthoReD_library;
use OrthoReD_v20180830::header_cleaner;
###

##########



#####Programs required#####
#AB-BLAST
if(!defined $ENV{'ABBLAST'} && $blast_type !~ m/NCBI/) {
	die "USAGE: <step-03-04.pl> requires environmental variable ABBLAST to be set as the path to the tools for AB-BLAST unless --blast_type is set to NCBI.\n";
}
my $blast_check1 = `which $ENV{'ABBLAST'}/xdformat`;
if (!$blast_check1 && $blast_type !~ m/NCBI/) {
	die "USAGE: <step-03-04.pl> requires AB-BLAST to be running unless --blast_type is set to NCBI.\n";
}
#NCBI-BLAST
if(!defined $ENV{'NCBIBLAST'} && $blast_type !~ m/AB/) {
	die "USAGE: <step-03-04.pl> requires environmental variable NCBIBLAST to be set as the path to the tools for NCBIBLAST unless --blast_type is set to AB.\n";
}
my $blast_check2 = `which $ENV{'NCBIBLAST'}/segmasker`;
if (!$blast_check2 && $blast_type !~ m/AB/) {
	die "USAGE: <step-03-04.pl> requires NCBI-BLAST to be running unless --blast_type is set to AB.\n";
}

##########



#####Generating required files#####
#Making the subject directory
my $subject_dir = "Step-03_Subjects_RAW";
rmtree($subject_dir);
mkdir $subject_dir;

#Making the database directory
my $product_dir = "Step-04_DATABASE";
rmtree($product_dir);
mkdir $product_dir;

#Making a log file
my $log = "step-03-04_LOG.txt";
open (LOG, ">$log") or die "cannot open $log.\n";
print "RUNNING SCRIPT: $pwd_ortho step-03-04.pl\n";
print LOG "RUNNING SCRIPT: $pwd_ortho step-03-04.pl\n";
print "(Procedure-00)Setting up the environment.\n";
print LOG "(Procedure-00)Setting up the environment.\n";

#Identify all the species included in the the species list.
my $full_name = 0;
my $short_name = 1;
my $acronym_name = 2;

my %spp1 = %{ &OrthoReD_library::species_lister ($spp_list, $acronym_name, $full_name) };
my %spp2 = %{ &OrthoReD_library::species_lister ($spp_list, $short_name, $full_name) };

print "Following species are included in the species list:\n";
print LOG "Following species were included in the species list:\n";
foreach my $sp (sort keys %spp1) {
	print "\t$spp1{$sp}\n";
	print LOG "\t$spp1{$sp}\n";
}

#Identifying all the raw subject files and copy them into directory Step-03_Subjects_RAW.
my @subjects = @{ $subject };
print "Following files are included in the database:\n";
print LOG "Following files were included in the database:\n";

foreach my $file (@subjects) {
	print "\t$file\n";
	print LOG "\t$file\n";
	
	my @file_name = split (/\//, $file);
	my $new_file = $pwd."/".$subject_dir."/".$file_name[-1];
	
	system "cp $file $new_file";
}
##########



#####(Procedure-01)Foreach raw subject file, a quality-checked single-lined FASTA file is generated.#####
print "(Procedure-01)Foreach raw subject file, generate a quality-checked single-lined FASTA file.\n";
print LOG "(Procedure-01)Foreach raw subject file, generate a quality-checked single-lined FASTA file.\n";

#Setting file names
my $procedure_01 = "01";
my $suffix_01 = "_CLEANED";
my $procedure_01_dir = $procedure_01.$suffix_01;

mkdir $procedure_01_dir;

#finding all subject files
my @multi_files;
#my $dir = $pwd."/".$subject_dir;
#opendir $dh, $dir or die "Couldn't open dir '$dir': $!";
#@multi_files = readdir $dh;
#closedir $dh;
#shift @multi_files;
#shift @multi_files;
my $multi_files = `ls $pwd/$subject_dir`;
@multi_files = split (/\n/, $multi_files);

print "\n\nDo the files listed below fit what you actually see in $pwd/$subject_dir?\n";
foreach my $file (@multi_files) {
	print "\t$file\n";
}
die "\n\n\n";

#processing subject files
print "processing:\n";
print LOG "processed:\n";

my @cleaned_files;
foreach my $file (@multi_files) {
	#multi file
	my $multi_file = $subject_dir."/".$file;
	
	#Recognizing the genome
	my @file_name = split (/\//, $multi_file);
	my $sp_name = $file_name[1];
	my @sp_name = split (/_/, $sp_name);
	$sp_name = $sp_name[0];
	
	#clean file
	my $cleaned_file;
	if (exists $spp2{$sp_name}) {
		$cleaned_file = $procedure_01_dir."/".$sp_name.$suffix_01."\.fas";
	}
	else {
		my @multi_file = split (/\//, $multi_file);
		$cleaned_file = $procedure_01_dir."/".$multi_file[1];
		$cleaned_file =~ s/\.fas/$suffix_01\.fas/;
	}
	push (@cleaned_files, $cleaned_file);
	
	print "\t$multi_file\n";
	print LOG "\t$multi_file\n";
	
	&OrthoReD_library::singler ($multi_file, "temp.fas");
	
	if ($s_clean =~ m/YES/) {
		open (TEMP, "<temp.fas") or die "cannot open temp.fas.\n";
		open (CLEAN, ">$cleaned_file") or die "$cleaned_file.\n";
		while (my $line = <TEMP>) {
			$line =~ s/\r//sig;
			$line =~ s/\n//sig;
			if ($line =~ m/^>/) {
				my $header = $line;
				my ($cleaned_header) = &header_cleaner::header_cleaner ($spp_list, $header);
				print CLEAN "$cleaned_header\n";
				next;
			}
			print CLEAN "$line\n";
		}
		close (TEMP);
		close (CLEAN);
	}
	if ($s_clean =~ m/NO/) {
		system "cp temp.fas $cleaned_file";
	}
	
	unlink ("temp.fas");
}
##########



#####(Procedure-02)From all the cleaned files, a database file is generated.#####
print "(Procedure-02)From all the cleaned files, generate a database file.\n";
print LOG "(Procedure-02)From all the cleaned files, generate a database file.\n";

#Setting file names
my $procedure_02 = "02";
my $suffix_02 = "_DATABASE";
my $procedure_02_dir = $procedure_02.$suffix_02;
my $abdir = $procedure_02_dir."/AB";
my $ncbidir = $procedure_02_dir."/NCBI";

my $database_file = $procedure_02_dir."/".$database_name."_".$s_seq_type.".fas";
my $removed_file = $procedure_02_dir."/".$database_name."_REMOVED.fas";

mkdir ($procedure_02_dir);
mkdir ($abdir);
mkdir ($ncbidir);

#Concatenating files.
print "Generating database file from:\n";
print LOG "Generated database file from:\n";
foreach my $file (@cleaned_files) {
	print "\t$file\n";
	print LOG "\t$file\n";
}

#Generating one file with all subjects.
my $subject_files = join (" ", @cleaned_files);
system "cat $subject_files > temp1.fas";

#Running quality check
if ($s_qc_handle =~ m/^KEEP$/) {
	system "mv temp1.fas temp2.fas";
}
else {
	if ($s_full_length =~ m/YES/) {
		print "Sequences that are not full-length CDS are removed.\n";
		print LOG "Sequences that are not full-length CDS were removed.\n";
	}
	elsif ($s_full_length =~ m/NO/) {
		print "Sequences with undetermined reading frame are removed.\n";
		print LOG "Sequences with undetermined reading frame were removed.\n";
	}
	else {
		die "ERROR: --s_full_length needs to be YES or NO.\n";
	}
	
	my ($value) = &OrthoReD_library::qcer_file ("temp1.fas", $s_seq_type, $s_full_length);
	
	if ($value =~ m/^TRUE$/) {
		system "touch $removed_file";
		system "mv temp1.fas temp2.fas";
		print "0% of the sequences were removed from database.\n";
		print LOG "0% of the sequences were removed from database.\n";
	}
	else {
		if ($s_qc_handle =~ m/^REPORT$/) {
			die "$value\n";
		}
		elsif ($s_qc_handle =~ m/^REMOVE$/) {
			my ($percent_removed) = &OrthoReD_library::separator ($spp_list, "temp1.fas", "temp2.fas", $removed_file, $s_seq_type, $s_full_length);
			unlink ("temp1.fas");
			print "$percent_removed\% of the sequences were removed from database.\n";
			print LOG "$percent_removed\% of the sequences were removed from database.\n";
		}
		else {
			die "ERROR: --s_qc_handle needs to be KEEP, REPORT, or REMOVE.\n";
		}
	}
}
if ($q_seq_type =~ m/AA/ && $s_seq_type =~ m/DNA/) {
	&OrthoReD_library::translator2("temp2.fas", "trans.fas");
	unlink ("temp2.fas");
	system "mv trans.fas temp2.fas";
}

#Adding all query sequences that has not been included into the database.
print "Unique sequences from the following query file is included in the database:\n";
print LOG "Unique sequences from the following query file was included in the database:\n";
print "\t$query\n";
print LOG "\t$query\n";

#Making a list of headers in the query.
my %querys = %{ &OrthoReD_library::fasta2hash($query) };
my $querys = \%querys;
my @q_headers = keys %querys;
my $q_headers = \@q_headers;

#Making a list of headers in the subject.
my %temp2 = %{ &OrthoReD_library::fasta2hash("temp2.fas") };
my @s_headers = keys %temp2;
my $s_headers = \@s_headers;

#Adding all the query sequences that were not found in the subject into the subject sequences.
my ($common, $mismatch1, $mismatch2) = &OrthoReD_library::mismatch_finder ($q_headers, $s_headers);
&OrthoReD_library::fasta_maker ($mismatch1, $querys, "temp3.fas");
if ($q_seq_type =~ m/DNA/ && $s_seq_type =~ m/AA/) {
	&OrthoReD_library::translator2("temp3.fas", "trans.fas");
	unlink ("temp3.fas");
	system "mv trans.fas temp3.fas";
}

system "cat temp2.fas temp3.fas > $database_file";
unlink ("temp2.fas");
unlink ("temp3.fas");

#Making a masked BLAST database file.
print "creating a database from:\n";
print LOG "created a database from:\n";
print "\t$database_file\n\n";
print LOG "\t$database_file\n\n";

my $dna_database_name = $database_name."_DNA";
my $aa_database_name = $database_name."_AA";
if ($blast_type =~ m/AB/ or $blast_type =~ m/BOTH/) {
	if ($q_seq_type =~ m/DNA/ and $s_seq_type =~ m/DNA/) {
		my $database_file_aa = $database_file;
		$database_file_aa =~ s/_DNA.fas/_AA.fas/;
		&OrthoReD_library::translator2($database_file, $database_file_aa);
		
		&OrthoReD_library::database_maker_AB ($database_file, "DNA", $dna_database_name, $ENV{'ABBLAST'});
		&OrthoReD_library::database_maker_AB ($database_file_aa, "AA", $aa_database_name, $ENV{'ABBLAST'});
	}
	else {
		&OrthoReD_library::database_maker_AB ($database_file, "AA", $aa_database_name, $ENV{'ABBLAST'});
	}
	
	system "mv $database_name* $abdir";
}
if ($blast_type =~ m/NCBI/ or $blast_type =~ m/BOTH/) {
	if ($q_seq_type =~ m/DNA/ and $s_seq_type =~ m/DNA/) {
		my $database_file_aa = $database_file;
		$database_file_aa =~ s/_DNA.fas/_AA.fas/;
		&OrthoReD_library::translator2($database_file, $database_file_aa);
		
		&OrthoReD_library::database_maker_NCBI ($database_file, "DNA", $dna_database_name, $ENV{'NCBIBLAST'});
		&OrthoReD_library::database_maker_NCBI ($database_file_aa, "AA", $aa_database_name, $ENV{'NCBIBLAST'});
	}
	else {
		&OrthoReD_library::database_maker_NCBI ($database_file, "AA", $aa_database_name, $ENV{'NCBIBLAST'});
	}
	
	system "mv $database_name* $ncbidir";
}
print "Database is generated.\n";
print LOG "Database was generated.\n";
##########



#####(Procedure-03)Formatting directories.#####
print "(Procedure-03)Formatting directories.\n";
print LOG "(Procedure-03)Formatting directories.\n";

#Formatting directories.
system "mv $procedure_01_dir $product_dir";
system "mv $procedure_02_dir $product_dir";

print "SCRIPT COMPLETE\n";
print LOG "SCRIPT COMPLETE\n";
close (LOG);
system "mv $log $product_dir";
##########
__END__