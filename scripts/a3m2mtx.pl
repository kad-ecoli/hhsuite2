#!/usr/bin/env perl
# a3m2mtx.pl
# Convert MSA to PSI-BLAST profile.
# This script is based on addss.pl from HHsuite version 2.0.16

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use File::Temp qw/ tempfile tempdir /;
use strict;


$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;              # verbose mode

my $informat="a3m";    # input format
my $neff = 0;          # use alignment with this diversity for PSIPRED prediction

my $help="
a3m2mtx.pl 
convert MSA to PSI-BLAST profile in MTX format.

Allowed input formats are A3M (default), A2M/FASTA (-fas, -a2m),
CLUSTAL (-clu), STOCKHOLM (-sto), PSICOV (-aln).

By default, the original, unfiltered MSA is directly converted to MTX. 
However, it was reported previously that if the MSA filtered to neff 7
before being converted to MTX, the PSIPRED secondary structure prediction
using this filtred MTX can be improved. This behavior can be specified by
-neff 7

Usage: a3m2mtx.pl <ali_file> [<outfile>] [-fas|-a3m|-clu|-sto|-aln] [-neff 7]
\n";

# Variable declarations
my $line;
my @seqs;              # sequences from infile (except >aa_ and >ss_pred sequences)
my $query_length;
my $header;            # header of MSA: everything before first '>'
my $name;              # query in fasta format: '>$name [^\n]*\n$qseq\n'
my $qseq;              # residues of query sequence
my $infile;
my $outfile;

###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

#Input format fasta?
if    ($options=~s/ -fas\s/ /g) {$informat="fas";}
elsif ($options=~s/ -a2m\s/ /g) {$informat="a2m";}
elsif ($options=~s/ -a3m\s/ /g) {$informat="a3m";}
elsif ($options=~s/ -clu\s/ /g) {$informat="clu";}
elsif ($options=~s/ -sto\s/ /g) {$informat="sto";}
elsif ($options=~s/ -aln\s/ /g) {$informat="aln";}

if ($options=~s/ -v\s+(\d+) / /g)  {$v=$1;}

# Set input and output file
if ($options=~s/ -i\s+(\S+) //) {$infile=$1;}
if ($options=~s/ -o\s+(\S+) //) {$outfile=$1;}
if ($options=~s/ -neff\s+(\S+) //) {$neff=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$infile=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$outfile=$1;}

# Warn if unknown options found or no infile/outfile
if ($options!~/^\s*$/) {$options=~s/^\s*(.*?)\s*$/$1/g; die("Error: unknown options '$options'\n");}
if (!$infile) {print($help); exit(1);}

my $v2 = $v-1;
if ($v2>2) {$v2--;}
if ($v2<0) {$v2=0;}

###############################################################################################
# Reformat input alignment to a3m and psiblast-readable format and generate file with query sequence
###############################################################################################

my $inbase; # $inbasename of infile: remove extension
my $inroot; # $inbasename of infile: remove path and extension
if ($infile=~/(.*)\..*/) {$inbase=$1;} else {$inbase=$infile;}  # remove extension
if ($inbase=~/.*\/(.*)/)  {$inroot=$1;} else {$inroot=$inbase;} # remove path 

# Create tmpfile
my $tmpdir;
if ($v<=3) {$tmpdir = tempdir( CLEANUP => 1);} else {$tmpdir = tempdir( CLEANUP => 0);}
my ($tmpf, $tmpfile) = tempfile( DIR => $tmpdir );
my $tmpfile_no_dir;
if ($tmpfile=~/.*\/(.*)/)  {$tmpfile_no_dir=$1;} else {$tmpfile_no_dir=$tmpfile;} # remove path 



############################################################################################
if (!$outfile) {$outfile="$inbase.mtx";}

# Use first sequence to define match states and reformat input file to a3m and psi
if ($informat eq "a3m") {
    &HHPaths::System("cp $infile $tmpfile.in.a3m");
} elsif ($informat eq "aln") {
    &HHPaths::System("sed = $infile |sed 'N;s/\\n/\\t/'|sed 's/^/>/g'|sed 's/\\t/\\n/g' > $tmpfile.in.a3m");
} else {
    &HHPaths::System("$hhscripts/reformat.pl -v $v2 -M first $informat a3m $infile $tmpfile.in.a3m");
}

# Read query sequence
open (INFILE, "<$tmpfile.in.a3m") or die ("ERROR: cannot open $tmpfile.in.a3m!\n");
$/=">"; # set input field separator
my $i=0;
$qseq="";
$header = <INFILE>;
$header =~s />$//; 
while ($line=<INFILE>) {
    $line=~s/>$//;
    if ($line=~/^ss_/ || $line=~/^aa_/) {next;}
    $seqs[$i++]=">$line";
    if(!$qseq) {
        $line=~s/^(.*)[^\n]*//;
        $name=$1;
        $qseq=$line;
        $qseq=~s/\n//g;
    }
}
close(INFILE);
$/="\n"; # set input field separator

if ($qseq =~ /\-/) {
    
    # First sequence contains gaps => calculate consensus sequence
    &HHPaths::System("$hhbin/hhconsensus -i $tmpfile.in.a3m -s $tmpfile.sq -o $tmpfile.in.a3m > /dev/null");
    
} else {
    
    $query_length = ($qseq=~tr/A-Z/A-Z/);
    $qseq=~tr/A-Z//cd; # remove everything except capital letters
    
    # Write query sequence file in FASTA format
    open (QFILE, ">$tmpfile.sq") or die("ERROR: can't open $tmpfile.sq: $!\n");
    printf(QFILE ">%s\n%s\n",$name,$qseq);
    close (QFILE);
}

# Filter alignment to diversity $neff 
if ($neff > 0) {
    if ($v>=1) {printf ("Filtering alignment to diversity $neff ...\n");}
    &HHPaths::System("$hhbin/hhfilter -v $v2 -neff $neff -i $tmpfile.in.a3m -o $tmpfile.in.a3m");
}

# Reformat into PSI-BLAST readable file for jumpstarting 
&HHPaths::System("$hhscripts/reformat.pl -v $v2 -r -noss a3m psi $tmpfile.in.a3m $tmpfile.in.psi");

# PSI-BLAST on dummy db
&RunPseudoPsiblast("$tmpfile.sq");
&HHPaths::System("cp $tmpfile.mtx $outfile");

if ($v<=3) {
    unlink("$tmpfile.mtx");
    unlink("$tmpfile.in.a3m");
    unlink("$tmpfile.in.psi");
} 

exit;
    
##############################################################################################
# Run blastpgp jump starting from alignment in $tmpfile.in.psi (called by BuildAlignment)
##############################################################################################
sub RunPseudoPsiblast() {
    # Note that it assumes that the
    # following programs are in the appropriate directories:
    # blastpgp - PSIBLAST executable (from NCBI toolkit)
    # makemat - IMPALA utility (from NCBI toolkit)
    # psipred - PSIPRED V2 program
    # psipass2 - PSIPRED V2 program
    
    my $infile=$_[0];
    my $basename;  #file name without extension
    my $rootname;  #basename without directory path
    if ($infile =~/^(.*)\..*?$/)  {$basename=$1;} else {$basename=$infile;}
    if ($basename=~/^.*\/(.*?)$/) {$rootname=$1;} else {$rootname=$basename;}

    # Does dummy database exist?
    if (!-e "$dummydb.phr") {
	if (!-e "$dummydb") {die "Error in a3m2mtx.pl: Could not find $dummydb\n";}

	&HHPaths::System("cp $infile $dummydb");
	&HHPaths::System("$ncbidir/formatdb -i $dummydb");
	if (!-e "$dummydb.phr") {die "Error in a3m2mtx.pl: Could not find nor create index files for $dummydb\n";}
    }

    # Start Psiblast from checkpoint file tmp.chk that was generated to build the profile
    &HHPaths::System("$ncbidir/blastpgp -b 1 -j 1 -h 0.001 -d $dummydb -i $infile -B $tmpfile.in.psi -C $tmpfile.chk 1> $tmpfile.blalog 2> $tmpfile.blalog");
    
    #print("Predicting secondary structure...\n");
    
    &HHPaths::System("echo "."$tmpfile_no_dir".".chk > $tmpfile.pn\n");
    &HHPaths::System("echo "."$tmpfile_no_dir".".sq  > $tmpfile.sn\n");
    &HHPaths::System("$ncbidir/makemat -P $tmpfile");
    
    # Remove temporary files
    if ($v<=3) { unlink(split ' ', "$tmpfile.pn $tmpfile.sn $tmpfile.mn $tmpfile.chk $tmpfile.blalog $tmpfile.aux $tmpfile.sq");}
    return;
}
