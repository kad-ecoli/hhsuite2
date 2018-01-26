#!/usr/bin/env perl
# hhb_psipred.pl
# Predict secondary structure using hhblits alignment and psipred.
# This script is based on addss.pl from HHsuite version 2.0.16


use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use File::Temp qw/ tempfile tempdir /;
use strict;

my $help="
hhb_psipred.pl
predict secondary structure from sequence, using hhblits and psipred

Usage: hhb_psipred.pl <query_file> <db> [<outfile_ss2>] [<outfile_mtx>] [<outfile_solv>]

<query_file>   single sequence fasta input.
<db>           hhsuite format sequence database, e.g. uniprot20 or uniclust
<outfile_ss2>  psipred secondary structure prediction output.
<outfile_mtx>  blastpgp (legacy PSI-BLAST) format sequence profile
<outfile_solv> solvpred solvent accessibility prediction output.

Note that if <outfile_mtx> already exists, it will be directly used instead
of being re-generated from hhblits alignment.
\n";

if (@ARGV<2) {die ($help);}

#### search parameters ####
my $cpu=1;
my $neff=7;

#### parse comand line argument ####
my $options="";
for (my $i=0; $i<@ARGV; $i++) {$options.=" $ARGV[$i] ";}

my $infile;
my $db;
my $outfile_ss2;
my $outfile_mtx;
my $outfile_solv;

if ($options=~s/^\s*([^-]\S*) //)  {$infile=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$db=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$outfile_ss2=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$outfile_mtx=$1;}
if ($options=~s/^\s*([^-]\S*) //)  {$outfile_solv=$1;}

# check if database is valid
if (! -s "$db\_a3m_db") {die "ERROR! No such file $db\_a3m_db";}

my $inbase; # $inbasename of infile: remove extension
if ($infile=~/(.*)\..*/) {$inbase=$1;} else {$inbase=$infile;}  # remove extension

#### make tmp folder ####
my $tmpdir = tempdir( CLEANUP => 1);
my ($tmpf, $tmpfile) = tempfile( DIR => $tmpdir );
&HHPaths::System("cp $infile $tmpfile.fasta");

#### make profile ####
if (!-s "$outfile_mtx")
{
    print "Creating sequence profile\n";
    #&HHPaths::System("$hhbin/hhblits -i $tmpfile.fasta -d $db -n 3 -o $tmpfile.hhr -oa3m $tmpfile.a3m -cpu $cpu");
    #&HHPaths::System("$hhscripts/a3m2mtx.pl $tmpfile.a3m $tmpfile.mtx -neff $neff");
    &HHPaths::System("$hhbin/hhblits -i $tmpfile.fasta -d $db -n 3 -o $tmpfile.hhr -oa3m $tmpfile.a3m -cpu $cpu -neff $neff");
    &HHPaths::System("$hhscripts/a3m2mtx.pl $tmpfile.a3m $tmpfile.mtx");
    if ($outfile_mtx) {&HHPaths::System("cp $tmpfile.mtx $outfile_mtx");}
}
else
{
    &HHPaths::System("cp $outfile_mtx $tmpfile.mtx");
}

#### perform psipred prediction ####
print "Predicting secondary structure...\nPass1\n";
&HHPaths::System("$execdir/psipred $tmpfile.mtx $datadir/weights.dat $datadir/weights.dat2 $datadir/weights.dat3 > $tmpfile.ss");

print "Pass2 ...\n";
# note that the parameters is the same as addss.pl but slightly different
# from that in standard psipred
&HHPaths::System("$execdir/psipass2 $datadir/weights_p2.dat 1 0.98 1.09 $tmpfile.ss2 $tmpfile.ss > $tmpfile.horiz\n");

if (!$outfile_ss2)
{
    $outfile_ss2="$inbase.ss2";
    print "PSIPRED output written to $outfile_ss2\n";
}
&HHPaths::System("cp $tmpfile.ss2 $outfile_ss2");

#### preform solvpred prediction
if ($outfile_solv)
{
    print "Predicting solvent accessibility...\n";
    &HHPaths::System("$execdir/solvpred $tmpfile.mtx $datadir/weights_solv.dat > $tmpfile.solv\n");
    &HHPaths::System("cp $tmpfile.solv $outfile_solv");
}
