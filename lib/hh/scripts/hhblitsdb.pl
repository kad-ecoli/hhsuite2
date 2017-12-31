#!/usr/bin/env perl
#
# hhblits.pl
# Creates HH-suite database files from A3M and HHM/HMMER-formatted files 
# Usage: Usage: perl hhblitsdb.pl -o <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>] [more_options]
#
#     HHsuite version 2.0.16 (Sept 2012)
#
#     Reference: 
#     Remmert M., Biegert A., Hauser A., and Soding J.
#     HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
#     Nat. Methods, epub Dec 25, doi: 10.1038/NMETH.1818 (2011).

#     (C) Johannes Soeding, 2012

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.

#     We are very grateful for bug reports! Please contact us at soeding@genzentrum.lmu.de

use lib $ENV{"HHLIB"}."/scripts";
use HHPaths;   # config file with path variables for nr, blast, psipred, pdb, dssp etc.
use strict;
#use File::Glob 'bsd_glob'; # splits patterns delimited by spaces into multiple patterns and applies them using OR

$|= 1; # Activate autoflushing on STDOUT

# Default values:
our $v=2;             # verbose mode
my $a_if_append = ""; # do not append by default (default: create new db)
my $remove = 0;       # do not remove by default (default: create new db)
my $hhmext = "hhm";   # default HHM-file extension
my $csext = "seq219";   # default HHM-file extension
my $cpu = 8;

# Variable declarations
my $line;
my $command;
my $a3mdir = "";     # name of input A3M directory
my $hhmdir = "";     # name of input HHM/HMM directory
my $csdir  = "";      # name of input cs directory
my $a3mfile = "";     # name of packed ouput A3M file 
my $hhmfile = "";     # name of packed ouput HHM file
my $csfile = "";      # name of cs sequence db file
my $dbname = "";      # output db name
my $logfile = "/dev/null"; # log file 
my $file;
my $numcsfiles= 0;
my $num_chars = 0;
my $numa3mfiles=0;
my $numhhmfiles=0;
my $fileglob="";
my $help="
hhblitsdb.pl from HHsuite $VERSION  
Builds HH-suite database from a3m formatted MSAs and/or from HMMs (-o).
MSAs and HMMs can also be added (-a) to or removed (-r) from an existing database. 

Usage: hhblitsdb.pl -o|-a|-r <db_name> [-ia3m <a3m_dir>] [-ihhm <hhm_dir>] [-ics <cs_dir>]... 

With option -o, the following HH-suite database files can be generated:
 <db_name>.cs219              column-state sequences, one for each MSA/HMM (for prefilter)
 <db_name>.cs219.sizes        number of sequences and characters in <db_name>.cs219
 <db_name>_a3m_db             packed file containing A3M alignments read from <a3m_dir>
 <db_name>_a3m_db.index       index file for packed A3M file
 <db_name>_a3m.db.index.sizes number of lines in <db_name>_a3m_db.index
 <db_name>_hhm_db             packed file containing HHM-formatted HMMs read from <hhm_dir>
 <db_name>_hhm_db.index       index file for packed HHM file
 <db_name>_hhm_db.index.sizes number of lines in <db_name>_hhm_db.index

Options:
 -o <db_name>    create database with this name
 -a <db_name>    append files to database with this name
 -r <db_name>    remove files from database with this name
 -ia3m <a3m_dir> input directory (or glob of directories) with A3M-formatted files
                 These files MUST have extension 'a3m'.
 -ihhm <hhm_dir> input directory (or glob of directories) with HHM (or HMMER) files 
                 These files MUST have extension 'hhm' (HHsuite) or 'hmm' (HMMER3). 
 -ics  <cs_dir>  input directory (or glob of directories) with column state sequences
 -log <logfile>  log file recording stderr stream of cstranslate and hhmake commands
 -csext <ext>    extension of column state sequences (default: $csext)
 -hmm            use HMMER-formatted files. These MUST have extension hmm
                 (WARNING! HMMER format results in decreased performance over HHM format)
 -v [1-3]        verbose mode (default: $v)
 -cpu <int>      number of threads to generate cs219 and hhm files (default = $cpu)
 -f 'file_glob'  string with list of glob expression of files to remove
 
Example 1: only -ia3m given; cs sequences and hhm files are generated from a3m files
   perl hhblitsdb.pl -o databases/mydb -ia3m mydb/a3ms/ 

Example 2: only -ihhm given; cs sequences are generated from hhm files, but no a3m db file 
   perl hhblitsdb.pl -o databases/mydb -ihhm mydb/hhms/ 

Example 3: -ia3m and -ihhm given; cs sequences are generated from a3m files
   perl hhblitsdb.pl -o databases/mydb -ia3m mydb/a3ms/ -ihhm mydb/hhms/   

Example 4: -ics, -ia3m, and -ihhm given; all db files are created 
   perl hhblitsdb.pl -o databases/mydb -ia3m mydb/a3ms/ -ihhm mydb/hhms/ -ics mydb/cs/  

Example 5: using glob expression to specify files (note the singe quotes)
   perl hhblitsdb.pl -o databases/mydb -ihhm 'mydbs*/hhms/*.hhm'  

Example 6: add files to database; cs sequences and hhm files are generated from a3m files
   perl hhblitsdb.pl -a databases/mydb -ia3m 'mydbs/a3ms/g1a*.a3m'  

Example 7: remove files from database
   perl hhblitsdb.pl -r databases/mydb -f 'mydbs/a3ms/g1a*.* mydbs2/'  
\n";


###############################################################################################
# Processing command line input
###############################################################################################

if (@ARGV<1) {die ($help);}

for (my $i=0; $i<@ARGV; $i++) {
    if ($ARGV[$i] eq "-ics") {
	if (++$i<@ARGV) {
	    $csdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ics option!\n");
	}
    } elsif ($ARGV[$i] eq "-ia3m") {
	if (++$i<@ARGV) {
	    $a3mdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ia3m option!\n");
	}
    } elsif ($ARGV[$i] eq "-ihhm") {
	if (++$i<@ARGV) {
	    $hhmdir=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing directory after -ihhm option!\n");
	}
    } elsif ($ARGV[$i] eq "-log") {
	if (++$i<@ARGV) {
	    $logfile=$ARGV[$i];
	    unlink $logfile;
	} else {
	    die ("$help\n\nERROR! Missing filename after -log option!\n");
	}
    } elsif ($ARGV[$i] eq "-csext") {
	if (++$i<@ARGV) {
	    $csext=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing extension after -csext option!\n");
	}
    } elsif ($ARGV[$i] eq "-hmm") {
	$hhmext="hmm";
	print("\nWARNING! HMMER format results in decreased performance over HHM format. We recommend to generate hhm files directly from multiple sequence alignments using hmake.\n");
    } elsif ($ARGV[$i] eq "-v") {
	if (++$i<@ARGV) {
	    $v=$ARGV[$i];
	} else {
	    $v = 2;
	}
    } elsif ($ARGV[$i] eq "-cpu") {
	if (++$i<@ARGV) {
	    $cpu=$ARGV[$i];
	}
    } elsif ($ARGV[$i] eq "-f") {
	if (++$i<@ARGV) {
	    $fileglob=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing expression after -f option!\n");
	}
   } elsif ($ARGV[$i] eq "-r") {
	if (++$i<@ARGV) {
	    if ($dbname!="") {die("$help\n\nERROR! options -o and -r not compatible!\n");}
	    $dbname=$ARGV[$i];
	    $remove=1;
	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } elsif ($ARGV[$i] eq "-a") {
	if (++$i<@ARGV) {
	    if ($remove==1) {die("$help\n\nERROR! options -r and -a not compatible!\n");}
	    if ($dbname!="") {die("$help\n\nERROR! options -o and -a not compatible!\n");}
	    $dbname=$ARGV[$i];
	    $a_if_append="a";
	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } elsif ($ARGV[$i] eq "-o") {
	if (++$i<@ARGV) {
	    if ($remove==1) {die("$help\n\nERROR! options -r and -o not compatible!\n");}
	    if ($a_if_append) {die("$help\n\nERROR! options -a and -o not compatible!\n");}
	    $dbname=$ARGV[$i];
	} else {
	    die ("$help\n\nERROR! Missing filename after -o option!\n");
	}
    } else {
	if ($dbname="") {
	    $dbname=$ARGV[$i];
	} else {
	    print "WARNING! Unknown option $ARGV[$i]!\n";
	}
    }
}

# Check input
if (!$dbname) {print($help); die("ERROR! Name of database is missing! Use -o <db_name>\n");}
$a3mdir=~s/\/$//; 
$a3mfile = $dbname."_a3m_db"; 
$hhmfile = $dbname."_hhm_db";
$csfile = $dbname.".cs219";
if ($hhmdir) {$hhmdir=~s/\/$//;} 

if ($a_if_append eq "" && $remove==0) {
    unlink $csfile, $a3mfile, $a3mfile.".index", $hhmfile, $hhmfile.".index"; 
}

if ($a3mdir eq "" && $hhmdir eq "" && $csdir eq "" && $remove==0) {
    print($help); print "ERROR! At least one input directory must be given!\n"; exit(1);
}

# If $csdir is simple directory instead of glob expression, turn it into glob expression
if ($csdir) {
    if ($csdir !~ /\*/ && $csdir !~ /\?/ && $csdir !~ / /) { 
	$csdir .= "/*.".$csext;
    }
}
if ($a3mdir) {
    if ($a3mdir !~ /\*/ && $a3mdir !~ /\?/ && $a3mdir !~ / /) { 
	$a3mdir .= "/*.a3m";
    }
}    
if ($hhmdir) {
    if ($hhmdir !~ /\*/ && $hhmdir !~ /\?/ && $hhmdir !~ / /) { 
	$hhmdir .= "/*.".$hhmext;
    }
}


# If in append mode, initialize size counters with present sizes
if ($a_if_append || $remove==1) {
    open (IN, "<$a3mfile.index.sizes") || die("Error: can't open $a3mfile.index.sizes: $!");
    $line = <IN>;
    close IN;
    $line =~ /^(\S*)/;
    $numa3mfiles = $1;

    open (IN, "<$hhmfile.index.sizes") || die("Error: can't open $hhmfile.index.sizes: $!");
    $line = <IN>;
    close IN;
    $line =~ /^(\S*)/;
    $numhhmfiles = $1;

    open (IN, "<$csfile.sizes") || die("Error: can't open $csfile.sizes: $!");
    $line = <IN>;
    close IN;
    $line =~ /(\S*)\s+(\S*)/;
    $numcsfiles = $1;
    $num_chars = $2;

    printf("Current number of a3m files in db:    %i\n",$numa3mfiles);
    printf("Current number of $hhmext files in db:    %i\n",$numhhmfiles);
    printf("Current number of $csext files in db: %i\n\n",$numcsfiles);

} else {
    $numa3mfiles = 0;
    $numhhmfiles = 0;
    $numcsfiles = 0;
    $num_chars = 0;
}

# Create tmp directory (plus path, if necessary)
my $tmpdir="/tmp/$ENV{USER}/$$";  # directory where all temporary files are written: /tmp/UID/PID
my $suffix=$tmpdir;
while ($suffix=~s/^\/[^\/]+//) {
    $tmpdir=~/(.*)$suffix/;
    if (!-d $1) {mkdir($1,0777);}
} 
unlink glob("$tmpdir/*"); # clean up directory if it already exists
unlink $logfile;
    


##############################################################################################
# Remove files?
##############################################################################################
if ($remove==1) {

    printf("Removing files from indices...\n");
    
    # Read numbers of sequences and characters in csfile
    open (IN, "<$csfile.sizes");
    $line = <IN>;
    close IN;
    $line =~ /(\S*)\s+(\S*)/;
    $numcsfiles = $1;
    $num_chars = $2;

    # Remove names from a3m and hhm index files
    my $files = " ".join(" ", glob($fileglob));
    $files =~ s/\S*\///g;
    &HHPaths::System("ffindex_modify -su $dbname"."_a3m_db.index ".$files); 
    &HHPaths::System("ffindex_modify -su $dbname"."_hhm_db.index ".$files); 
    
    # Adjust number of files in $a3mfile.index.sizes
    $numa3mfiles = 0;
    open (IN, "<$a3mfile.index") || die("Error: can't open $a3mfile.index: $!");
    while(<IN>) {$numa3mfiles++;}
    close IN;
    open (OUT, ">$a3mfile.index.sizes") || die("Error: can't open $a3mfile.index.sizes: $!");
    printf(OUT "%i\n",$numa3mfiles);
    close(OUT);

    # Adjust number of files in $hhmfile.index.sizes
    $numhhmfiles = 0;
    open (IN, "<$hhmfile.index") || die("Error: can't open $hhmfile.index: $!");
    while(<IN>) {$numhhmfiles++;}
    close IN;
    open (OUT, ">$hhmfile.index.sizes") || die("Error: can't open $hhmfile.index.sizes: $!");
    printf(OUT "%i\n",$numa3mfiles);
    close(OUT);


    # Remove sequences of globbed files from cs file
    my 	$skipseq=0;
    $numcsfiles = 0;
    $num_chars = 0;
    open (IN,  "<$csfile");
    open (OUT, ">$csfile".".tmp");
    foreach my $line (<IN>) {
	if ($line =~ /^>(\w*)/) {
	    my $name = $1;
	    if ($files =~ / $name\./) { # found name in list of globbed file names?
		$skipseq=1;
	    } else {
		$skipseq=0;
		printf(OUT "%s",$line); 
		$numcsfiles++;
	    }
	} else {
	    if (!$skipseq) {
		printf(OUT "%s",$line); 
		$num_chars += length($line);
	    }
	}
    }
    close(OUT);	
    close(IN);	
    unlink($csfile);
    &HHPaths::System("mv $csfile".".tmp ".$csfile); 

    # Adjust $csfile.sizes
    open (OUT, ">$csfile.sizes");
    print OUT "$numcsfiles $num_chars\n";
    close OUT;

} else {

##############################################################################################
# Generate new db or append to old
##############################################################################################


##############################################################################################
# Generate column-state database file
##############################################################################################

    # Generate column-state sequences in $tmpdir if no -ics directory given
    if (!$csdir) 
    {
	my $x = 0.3;    # parameters for cstranslate
	my $c = 4;      # parameters for cstranslate
	
	if ($a3mdir) {
	    print("Generating seq219 files in $tmpdir/ from a3m files $a3mdir\n\n");
	    $command = "$hhbin/cstranslate -i \$file -o $tmpdir/\$base.seq219 -D $context_lib -A $cs_lib -x $x -c $c 1>>$logfile 2>>$logfile";
	    &HHPaths::System("$hhscripts/multithread.pl '".$a3mdir."' '$command' -cpu $cpu");
	    
	} elsif ($hhmdir) {
	    
	    if ($hhmext eq "hmm") {
		print("\nGenerating prf profile files in $tmpdir/ from hmm files $hhmdir/\n\n");
		$command = "$hhscripts/create_profile_from_hmmer.pl -i \$file -o $tmpdir/\$base.prf 1>/dev/null 2>>$logfile";
		&HHPaths::System("$hhscripts/multithread.pl '".$hhmdir."' '$command' -cpu $cpu");
	    } else { # $hhmext eq "hhm"
		print("\nGenerating prf profile files in $tmpdir/ from hhm files $hhmdir/\n\n");
		$command = "$hhscripts/create_profile_from_hhm.pl -i \$file -o $tmpdir/\$base.prf 1>/dev/null 2>>$logfile";
		&HHPaths::System("$hhscripts/multithread.pl '".$hhmdir."' '$command' -cpu $cpu");
	    }
	    
	    print("\nGenerating seq219 files in $tmpdir/ from prf files in $tmpdir/\n\n");
	    if ($hhmext eq "hmm") {
		$command = "$hhbin/cstranslate -i \$file -o \$name.seq219 -A $cs_lib 1>>$logfile 2>>$logfile";	    
	    } else { # $hhmext eq "hhm"
		$command = "$hhbin/cstranslate -i \$file -o \$name.seq219 -A $cs_lib -D $context_lib -x $x -c $c 1>>$logfile 2>>$logfile";
	    }
	    &HHPaths::System("$hhscripts/multithread.pl '".$tmpdir."/*.prf' '$command' -cpu $cpu");
	}
	
	$csdir = $tmpdir."/*.$csext";
    }
    
    
    # Write columns state sequences into cs database file, 
    # replace names in cs sequences with filenames: ">name+description" => ">filename"
    if ($a_if_append) {
	open (OUT, ">>$csfile");
    } else {
	open (OUT, ">$csfile");
    }
    foreach my $seq219file (glob($csdir)) { 
	open (IN, "<$seq219file");
	my @lines = <IN>;
	close(IN);
	$seq219file =~ s/.*?([^\/]*)\.$csext\s*/$1/ or die ("Error: $seq219file does not have the extension $csext!?\n");
	foreach my $line (@lines) {
	    if ($line =~ /^>/) {
		$line = ">".$seq219file."\n";
		$numcsfiles++;
	    } else {
		$num_chars += length($line);
	    }
	    printf(OUT "%s",$line); 	
	}
    } 
    close(OUT);	
    
    open (OUT, ">$csfile.sizes");
    print OUT "$numcsfiles $num_chars\n";
    close OUT;
    
    
##############################################################################################
# Generate hhm files with hhmake from a3m files if no -ihhm directory given
##############################################################################################
    
    if (!$hhmdir) 
    {
	if ($a3mdir) {
	    print("\nGenerating hhm files in $tmpdir/ from a3m files $a3mdir/\n\n");
	    $command = "hhmake -i \$file -o $tmpdir/\$base.hhm  1>/dev/null 2>>$logfile";
	    &HHPaths::System("$hhscripts/multithread.pl '".$a3mdir."' '$command' -cpu $cpu");	
	    $hhmdir = $tmpdir."/*.$hhmext";;
	    $numhhmfiles += scalar(glob("$hhmdir"));
	}
    }
    
    
##############################################################################################
# Generate packed A3M and HMM files and index files
##############################################################################################
    
    # Generate packed A3M file and index file?
    if ($a3mdir ne "") {
	print "Creating packed A3M database file $a3mfile ...\n";
	
	open (OUT, ">$tmpdir/a3m.filelist");
	my @files = glob("$a3mdir");
	$numa3mfiles += scalar(@files);
	foreach $file (@files) {
	    print OUT "$file\n";
	}
	close OUT;
	
	# Build packed file (concatenated with '\0' as delimiters) and index file from files in file list
	# The ffindex binaries are contained in <install_dir>/lib/ffindex/bin/
	$command = "ffindex_build -".$a_if_append."s -f $tmpdir/a3m.filelist $a3mfile $a3mfile.index";
	&HHPaths::System($command);
	
	open (OUT, ">$a3mfile.index.sizes");
	print OUT "$numa3mfiles\n";
	close OUT;
    } 
    
    # Generate packed HHMM file and index file?
    if ($hhmdir ne "") {
	print "Creating packed HHM database file $hhmfile ...\n";
	
	open (OUT, ">$tmpdir/hhm.filelist");
	my @files = glob("$hhmdir");
	$numhhmfiles += scalar(@files);
	foreach $file (@files) {
	    print OUT "$file\n";
	}
	close OUT;
	
	# Build packed file (concatenated with '\0' as delimiters) and index file from files in file list
	# The ffindex binaries are contained in <install_dir>/lib/ffindex/bin/
	$command = "ffindex_build -".$a_if_append."s -f $tmpdir/hhm.filelist $hhmfile $hhmfile.index";
	&HHPaths::System($command);
	
	open (OUT, ">$hhmfile.index.sizes");
	print OUT "$numhhmfiles\n";
	close OUT;
    }
} # end if $remove==0

print("\n");
printf("New number of a3m files in db:    %i\n",$numa3mfiles);
printf("New number of $hhmext files in db:    %i\n",$numhhmfiles);
printf("New number of $csext files in db: %i\n\n",$numcsfiles);

my $err=0;
if ($numa3mfiles && $numhhmfiles && $numa3mfiles != $numhhmfiles) {
    print("**************************************************************************
WARNING: Number of a3m files not equal to number of $hhmext files\n"); $err=1;
}
if ($numcsfiles && $numhhmfiles && $numcsfiles != $numhhmfiles) {
    print("**************************************************************************
WARNING: Number of $csext files not equal to number of $hhmext files\n"); $err=1;
}
if ($numcsfiles && $numa3mfiles && $numcsfiles != $numa3mfiles) {
    print("**************************************************************************
WARNING: Number of $csext files not equal to number of a3m files\n"); $err=1;
}

if ($err==1) {print("**************************************************************************
$tmpdir will not be removed to check for missing files
**************************************************************************\n");}
elsif ($v<3) {
    $command = "rm -rf $tmpdir";
#    &System($command);
}
wait;
exit;

