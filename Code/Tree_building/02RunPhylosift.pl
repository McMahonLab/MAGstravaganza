###############################################################################
# RunPhylosift.pl
# Copyright (c) 2014, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Takes an input directory of fasta files. For each file, runs phylosfit search
# and align steps to identify and align highly-conserved marker genes to the
# reference tree.
################################################################################
# Date of Last Revision: May 2, 2015 (update file description)
################################################################################

use strict;
 
# Bring in input directory
# usage statement if one or both arguments are missing.
my $usage  = "Command sequence: perl 02RunPhylosift.pl inputDirectory \n";
my $inputDirectory  = shift or die $usage;

# Read in the list of files in the input directory into an array
opendir(DIR, $inputDirectory) || die "Can't open directory: $inputDirectory\n";
my @dirList = grep !/^\.+/, readdir(DIR);

my $arrayLen = scalar @dirList;

for my $i (0 .. $arrayLen-1) {
  # Run the appropriate phylosift commands
  my $index = $i + 1;
  print "Running search step on @dirList[$i] ($index of $arrayLen) \n";
  system("/usr/local/bioinf/phylosift_v1.0.1/phylosift search --besthit $inputDirectory/@dirList[$i]");

  print "Running align step on @dirList[$i] ($index of $arrayLen) \n";
  system("/usr/local/bioinf/phylosift_v1.0.1/phylosift align --chunk_size 100000 --besthit $inputDirectory/@dirList[$i]");
}
