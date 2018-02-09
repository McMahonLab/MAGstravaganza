###############################################################################
# FastaSequenceMerger.pl
# Copyright (c) 2014, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Takes an input directory of fasta files. For each file, merges all of the
# sequences into a single contig. Writes to the specified output directory.
################################################################################
# Date of Last Revision: November 12, 2014
################################################################################

use strict;

# Bring in input and output directories
# usage statement if one or both arguments are missing.
my $usage  = "Command sequence: perl 01FastaSequenceMerger.pl inputDirectory outputDirectory \n";
my $inputDirectory  = shift or die $usage;
my $outputDirectory  = shift or die $usage;

# Read in the list of files in the input directory into an array
opendir(DIR, $inputDirectory) || die "Can't open directory: $inputDirectory\n";
my @dirList = grep !/^\.+/, readdir(DIR);

foreach (@dirList) {
  # Establish file for writing
  open(inputFile,'<', "$inputDirectory/$_");
  open(outputFile,'>', "$outputDirectory/$_");
  select(outputFile);

  # Write out header information for the fasta file
  print ">$_ \n";
  
  # Read in each line and write if appropriate
  while(<inputFile>) {
    chomp($_);
    if ($_ !~ />(.+)/) {
      if ($_ !~ /^\s*$/) {
	print "$_";
      }
    }
  } 
}

