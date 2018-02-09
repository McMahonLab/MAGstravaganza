###############################################################################
# CreateAlignmentFile.pl
# Copyright (c) 2014, Joshua J Hamilton and Katherine D McMahon
# Affiliation: Department of Bacteriology
#              University of Wisconsin-Madison, Madison, Wisconsin, USA
# URL: http://http://mcmahonlab.wisc.edu/
# All rights reserved.
################################################################################
# Resides in the phylosift folder. Uses FastTree to construct a phylogenetic
# tree. Command takes as arguments names for the alignment and tree files.
################################################################################
# Date of Last Revision: May 2, 2015 (update file description)
################################################################################

use strict;
 
# Bring in input directory
# usage statement if one or both arguments are missing.
my $usage  = "Command sequence: perl 03CreateAlignmentFile.pl alignmentFileName.aln \n";
my $alignmentFileName  = shift or die $usage;

# Create the output file
open(outputFile,'>', "$alignmentFileName");
#select(outputFile);

# Read in the list of files in the input directory into an array
opendir(DIR, "/usr/local/bioinf/phylosift_v1.0.1/PS_temp") || die "Can't open directory: PS_temp\n";

my @dirList = grep {-d "/usr/local/bioinf/phylosift_v1.0.1/PS_temp/$_" && ! /^\.{1,2}$/} readdir(DIR);
my $arrayLen = scalar @dirList;

for my $i (0 .. $arrayLen-1) {
  my $index = $i + 1;
  select(STDOUT);
  print "Processing $dirList[$i] ($index of $arrayLen) \n";

  # Do the actual processing
  open(inputFile,'<', "/usr/local/bioinf/phylosift_v1.0.1/PS_temp/$dirList[$i]/alignDir/concat.updated.1.fasta");
  select(outputFile);

  # Rewrite the header and leave the rest alone
  print ">$dirList[$i]\n";
  while(<inputFile>) {
    if ($_ !~ />(.+)/) {
      print "$_";
    }
  }
}
