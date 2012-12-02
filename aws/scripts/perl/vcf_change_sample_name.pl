#!/usr/bin/perl

# take input of VCF file from mpileup 
# move DP4 values into samples col
# NB will only work on single-sample VCFs
# doesn't perform error checking

while(<>) {
  unless(m/^\#CHROM/) {print $_; next;}
  chomp;
  my @F = split("\t",$_);
  $sample = $ARGV;
  $sample =~ s/_part-\d+//gi;
  # only change if this is the only sample in the file
  $F[9] = $sample if $#F == 9; 
  print join("\t",@F,"\n");
}
