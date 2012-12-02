#!/usr/bin/perl -i 

# take input of VCF file from mpileup 
# move DP4 values into samples col
# NB will only work on single-sample VCFs
# doesn't perform error checking

while(<>) {
  if(m/^\#/) {print $_; next;}
  chomp;
  my @F = split("\t",$_);
  $DP4 = '0,0,0,0'; 
  if($F[7] =~ m/DP4=(\d+,\d+,\d+,\d+);/gi) { 
    $DP4 = $1; }
 
  $F[8] .= ":DP4"; $F[9] .=":".$DP4; 
  print join("\t",@F,"\n");
}
