#!/usr/bin/perl 

use Getopt::Long;
# take input of VCF file from mpileup 
# move DP4 values into samples col
# NB will only work on single-sample VCFs
# doesn't perform error checking

my @removes;
my @adds;

GetOptions(
           'removes|r=s' => \@removes,
	   'add|a=s' => \@adds
          );

print STDERR "@removes \n";

while(<>) {
  if(m/^\#\#/) {print $_;
		if ($_ =~ m/\#\#INFO\=\<ID\=DP4/gi) {
		  my $format = $_; $format =~ s/\#\#INFO/\#\#FORMAT/gi;
		  print $format;
		}
		next;}
  chomp;
  my @F = split("\t",$_);

  if(m/^\#CHROM/) {
    $sample = $ARGV;
    $sample =~ s/_part-\d+//gi;
    # only change if this is the only sample in the file
    $F[9] = $sample if $#F == 9; 
  }
  else {
    
    #REMOVE UNWANTED FEATURES FROM INDIVIDUALS
    my @keys = split(":",$F[8]);
    my @vals = split(":",$F[9]);
    for (my $i=0; $i<=$#keys; $i++) {
      foreach $remove (@removes) {
	if($keys[$i] eq $remove) {
	  splice(@keys, $i, 1);
	  splice(@vals, $i, 1);
	}
      }
     }
    
    #SHIFT WANTED VALUES (curr only DP4) FROM INFO -> SAMPLES
    $DP4_val = '0,0,0,0'; 
    $DP4_key = 'DP4';
    if($F[7] =~ m/DP4=(\d+,\d+,\d+,\d+);/gi) { 
      $DP4_val = $1; }
    #    $F[8] .= ":DP4"; $F[9] .=":".$DP4; 
    push(@keys, $DP4_key);
    push(@vals, $DP4_val);
    
    # REMAKE KEYS + VALS
    $F[8] = join(":",@keys); 
    $F[9] = join(":",@vals); 
  }

  print join("\t",@F)."\n";
}
