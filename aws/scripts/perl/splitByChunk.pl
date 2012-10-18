#!/usr/bin/perl
use strict;
use Getopt::Long;

my $ref;
my $samhead;
my $flank = 10000;
my $chunk = 1000000;
my $no_jobs=30;
GetOptions(
	   'flank=s' => \$flank,
	   'chunk=s' => \$chunk
          );

sub cleanNo {
  my $val = shift;
  $val =~ s/,//g;
  $val =~ s/['k','kb']/000/gi;
  $val =~ s/['m','mb']/000000/gi;
  return $val;
}

$flank = cleanNo($flank);
$chunk = cleanNo($chunk);
die("flank must be smaller than chunk") if $flank >= $chunk;
print STDERR "chunk = ".$chunk."\tflank = ".$flank."\n";

while(<>) {
  my ($chr,$pos,$val);
  if ($_ =~ m/^\w*\:*(\d)\.(\d+)\s+(.+)$/gi) {
    ($chr,$pos,$val) = ($1, $2, $3);
  }
  elsif ($_ =~ m/^(\w+)\t(\d+)\t(\S+)$/gi) {
    ($chr,$pos,$val) = ($1, $2, $3);
  }
  else {
  print STDERR "ERROR - input format not recognised: ".$_."\n";
  }

  my $rem = int($pos / $chunk);            #get whole divisor of posn by split 

  #print $chr.block\tpos as keys
  print STDOUT $chr.".".$rem."\t".$pos."\t".$val."\n";
  #if start in lower flank of next one up... (i.e. might run into next block) 
  #also print in next key
  if(($rem +1)  * $chunk <= $pos + $flank){
    print STDOUT $chr.".".($rem+1)."\t".$pos."\t".$val."\n";
  }
}
