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
  print STDERR "cleaning ".$val;
  $val =~ s/,//g;
  if ($val =~ m/^([\d,\.]+)['k','kb']/gi) {
    $val = $1 * 1000;
  }
  elsif($val =~ m/^([\d,\.]+)['m','mb']/gi) {
    $val = $1 * 1000000;
  }
  print STDERR " -> ".$val."\n";
  return $val;
}

sub parseHeader {
  my $head_file = shift;
  my $no_jobs = shift;
  my $split = -1;
  my $base_count = 0;
  open(HEADER, "<$head_file");
  while(<HEADER>) {
    my @F = split; 
    if ($F[0] =~ m/^\@SQ/) {
      $F[2] =~ s/LN://;
      $base_count += $F[2];
    }
    if ($F[0] =~ m/^\@/) {
      next;
    }
    elsif ($split == -1) {
    }
  }
  $split = $base_count / $no_jobs;
  print STDERR "SPLIT = ".$split." = ".$base_count." / ".$no_jobs."\n";
  return $split;
}

sub getKey {
  my $rem = shift;
  my $chr = shift;
  return $chr.':'.($rem * $chunk)."-".(($rem+1)*$chunk);
}


$flank = cleanNo($flank);
$chunk = cleanNo($chunk);
die("flank must be smaller than chunk") if $flank >= $chunk;
print STDERR "chunk = ".$chunk."\tflank = ".$flank."\n";

while(<>) {
  my @F = split; 
  if ($F[0] =~ m/^\@/) {
    next;
  }


  my ($chr,$pos) = ($F[2], $F[3]);
  next if $chr eq '*';
  my $rem = int($pos / $chunk);            #get whole divisor of posn by split 

#  print STDOUT $chr.".".$rem."\t".$_;
  print STDOUT getKey($rem,$chr)."\t".$_;
  #if start in lower flank of next one up... (i.e. might run into next block) 
  #also print in next key
  if(($rem +1)  * $chunk <= $pos + $flank){
    print STDOUT getKey(($rem+1),$chr)."\t".$_;
  }
}


