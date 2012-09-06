#!/usr/bin/perl
use strict;
use Getopt::Long;

my $ref;
my $samhead;
my $flank = 80;
my $no_jobs=30;
GetOptions(
           'jobs=s' => \$no_jobs,
	   'header=s' => \$samhead,
	   'flank=s' => \$flank,
          );


my $i = 0;
my %minmax;
my $split = -1;
my $base_count = 0;
if($samhead) {
  $split = parseHeader($samhead, $no_jobs);
}

while(<>) {
  my @F = split; 
  if ($F[0] =~ m/^\@SQ/ && !$samhead) {
    $F[2] =~ s/LN://;
    $base_count += $F[2];
  }
  if ($F[0] =~ m/^\@/) {
    next;
  }
  elsif ($split == -1) {
    if($base_count > 0) {$split = 100000;}
    else {  print STDERR "SPLIT = ".$split." = ".$base_count." / ".$no_jobs."\n";
	    $split = $base_count / $no_jobs;}
  }

#  print STDERR $_ if ($i < 10);

  my $posn = $F[3];
  my $chr = $F[2];
  
  my $rem = int($posn / $split);            #get whole divisor of posn by split 
  #print $chr\t$block as keys
#  print STDOUT $chr."\t".$rem."\t".$_;
  #print $chr.block\tpos as keys
  print STDOUT $chr.".".$rem."\t".$posn."\t".$_;
  #if start in lower flank of next one up... (i.e. might run into next block) 
  #also print in next key
  if(($rem +1)  * $split < $posn + $flank){
    print STDOUT $chr.".".($rem+1)."\t".$posn."\t".$_;
  }
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
