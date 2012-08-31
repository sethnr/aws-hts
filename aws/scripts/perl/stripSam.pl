#!/usr/bin/perl
use strict;
use Getopt::Long;

my $samhead;
my $aligned = undef;
my $unaligned = undef;
my $cliptag = undef;
GetOptions(
	   'header=s' => \$samhead,
	   'aligned' => \$aligned,
	   'unaligned' => \$unaligned,
	   'cliptag' => \$cliptag,
          );

if((!$aligned) && (!$unaligned)) {
  $aligned = 1;
  $unaligned = 1;
}

my $header = parseHeader($samhead) if $samhead;

while(<>) {
  my @F = split;
  next if (m/^\@/);
  next if (m/FAKE/g);
  next if ((!$unaligned) && $F[4] eq '*');
  next if ((!$aligned) && $F[4] ne '*');

  print $header if $header; 
  $header = undef;

  print join("\t",@F[2..$#F])."\n" if $cliptag; 
  print join("\t",@F)."\n" unless $cliptag; 
}



sub parseHeader {
  my $head_file = shift;
  open(HEADER, "<$head_file");
  while(<HEADER>) {
    $header .= $_;
  }
  return $header;
}
