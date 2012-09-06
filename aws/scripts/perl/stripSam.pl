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

my $chrpos = 4;
$chrpos = 2 if (!$cliptag); 


while(<>) {
  my @F = split;
  next if (m/^\@/);
  next if (m/FAKE/g);
  next if ((!$unaligned) && $F[$chrpos] eq '*');
  next if ((!$aligned) && $F[$chrpos] ne '*');

  die "possible malformed sam-tag file\n".$_."\n" if ($F[2] ne 'r' && $cliptag);
  die "possible malformed sam file\n".$_."\n" if ($F[0] ne 'r' && !$cliptag);

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
