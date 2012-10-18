#!/usr/bin/perl
use strict;
use Getopt::Long;

#this version is optimised for large feature files


my %feats;
my $feat_file;
my $flank = 0;

GetOptions(
	   'flank=s' => \$flank,
	   'features|feats=s' => \$feat_file
          );

sub cleanNo {
  my $val = shift;
  $val =~ s/,//g;
  $val =~ s/['k','kb']/000/gi;
  $val =~ s/['m','mb']/000000/gi;
  return $val;
}

$flank = cleanNo($flank);



open(FEAT,"<$feat_file") or die "could not open feat_file: ".$feat_file."\n";


# my %blocks;
# my $last_chr = -1;
# my $last_pos = -1;

my %vals;
my $i=0;
while(<>) {
  $i++;
  print STDERR $i." lines parsed\n" if ($i % 100000 == 0);
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
  $vals{$chr}{$pos} = $val;
  
#  if($chr != $last_chr || $pos != $last_pos+1)
#  $last_pos = $pos; 
#  $last_chr = $chr; 
}
print STDERR "stdin parsed: ".$i." lines\n";


print STDERR "parsing feature blocks\n";
my $no_lines = `wc -l $feat_file`; my $i=0;
my $milestone = $no_lines/10;
while(<FEAT>) {
  $i ++;
  if ($i >= $milestone) {print STDERR "parsed ".$milestone." blocks\n"; $milestone += $no_lines/10;} 
  my ($feat, $chr, $start, $end) = split;

  for (my $base = $start; $base <= $end; $base++) {
    if($vals{$chr}{$base}) { print STDOUT $feat."\t".$chr.".".$base."\t".$vals{$chr}{$base}."\n";}
  }
}
print STDERR "feature blocks parsed\n";

