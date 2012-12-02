#!/usr/bin/perl
use strict;
use Getopt::Long;

#this version is optimised for large sequence files

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

print STDERR "parsing feature blocks\n";
my $no_lines = `wc -l $feat_file`; my $i=0;
my $milestone = $no_lines/10;
while(<FEAT>) {
  $i ++;
  if ($i >= $milestone) {print STDERR "parsed ".$milestone." blocks\n"; $milestone += $no_lines/10;} 
  my ($feat, $chr, $start, $end) = split;

  push(@{$feats{$chr}{"START"}},$start);
  push(@{$feats{$chr}{"END"}},$end);
  push(@{$feats{$chr}{"FEAT"}},$feat);
  $feats{$chr}{"COUNT"} = 0 if(!$feats{$chr}{"COUNT"});
  $feats{$chr}{"COUNT"}++;
}
print STDERR "feature blocks parsed\n";

my $in_feat = "";
my $feat_end;
my $in_chr;

my $i=0;
while(<>) {
  $i++;
  print STDERR $i." lines parsed\n" if ($i % 10000 == 0);
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

  # if in gene:
  if ($in_feat) {
    # check haven't left gene:
    if($pos >= $feat_end) {
      $feat_end = undef;
      $in_feat = undef;
    }
    elsif ($chr != $in_chr) {
      $in_chr = undef;    
      $in_feat = undef;
    }
    else {
      print STDOUT $in_feat."\t".$chr.".".$pos."\t".$val."\n";
    }
  }
  for (my $i=0; $i <= $feats{$chr}{"COUNT"}; $i++) {
    if ($pos >= $feats{$chr}{"START"}->[$i] && $pos <= $feats{$chr}{"END"}->[$i]) {
      $in_feat = $feats{$chr}{"FEAT"}->[$i];
      $feat_end = $feats{$chr}{"END"}->[$i];
      print STDOUT $in_feat."\t".$chr.".".$pos."\t".$val."\n";
    }

  }
#  if($feats{$chr}{$pos}) {
#    foreach my $feat (@{$feats{$chr}{$pos}}) {
#      print STDOUT $feat."\t".$chr.".".$pos."\t".$val."\n";
#    }
#  }
}
