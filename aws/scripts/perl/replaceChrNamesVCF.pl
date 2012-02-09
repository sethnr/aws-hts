#!/usr/bin/perl
use strict;
use Getopt::Long;

my $cmap = "./cmap.txt";


GetOptions(
	   'cmap=s' => \$cmap
	  );

open(CMAP,$cmap);
my @cmap;
while(<CMAP>) {
  my($name, $index) = split;
  print STDERR $index." => ".$name."\n";
  $cmap[$index] = $name;
}

while(<>) {
  if(m/^\#/gi) {
    print STDOUT $_;
  }
  else {
    my @F = split;
    $F[0] = $cmap[$F[0]];
    print STDOUT join("\t",@F)."\n";
  }
}

1;
