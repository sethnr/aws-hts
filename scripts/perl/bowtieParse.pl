#!/usr/bin/perl

print STDERR "parsing file @ARGV\n";

while(<>) {
my @F = split;
print STDOUT join("\t",$F[9],$F[3],$F[0],$F[2],@F[4..7])."\n";
}
