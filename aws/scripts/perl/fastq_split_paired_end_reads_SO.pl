#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
my ($e1, $e2);
GetOptions(
	   '1' => \$e1,
	   '2' => \$e2
          );
while(<>){
    chomp;
    print STDOUT "$_\/1\n" if $e1; 
    print STDOUT "$_\/2\n" if $e2;
    my $newline = <>; chomp($newline);
    print STDOUT substr($newline, 0, length($newline)/2)."\n" if $e1;
    print STDOUT substr($newline, length($newline)/2, length($newline)/2)."\n" if $e2;
    $newline = <>; chomp($newline);
    print STDOUT "$newline\/1\n" if $e1;
    print STDOUT "$newline\/2\n" if $e2;
    $newline = <>; chomp($newline);
    print STDOUT substr($newline, 0, length($newline)/2)."\n" if $e1;
    print STDOUT substr($newline, length($newline)/2, length($newline)/2)."\n" if $e2;
}
