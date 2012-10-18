#!/usr/bin/perl

while(<>) {
  $_ =~ s/^(\d+)\.(\d+)/$1/gi;
  print $_;
}
