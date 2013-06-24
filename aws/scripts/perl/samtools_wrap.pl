#!/usr/bin/perl
use strict;

my @ARGV = shift;

my $tmpfile = "samtools_tmp.sam";

open(TMP,">".$tmpfile);

print STDERR "dumping input into file\n";
my $i = 0;
while(<>) {
  $i++;
  print TMP "$_";
}
close(TMP);
print STDERR $i." lines to process\n";

unless($i > 0) {
  warn("WARN: input file empty, check inputs for errors!");
  exit 0;
}

runcheck("samtools @ARGV $tmpfile");


exit 0;

sub runcheck() {
  my $command = shift;
  system($command) == 0
    or die "command:\n".$command."\n\tfailed: $? \n" .  `find .`; 
}
