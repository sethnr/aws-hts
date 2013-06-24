#!/usr/bin/perl -w

open(FF,">fwd.fastq");
open(FR,">rev.fastq");


while(<STDIN>) {
  ($name, $rf, $qf, $rr, $qr) = split;
  print FF '@'.$name."\n".$rf."\n+\n".$qf."\n\n";
  print FR '@'.$name."\n".$rr."\n+\n".$qr."\n\n";
}
close(FF);
close(FR);

$btindex = shift(@ARGV);
$ENV{'PATH'} = '.:./bin/:'.$ENV{'PATH'};
$command = "bowtie2 ".join(" ",@ARGV)." -x ".$btindex." -1 fwd.fastq -2 rev.fastq";
print STDERR $command."\n";
system($command);

exit 0;
