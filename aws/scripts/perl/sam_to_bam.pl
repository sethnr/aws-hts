#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "sam_to_bam_tmp";

my $ref = "AgamP3.fa";
my $bowtie;
my $samhead;
my $cliptag;


GetOptions(
           'fasta|ref=s' => \$ref,
	   'tmp=s' => \$tmpfile,
	   'bowtie' => \$bowtie,
	   'header=s' => \$samhead,
	   'cliptag' => \$cliptag
          );


open(TMP,">".$tmpfile.".sam");

print STDERR "dumping input into file\n";
my $i = 0;
while(<>) {
  my @F = split;
  next if m/FAKE/g;
  if($bowtie) {
    print TMP join("\t",$F[9],$F[3],$F[0],$F[2],@F[4..7])."\n";
  }
  elsif($cliptag) {
    if($F[0] =~ m/^\@/) {print TMP $_;}
    else {
      print TMP join("\t",@F[2..$#F])."\n";
}
  }
  else {
    print TMP "$_";
  }
    $i++;
}
close(TMP);
print STDERR $i." lines to process\n";

unless($i > 0) {
  warn("WARN: input file empty, check inputs for errors!");
  exit 0;
}

my $command;

if($samhead) {
  $command = "cat ".$samhead." ".$tmpfile.".sam > ".$tmpfile.".sam.h"; 
  print STDERR "\t".$command."\n" ;
  runcheck($command);
  $command = "mv ".$tmpfile.".sam.h ".$tmpfile.".sam";   
}

$command = "head ".$tmpfile.".sam";
print STDERR `$command`;


#print STDERR "converting SAM to BAM\n";
$command = "./samtools view -Sb ".$tmpfile.".sam > ".$tmpfile.".bam";
 print STDERR "\t".$command."\n";
runcheck($command);

print STDERR "sorting bam files\n";
$command = "./samtools sort ./".$tmpfile.".bam ".$tmpfile.".s";
print STDERR "\t".$command."\n";
runcheck($command);


$command = "cat ./".$tmpfile.".s.bam";
runcheck($command);  #print to stdout

exit 0;

sub runcheck() {
  my $command = shift;
  system($command) == 0
    or die "command:\n".$command."\n\tfailed: $? \n" .  `find .`; 
}
