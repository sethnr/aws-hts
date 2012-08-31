#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "sam_to_vcf_tmp";

my $ref = "AgamP3.fa";
my $bowtie;
my $samhead;
my $cliptag;


GetOptions(
           'fasta|ref=s' => \$ref,
	   'tmp=s' => \$tmpfile,
	   'bowtie' => \$bowtie,
	   'header=s' => \$samhead,
	   'cliptag=s' => \$cliptag
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
    print TMP join("\t",@F[2..$#F])."\n";
  }
  else {
    print TMP "$_";
  }
    $i++;
}
close(TMP);
print STDERR $i." lines to process\n";

my $command;

$command = "cat ".$samhead." ".$tmpfile.".sam > ".$tmpfile.".sam.h";
print STDERR "\t".$command."\n";
system($command);

$command = "head ".$tmpfile.".sam.h";
print STDERR `$command`;

#print STDERR `ls -l`;
#print STDERR `which samtools`;
# system("tar zxvf ./genomes.tgz");

  

#print STDERR "converting SAM to BAM\n";
$command = "./samtools view -Sb ".$tmpfile.".sam.h > ".$tmpfile.".bam";
 print STDERR "\t".$command."\n";
system($command);

print STDERR "sorting bam files\n";
$command = "./samtools sort ./".$tmpfile.".bam ".$tmpfile.".s";
print STDERR "\t".$command."\n";
system($command);

#print STDERR "converting BAM to pileup\n";
#$command = "samtools mpileup -f ".$ref."  ".$tmpfile.".s.bam > ".$tmpfile.".pileup";
#print STDERR "\t".$command."\n";
#system($command);

print STDERR "converting BAM to BCF\n";
#options: -f = reference file
# -g = call genotypes (output bcf)
# -D = output per-sample DP depth (high scoring bases)
$command = "./samtools mpileup -Dgf ".$ref."  ".$tmpfile.".s.bam > ".$tmpfile.".bcf";
print STDERR "\t".$command."\n";
system($command);

my $command = "./bcftools view -g ".$tmpfile.".bcf 2> /dev/null";
print STDERR "converting BCF to VCF\n";
print STDERR "\t".$command."\n";
system($command);  #print to stdout

$command = "head ".$tmpfile.".vcf";
print STDERR `$command`;

print STDERR "complete:\n".`wc $tmpfile.bcf`."\n";

exit 0;
