#!/usr/bin/perl

my $tmpfile = "./sam_to_vcf_tmp";

open(TMP,">".$tmpfile.".sam");

print STDERR "dumping input into file\n";
while(<>){
  print TMP "$_";
}
close(TMP);

print STDERR "converting SAM to BAM\n";
$command = "./samtools view -Sb ".$tmpfile.".sam > ".$tmpfile.".bam";
system($command);

print STDERR "converting BAM to BCF\n";
$command = "./samtools mpileup -g ".$tmpfile.".bam > ".$tmpfile.".bcf";
# print STDOUT `./samtools`;
system($command);

my $command = "./bcftools view -g ".$tmpfile.".bcf";
print STDERR "converting BCF to VCF\n$command\n";
system($command);
print STDERR "complete:\n".`wc $tmpfile.*`."\n";

exit 0;
