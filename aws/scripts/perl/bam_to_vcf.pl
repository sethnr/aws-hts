#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "bam_to_vcf_tmp";

my $ref = "AgamP3.fa";

GetOptions(
           'fasta|ref=s' => \$ref,
	   'tmp=s' => \$tmpfile
          );

sub _run_command {
  my $command = shift;
  print STDERR "\t".$command."\n";
  system($command);
}


<<<<<<< HEAD
while(<>) {
  print STDERR $_;
  my @files = split;
  my @samples = ();
  foreach my $file (@files) {
    my ($tmpfile, $getfile) = split('::',$file);
    print STDERR $getfile." -> ".$tmpfile."\n";
    _run_command("s3cmd get ".$getfile." ".$tmpfile);    
    _run_command("samtools view -Sb ".$tmpfile." > ".$tmpfile.".bam");
    _run_command("samtools sort ./".$tmpfile.".bam ".$tmpfile.".s");
    _run_command("samtools mpileup -Dgf ".$ref."  ".$tmpfile.".s.bam > ".$tmpfile.".bcf");
    _run_command("bcftools view -g ".$tmpfile.".bcf 2> /dev/null > ".$tmpfile.".vcf");
    $tmpfile .= ".vcf";
    push(@samples, $tmpfile);
  }

  _run_command("perl vcfs_merge.pl ".join(" ",@samples));
}
  exit 0;
=======
print STDERR "sorting bam files\n";
$command = "samtools sort ./".$tmpfile.".bam ".$tmpfile.".s";
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
$command = "samtools mpileup -Dgf ".$ref."  ".$tmpfile.".s.bam > ".$tmpfile.".bcf";
print STDERR "\t".$command."\n";
system($command);

# my $command = "bcftools view -g ".$tmpfile.".bcf 2> /dev/null";
my $command = "bcftools view -gv ".$tmpfile.".bcf 2> /dev/null";
print STDERR "converting BCF to VCF\n";
print STDERR "\t".$command."\n";
system($command);
print STDERR "complete:\n".`wc $tmpfile.bcf`."\n";

exit 0;
>>>>>>> 3e00d747d99b01163e873ce4aaa3c7e94d1a18d7
