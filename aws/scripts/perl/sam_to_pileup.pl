my $tmpfile = "/tmp/sam_to_pileup_output";

open(TMP,">".$tmpfile.".sam");

print STDERR "dumping input into file\n";

while(<>){
  print TMP "$_";
}
close(TMP);

print STDERR "converting SAM to BAM\n";
my $command = "samtools view -Sb ".$tmpfile.".sam > ".$tmpfile.".bam";
system($command);

$command = "samtools mpileup ".$tmpfile.".bam";
print STDERR "converting BAM to pileup\n$command\n";
system($command);

exit 0;
