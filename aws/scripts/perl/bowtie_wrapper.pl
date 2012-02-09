my $tmpfile = "./tmp/run_bowtie";

system("jar -xf ./AgamP3.jar");

open(TMP,">".$tmpfile.".fastq");

print STDERR "dumping input into file\n";

while(<>){
  print TMP "$_";
}
close(TMP);

print STDERR `ls -la`; 
my $command = "bowtie -S -M 1 -t  --12 ".$tmpfile.".fastq ./index/Agam";
print STDERR "running_bowtie\n$command";
print STDOUT `$command`;


exit 0;
