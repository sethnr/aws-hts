#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "bam_to_vcf_tmp";

my $ref = "AgamP3.fa";
my $noclean;
my $vars_only;

GetOptions(
           'fasta|ref=s' => \$ref,
	   'tmp=s' => \$tmpfile,
	   'vars-only|vars_only' => \$vars_only,
	   'noclean' => \$noclean
          );

sub _run_command {
  my $command = shift;
  print STDERR "\t".$command."\n";
  system($command);
  if ( $? == -1 ) { die "  command failed: $!\n"; }
}

my @archives =  `ls *.tgz`;
foreach my $archive (@archives) {
  `tar zxvf $archive`;
  }
_run_command("chmod a+x ./s3cmd/s3cmd ./bgzip ./tabix ");    
$ENV{'PATH'} = '.:./bin/:'.$ENV{'PATH'};
print STDERR $ENV{'PATH'};
print STDERR `ls -l`;

while(<>) {
  print STDERR $_;
  my @files = split;
  my @vcfs = ();
  my @samples = ();
  foreach my $file (@files) {
    my ($tmpfile, $getfile) = split('::',$file);
    print STDERR $getfile." -> ".$tmpfile."\n";

#    print STDERR `ls -l tabix`;
    _run_command("./s3cmd/s3cmd --config ./s3cmd/s3cmd_config get ".$getfile." ".$tmpfile." > /dev/null");    
#    my $s3command = `./s3cmd/s3cmd --config ./s3cmd/s3cmd_config get $getfile $tmpfile `;    
    _run_command("./samtools view -Sb ".$tmpfile." > ".$tmpfile.".bam");
    _run_command("./samtools sort ./".$tmpfile.".bam ".$tmpfile.".s");
    _run_command("./samtools mpileup -Dgf ".$ref."  ".$tmpfile.".s.bam > ".$tmpfile.".bcf");
    _run_command("./bcftools view -g ".$tmpfile.".bcf > ".$tmpfile.".vcf");

    _run_command("perl vcf_prep_for_merge.pl -r GT  ".$tmpfile.".vcf");
    _run_command("./bgzip ".$tmpfile.".vcf");
    _run_command("./tabix -p vcf ".$tmpfile.'.vcf.gz');
    $tmpfile;
    push(@samples, $tmpfile);
  }

  my $cleanup = "combined.vars "; my $vcfs = "";
  foreach my $sample (@samples) {
    $vcfs .= $sample.".vcf.gz ";
    $cleanup .= $sample."* ";
#    push(@vcfs, $sample.".vcf.gz");
  }

#  _run_command("perl vcfs_merge.pl ".$vcfs);
  _run_command("./vcf-merge -s ".$vcfs) unless $vars_only;
  if ($vars_only) {
    _run_command("./vcf-merge -s ".$vcfs." > combined.vars");
    open(COMBINED,"<combined.vars");
    while(<COMBINED>) {
      my @F = split;
      next if $F[0] eq 'File';
      next if $F[4] eq '.';
      print STDOUT $_;
    }
  }

  print STDERR "cleaning ".`ls -l $cleanup`;
  _run_command("rm ".$cleanup) unless $noclean;
}

exit 0;
