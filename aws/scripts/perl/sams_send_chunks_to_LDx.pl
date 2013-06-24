#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "LDx.input";

my $ref = "AgamP3.fa";
my $noclean;
my $vars_only;
my $noDups;
my $out_dir;
my $step;
my $chunk;

GetOptions(
           'fasta|ref=s' => \$ref,
	   'tmp=s' => \$tmpfile,
	   'out|outdir=s' => \$out_dir,
	   'vars-only|vars_only' => \$vars_only,
	   'noclean' => \$noclean,
	   'step=s' => \$step,
	   'chunk=s' => \$chunk,
	   'no-dups|nd' => \$noDups
          );

if($noDups) {$noDups = " -d ";} else{$noDups = "";}

sub cleanNo {
  my $val = shift;
  print STDERR "cleaning ".$val;
  $val =~ s/,//g;
  if ($val =~ m/^([\d,\.]+)['k','kb']/gi) {
    $val = $1 * 1000;
  }
  elsif($val =~ m/^([\d,\.]+)['m','mb']/gi) {
    $val = $1 * 1000000;
  }
  print STDERR " -> ".$val."\n";
  return $val;
}

$step = cleanNo($step);
$chunk = cleanNo($chunk);



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

    _run_command("perl vcf_prep_for_merge.pl  ".$tmpfile.".vcf");
    _run_command("./bgzip ".$tmpfile.".vcf");
    _run_command("./tabix -p vcf ".$tmpfile.'.vcf.gz');
    $tmpfile;
    push(@samples, $tmpfile);
  }

  my $cleanup = "combined.vcf combined.tmp.vcf "; my $vcfs = "";
  foreach my $sample (@samples) {
    $vcfs .= $sample.".vcf.gz ";
    $cleanup .= $sample."* ";
#    push(@vcfs, $sample.".vcf.gz");
  }

#  _run_command("perl vcfs_merge.pl ".$vcfs);
  _run_command("./vcf-merge ".$noDups." -s ".$vcfs." > combined.tmp.vcf");
my @header;
my ($firstPos, $lastPos, $lastChr);
$lastPos = -1;
    open(OUTFILE,">combined.vcf");
    open(COMBINED,"<combined.tmp.vcf");

    while(<COMBINED>) {
      my @F = split;
      next if $F[0] eq 'File';
      if ($F[0] =~ m/^\#/gi) {
	chomp;
	push(@header,$_);
 	print OUTFILE $_."\n";
	next;
      }

      if ($vars_only) {
	next if $F[4] eq '.';
      }
      unless($firstPos || ($F[0] =~ m/^#.*/)) {$firstPos = $F[1];}
      if ($lastPos > -1) {
	if ((($lastPos + $chunk) < $F[1]) || ($lastChr ne $F[0])) {
	  print STDERR "lastPos = ".$lastPos." chunk = ".$chunk." F1 = ".$F[1]." lastChr = ".$lastChr." F0 = ".$F[0]."\n";
	  my $out_file = $lastChr.".".$firstPos.".".$lastPos.".vcf";
	  _run_command("./s3cmd/s3cmd --config ./s3cmd/s3cmd_config put ./combined.vcf  s3:/".$out_dir."/".$out_file);    
	  close(OUTFILE);
	  print STDERR "found break larger than ".$chunk." splitting into two files\n";
	  open(OUTFILE,">combined.vcf");
	  print OUTFILE join("\n",@header)."\n";
	  $firstPos = $F[1];
	  $lastChr = $F[0];
	}
      }
      $lastPos = $F[1];
      $lastChr = $F[0];
      
      print OUTFILE $_;
    }
  


 my $out_file = $lastChr.".".$firstPos.".".$lastPos.".vcf";
  _run_command("./s3cmd/s3cmd --config ./s3cmd/s3cmd_config put ./combined.vcf  s3:/".$out_dir."/".$out_file);    


  print STDERR "cleaning ".`ls -l $cleanup`;
  _run_command("rm ".$cleanup) unless $noclean;
}

exit 0;
