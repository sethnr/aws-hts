#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "sam_tmp";

my $ref = "reference.vcf.gz";
my $step;
my $chunk;
my $ld_script = "LDx.pl";


GetOptions(
           'ref=s' => \$ref,
	   'tmp=s' => \$tmpfile,
	   'step=s' => \$step,
	   'chunk=s' => \$chunk,
	   'ld=s', => \$ld_script
          );

sub _run_command {
  my $command = shift;
  print STDERR "\t".$command."\n";
  system($command);
  if ( $? == -1 ) { die "  command failed: $!\n"; }
}

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

my ($firstPos, $lastPos, $lastChr);
$lastPos = -1;
my @header;
my @F;
print STDERR "starting new outfile\n";
open (OUTFILE , ">$tmpfile");

while(<>) {
  @F = split;
  next if $F[0] eq 'File';
  if ($F[0] =~ m/^\@/gi) {
    chomp;
    push(@header,$_);
    print OUTFILE $_."\n";
    next;
  }

  unless($firstPos || ($F[0] =~ m/^#.*/)) {$firstPos = $F[3];}
  if ($lastPos > -1) {
    if ((($lastPos + $chunk) < $F[3]) || ($lastChr ne $F[2])) {
      print STDERR "found break larger than ".$chunk." \n";
      print STDERR "lastPos = ".$lastPos." chunk = ".$chunk." F3 = ".$F[3]." lastChr = ".$lastChr." F2 = ".$F[2]."\n";
      close(OUTFILE);
      print STDERR "LAST FILE = ".$out_file
      my $out_file = $lastChr.".".$firstPos.".".$lastPos.".sam";
      _run_command("samtools view -Sb ".$tmpfile." > ".$tmpfile.".bam");    
      _run_command("samtools sort ".$tmpfile.".bam ".$tmpfile.".s");    
      _run_command("samtools view ".$tmpfile.".s.bam > ".$out_file);
      _run_command("perl ".$ld_script." ".$out_file." ".$ref);
      
      print STDERR "starting new file\n";
      
      open(OUTFILE,">$tmpfile");
      print OUTFILE join("\n",@header)."\n";
      $firstPos = $F[3];
      $lastChr = $F[2];
    }
  }
  $lastPos = $F[3];
  $lastChr = $F[2];
  
  print OUTFILE $_;
}

print STDERR "end of infile\n";
print STDERR "lastPos = ".$lastPos." chunk = ".$chunk." F3 = ".$F[3]." lastChr = ".$lastChr." F2 = ".$F[2]."\n";
close(OUTFILE);
my $out_file = $lastChr.".".$firstPos.".".$lastPos.".sam";
_run_command("samtools view -Sb ".$tmpfile." > ".$tmpfile.".bam");    
_run_command("samtools sort ".$tmpfile.".bam ".$tmpfile.".s");    
_run_command("samtools view ".$tmpfile.".s.bam > ".$out_file);
_run_command("perl ".$ld_script." ".$out_file." ".$ref);

exit 0;
