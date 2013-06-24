#!/usr/bin/perl
use strict;
use Getopt::Long;

my $tmpfile = "sam_tmp";

my $ref = "reference.vcf.gz";
my $step = "10k";
my $chunk = "0.5m";
my $ld_script = "LDx.pl";
my $header;

GetOptions(
           'ref=s' => \$ref,
	   'tmp=s' => \$tmpfile,
	   'step=s' => \$step,
	   'chunk=s' => \$chunk,
	   'ld=s', => \$ld_script,
	   'replace_header|header=s' => \$header
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

print STDERR `ls`;
system('chmod 755 samtools');
$ENV{'PATH'} = '.:./bin/:'.$ENV{'PATH'};
print STDERR $ENV{'PATH'};
print STDERR `ls -l`;

print STDERR "dumping infile\n";
open (INFILE , ">$tmpfile");
if($header) {
  open (HEADER , "$header");
  while(<HEADER>) {print INFILE $_;}
  close(HEADER);
}
while(<>) {
  next if($header && m/^\@/gi);
  @F = split;
  print INFILE $_;
  }
close(INFILE);

print STDERR `head -n 10 $tmpfile`;

print STDERR "sorting infile\n";
_run_command("samtools view -Sb ".$tmpfile." > ".$tmpfile.".bam");    
_run_command("samtools sort ".$tmpfile.".bam ".$tmpfile.".s");    
_run_command("samtools view ".$tmpfile.".s.bam > ".$tmpfile);
print STDERR "sorted infile\n";

my $fi = 1;

print STDERR "starting new outfile $fi\n";
open (OUTFILE , ">$tmpfile.$fi");
open (INFILE, "$tmpfile");

while(<INFILE>) {
  @F = split;
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
      print STDERR "LAST FILE = ".$lastChr.".".$firstPos.".".$lastPos."\n";
      _run_command("perl ".$ld_script." ".$tmpfile.".".$fi." ".$ref);
      $fi++;
      print STDERR "starting new file $fi\n";
      
      open(OUTFILE,">$tmpfile.$fi");
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
print STDERR "LAST FILE = ".$lastChr.".".$firstPos.".".$lastPos."\n";
_run_command("perl ".$ld_script." ".$tmpfile.".".$fi." ".$ref);
print STDERR "finished!\n";

exit 0;
