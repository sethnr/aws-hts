#!/usr/bin/perl

# take list of VCF files, copy across from S3 
# move DP4 values into samples col
# NB will only work on single-sample VCFs
# doesn't perform error checking

sub _run_command {
  my $command = shift;
  print STDERR $command."\n";
  system($command);
  if ( $? == -1 ) { die "  command failed: $!\n"; }
}


my @vcfs;

foreach $vcf (@ARGV) {
#  print STDERR $vcf;
  if ($vcf =~ m/^s3:/) {
    @filename = split('/',$vcf);
    my $sample = $filename[-3];    
    my $getfile = $filename[-3]."_".$filename[-1];    
    system("s3cmd get $vcf $getfile");
    $vcf = $getfile;
  }
  _run_command("perl vcf_prep_for_merge.pl $vcf");
  _run_command("./tabix/bgzip $vcf");
  _run_command("./tabix/tabix -p vcf $vcf".'.gz');
  push(@vcfs, $vcf.".gz");
}
_run_command("./vcf-merge -s ".join(" ",@vcfs));
