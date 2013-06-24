#!/usr/bin/perl

# make manifest from bucket listings

my %files;
my %sizes;
my @samples; 

foreach $arg (@ARGV) {
  my ($sample, $dir) = split("::",$arg);
  print STDERR "getting $sample samples from $dir \n";

  my @responses = `s3cmd ls $dir`;
  print STDERR $#responses. " responses returned\n";

  foreach my $prev_sample (@samples) {
    die "different nos of files found:\n\t".
      $prev_sample." ".$#{$files{$prev_sample}}."\n\t".
      $sample." ".$#responses  if $#responses != $#{$files{$prev_sample}};
  }

  foreach my $response (@responses) {

    my (undef, undef, $size, $file) = split('\s+', $response);
    push (@{$files{$sample}},$file);
    push (@{$sizes{$sample}},$size);    
  }
  push (@samples, $sample);
}


for ($i=0; $i<= $#{$files{$samples[0]}}; $i++) {
  my @line;
  $sizes = 0;
  for ($j=0; $j<= $#samples; $j++) {
    $sample = $samples[$j];
    push @line, $sample."::".$files{$sample}[$i];
    $sizes += $sizes{$sample}[$i];
  }
  print join("\t",@line)."\n" unless $sizes == 0;
}
