
my $i = 0;
my %minmax;

while(<>) {
  @F = split; 
  next if $F[0] =~ m/^\@/;
  print $_ if ($i < 100);
  $i++;
  $posn = $F[5];
  $chr = $F[0];
  $key = $F[1];
  $minmax{$chr}{$key}{'min'} = $posn if (!$minmax{$chr}{$key}{'min'});

  $minmax{$chr}{$key}{'min'} = $posn unless($minmax{$chr}{$key}{'min'} < $posn) ; 
  $minmax{$chr}{$key}{'max'} = $posn unless($minmax{$chr}{$key}{'max'} > $posn);
  $minmax{$chr}{$key}{'sum'} += $posn;
  $minmax{$chr}{$key}{'count'}++;

}

print "...\n".$i." lines in file\n";
foreach my $chr(keys %minmax) {
  foreach my $key(keys %{$minmax{$chr}}) {
    print "#".$chr."\t".$key."\t".$minmax{$chr}{$key}{'min'}."\t".$minmax{$chr}{$key}{'max'}.
      "\t".$minmax{$chr}{$key}{'count'}." ".($minmax{$chr}{$key}{'sum'} / $minmax{$chr}{$key}{'count'})."\n";
  }
}
