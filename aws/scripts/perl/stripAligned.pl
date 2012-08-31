while(<>) {
  @F = split; 
  print $_ if $F[2] eq '*';
}
