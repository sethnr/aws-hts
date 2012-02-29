# checkout ensembl
echo "CVS password: CVSUSER"
cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl  login
cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl co -r branch-ensembl-62 ensembl-api 
cvs -d :pserver:cvsuser@cvs.sanger.ac.uk:/cvsroot/ensembl co -r branch-ensembl-62 ensembl-tools 

#checkout bioperl
echo "CVS password: cvs"
cvs -d :pserver:cvs@cvs.open-bio.org:/home/repository/bioperl login
cvs -d :pserver:cvs@cvs.open-bio.org:/home/repository/bioperl co bioperl-live