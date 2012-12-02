#get gpg key
#gpg --keyserver subkeys.pgp.net --recv-key 381BA480
#gpg -a --export 381BA480 > jranke_cran.asc
#sudo apt-key add jranke_cran.asc

# update all rpms
#echo "deb http://cran.ma.imperial.ac.uk/bin/linux/debian lenny-cran/" | sudo tee -a /etc/apt/sources.list
#sudo apt-get update

#fortran not installing properly, remove and replace:
#sudo dpkg --purge r-base-dev
#sudo dpkg --purge gfortran

# update R & dependencies (fucking GCC?!)
#sudo apt-get -t lenny-cran install --yes --force-yes gfortran libgfortran3 pkg-config 
#sudo apt-get -t lenny-cran install --yes --force-yes libxml2-dev libcurl4-gnutls-dev
#sudo apt-get -t lenny-cran install --yes --force-yes r-base r-base-dev

# update GCC! WTF?!
#sudo apt-get -t lenny-cran install --yes --force-yes gcc-4.3
#sudo ln -fs `which gcc-4.3` `which gcc`

# install bioC packages
sudo R -e "options(repos=structure(c(CRAN=\"http://cran.ma.imperial.ac.uk/\"))); install.packages(\"zoo\")"
sudo R -e "source('http://bioconductor.org/biocLite.R');biocLite(character(0), ask=FALSE);biocLite(c(\"fastseg\",\"snpStats\",\"ShortRead\",\"GGtools\"))"