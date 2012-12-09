#!/usr/bin/python

from string import Template

#from oxbow.config import config
#print "HADOOP JAR = " + config.HADOOP_JAR

"""
json text strings for all commands, and required values to be passed (for error reporting)
"""

json = {}
descs = {}
default_args = {}

# PREPROCESS
descs['preprocess'] = '''take manifest of sam output files, split into smaller chunks and tgz for processing'''
json['preprocess'] = '''{
  "Name": "Preprocess short reads",
  "ActionOnFailure": "$ACTION_ON_FAILURE",
  "HadoopJarStep": {
    "Jar": "$HADOOP_JAR",
    "Args": [
      "-D", "mapred.reduce.tasks=1",
      "-input",       "$INPUT",
      "-output",      "${S3DIR}null",
      "-mapper",      "$CROSSBOW_EMR/Copy.pl  --compress=gzip --stop=0 --maxperfile=500000 --s --push=$OUTPUT",
      "-inputformat", "org.apache.hadoop.mapred.lib.NLineInputFormat",
      "-cacheFile",   "$CROSSBOW_EMR/Get.pm#Get.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Counters.pm#Counters.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Util.pm#Util.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Tools.pm#Tools.pm",
      "-cacheFile",   "$CROSSBOW_EMR/AWS.pm#AWS.pm",
      "-cacheFile",   "$CROSSBOW_EMR/fastq-dump64#fastq-dump",
      "-cacheFile",   "$CROSSBOW_EMR/samtools64#samtools"
      $OPTS
    ]
  }
}'''

# ALIGN WITH Bowtie (sam)
descs['align'] = '''align seqs in tgz files using align.pl script from crossbow'''
json['align'] = '''{
  "Name": "Align with Bowtie", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-input",       "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "$CROSSBOW_EMR/Align.pl  --discard-reads=0 --ref=$GENOME_REF_JAR --destdir=/mnt/16326 --partlen=1000000 --qual=phred33 --truncate=0  -- --partition 1000000 -t --hadoopout --startverbose -M 1",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "$CROSSBOW_EMR/Get.pm#Get.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Counters.pm#Counters.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Util.pm#Util.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Tools.pm#Tools.pm",
      "-cacheFile",   "$CROSSBOW_EMR/AWS.pm#AWS.pm"
      $OPTS
    ]
  }
}'''


# ALIGN WITH Bowtie (sam)
descs['sam'] = '''align seqs in tgz files using align.pl script from crossbow'''
json['sam'] = '''{
  "Name": "Align with Bowtie (sam)", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-input",       "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "$CROSSBOW_EMR/Align.pl  --discard-reads=0 --ref=$GENOME_REF_JAR --destdir=/mnt/16326 --partlen=1000000 --qual=phred33 --truncate=0  -- --partition 1000000 -t --hadoopout --startverbose -M 1 --sam",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "$CROSSBOW_EMR/Get.pm#Get.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Counters.pm#Counters.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Util.pm#Util.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Tools.pm#Tools.pm",
      "-cacheFile",   "$CROSSBOW_EMR/AWS.pm#AWS.pm"
      $OPTS
    ]
  }
}'''



descs['samsplit']='''takes output of sam splits based on chr & keys, remakes bam file'''
json['samsplit']='''{
  "Name": "field split & remake sam $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=$NO_SAM_TASKS",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl sam_split_by_chunk.pl --chunk $CHUNK_SIZE -flank $FLANK_SIZE ",
      "-reducer",     "perl sam_strip.pl --header sam.header --aligned --cliptag",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_strip.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_strip.pl#sam_strip.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_split_by_chunk.pl#sam_split_by_chunk.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['bamsplit']='''takes bam, splits based on chr & keys, remakes sam file'''
json['bamsplit']='''{
  "Name": "field split & remake sam $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.map.tasks=1",
      "-D", "mapred.reduce.tasks=$NO_SAM_TASKS",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl bam_split_by_chunk.pl --chunk $CHUNK_SIZE -flank $FLANK_SIZE ",
      "-reducer",     "perl sam_strip.pl --header sam.header --aligned --cliptag",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_strip.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_strip.pl#sam_strip.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/bam_split_by_chunk.pl#bam_split_by_chunk.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''
descs['gbamsplit']='''takes bam, splits based on chr & keys, remakes sam file'''
json['gbamsplit']='''{
  "Name": "field split & remake sam $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=$NO_SAM_TASKS",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl bam_split_by_chunk.pl --chunk $CHUNK_SIZE -flank $FLANK_SIZE -gz",
      "-reducer",     "perl sam_strip.pl --header sam.header --aligned --cliptag",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_strip.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_strip.pl#sam_strip.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/bam_split_by_chunk.pl#bam_split_by_chunk.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''



descs['samtobam']='''takes output of sam / splitsam, makes bam file (no reduce)'''
json['samtobam']='''{
  "Name": "sam to (sorted) bam - samtools", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl sam_to_bam.pl",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheFile", "$GENOME_REF_I#ref_genome.fa",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['samdepth']='''takes output of sam / splitsam, calculates depth (samtools), sum by posn (mapred:aggregate) '''
json['samdepth']='''{
  "Name": "read depth", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl sam_depth.pl --aggregate",
      "-reducer",     "aggregate",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_depth.pl#sam_depth.pl",

      "-cacheFile", "$GENOME_REF_I#ref_genome.fa",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['samdepthgrep']='''takes output of sam / splitsam, calculates depth (samtools), sum by posn (mapred:aggregate) '''
json['samdepthgrep']='''{
  "Name": "read depth", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl sam_depth.pl --aggregate",
      "-reducer",     "grep -e 0.10076461 -e 0.10078610 -e 6.222023",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_depth.pl#sam_depth.pl",

      "-cacheFile", "$GENOME_REF_I#ref_genome.fa",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['depthseg']='''segment depth values via bioconductor / fastseg '''
json['depthseg']='''{
  "Name": "fastseg depth", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.map.tasks=400",
      "-D", "mapred.reduce.tasks=1",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "R --slave --vanilla -f fastseg_depth.R",
      "-reducer",     "R --slave --vanilla -f fastseg_depth.R --args minseg=1000",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/scripts/R/fastseg_depth.R#fastseg_depth.R",

      "-cacheFile", "$GENOME_REF_I#ref_genome.fa",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['depthwindow']='''chunk values by bases and run sliding window '''
json['depthwindow']='''{
  "Name": "sliding window mean (depth) $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl splitByChunk.pl --chunk 1m -flank 10k",
      "-reducer",     "R --slave --vanilla -f slidingWin_depth.R --args window=10000 step=1000 fillBases=T " ,
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/scripts/R/slidingWin_depth.R#slidingWin_depth.R",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/cpv_split_by_chunk.pl#splitByChunk.pl",

    ] 
  }
}'''

descs['meanbyfeat']='''chunk values by bases and run sliding window '''
json['meanbyfeat']='''{
  "Name": "feature mean $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl splitByFeature.pl --feats gene_groups.txt",
      "-reducer",     "R --slave --vanilla -f meanByBlock_depth.R" ,
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/scripts/R/meanByBlock_depth.R#meanByBlock_depth.R",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/cpv_split_by_feature2.pl#splitByFeature.pl",
      "-cacheFile",   "s3n://hts-analyses/resources/ortho_gene_groups.txt#gene_groups.txt",

    ] 
  }
}'''


descs['uniqwin']='''chunk values by bases and run sliding window '''
json['uniqwin']='''{
  "Name": "remove tag from chr name, get unique windows $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=6",
      "-input",	      "$INPUT",
      "-output",      "$INPUTuniq",
      "-mapper",      "perl clipChunkTag.pl'",
      "-reducer",     "uniq",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/misc/clip_chunk_tag.pl#clipChunkTag.pl",
    ] 
  }
}'''



descs['vcf']='''use samtools (mpileup) to create vcf file; takes output of 'sam' as input'''
json['vcf']='''{
  "Name": "sam to vcf (samtools) $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl sam_to_vcf.pl -ref ./ref_genome.fa",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheFile", "$GENOME_REF_I#ref_genome.fa",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['strip']='''strip nulls or unaligned from SAM; takes output of 'sam' as input'''
json['strip']='''{
  "Name": "stripsam", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['unaligned']='''split bams taking only unaligned (random partitioning)'''
json['unaligned']='''{
  "Name": "parse out only unaligned seqs $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=1",
      "-D", "stream.num.map.output.key.fields=2",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --unaligned",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['checkkeys']='''for each file, show unique keys'''
json['checkkeys']='''{
  "Name": "check keys", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "perl headCount.pl",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/headCount.pl#headCount.pl",
    ] 
  }
}'''


descs['make_manifest']='''make manifest of files for each sample folder'''
json['make_manifest']='''{
  "Name": "make_manifest", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=3",
      "-input",	      "$S3DIR/control/sample_folder_list.txt",
      "-output",      "$S3DIR/control/sample_manifest.txt",
      "-mapper",      "perl s3_make_sample_manifest.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/s3_make_sample_manifest.pl#s3_make_sample_manifest.pl"
    ] 
  }
}'''

descs['combined_vcf']='''make manifest of files for each sample folder'''
json['combined_vcf']='''{
  "Name": "combined_vcf $S3DIR", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.map.tasks=100",
      "-D", "mapred.reduce.tasks=0",
      "-input",      "$INPUT",
      "-output",     "$S3DIR/null",
      "-mapper",      "perl sams_to_vcf.pl -ref ./ref_genome.fa -vars-only -outdir /hts-outputs/Fd09/combined_vcf -chunk $CHUNK_SIZE -no-dups",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/vcf_prep_for_merge.pl#vcf_prep_for_merge.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/vcfs_merge.pl#vcfs_merge.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sams_to_vcf.pl#sams_to_vcf.pl",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/vcf-merge#vcf-merge",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile","s3n://hts-analyses/software/s3cmd-1.1.0.tgz#s3cmd.tgz",
      "-cacheFile","s3n://hts-analyses/software/tabix/tabix#tabix",
      "-cacheFile","s3n://hts-analyses/software/tabix/bgzip#bgzip",
      "-cacheFile","s3n://hts-analyses/lib/vcfPerlMods.tgz#vcfPerlMods.tgz",
      "-cacheFile", "$GENOME_REF#ref_genome.fa"
    ] 
  }
}'''

descs['theta']='''make manifest of files for each sample folder'''
json['theta']='''{
"Name": "vcf theta $S3DIR ", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-D", "mapred.tasktracker.map.tasks.maximum=$MAX_TASKS",
      "-input",      "$INPUT",
      "-output",     "$OUTPUT",
      "-mapper",      "R --slave --vanilla -f vcf_slidingWin_theta.R --args window=10000 step=1000 fillBases=T pool_size=20" ,
      "-cacheFile",   "s3n://hts-analyses/scripts/R/vcf_slidingWin_theta.R#vcf_slidingWin_theta.R",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/vcf-merge#vcf-merge",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile","s3n://hts-analyses/software/s3cmd-1.1.0.tgz#s3cmd.tgz",
      "-cacheFile","s3n://hts-analyses/software/tabix/tabix#tabix",
      "-cacheFile","s3n://hts-analyses/software/tabix/bgzip#bgzip",
      "-cacheFile","s3n://hts-analyses/lib/vcfPerlMods.tgz#vcfPerlMods.tgz",
      "-cacheFile", "$GENOME_REF_I#ref_genome.fa"
    ] 
  }
}'''

