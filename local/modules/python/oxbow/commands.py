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

#split paired-end reads
descs['splitpe'] = '''split paired-end reads into separate files'''
json['splitpe'] = '''{
  "Name": "split paired-end reads",
  "ActionOnFailure": "$ACTION_ON_FAILURE",
  "HadoopJarStep": {
    "Jar": "$HADOOP_JAR",
    "Args": [
      "-D", "mapred.reduce.tasks=1",
      "-input",       "$INPUT",
      "-output",      "${S3DIR}null",
      "-mapper",      "./splitPairedEndReads.pl $OUTPUT",
      "-inputformat", "org.apache.hadoop.mapred.lib.NLineInputFormat",
   
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/splitPairedEndReads.pl#splitPairedEndReads.pl",
      $OPTS
    ]
  }
}'''

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


# ALIGN WITH Bowtie (sam)
#descs['samtag'] = '''aligns seqs in tgz files using samtools'''
json['samtag'] = '''{
  "Name": "Align with Bowtie (sam) and tag", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=20",
      "-D", "mapred.text.key.partitioner.options=-k3,4",
      "-D", "stream.num.map.output.key.fields=4",
      "-input",       "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "$CROSSBOW_EMR/Align.pl  --discard-reads=0 --ref=$GENOME_REF_JAR --destdir=/mnt/16326 --partlen=1000000 --qual=phred33 --truncate=0  -- --partition 1000000 -t --hadoopout --startverbose -M 1 --sam",
      "-reducer",     "perl splitAlignTags.pl -jobs 50 -header ./sam.header",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "$CROSSBOW_EMR/Get.pm#Get.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Counters.pm#Counters.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Util.pm#Util.pm",
      "-cacheFile",   "$CROSSBOW_EMR/Tools.pm#Tools.pm",
      "-cacheFile",   "$CROSSBOW_EMR/AWS.pm#AWS.pm",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/splitAlignTags.pl#splitAlignTags.pl",
    "-cacheFile",     "$SAM_HEADER#sam.header"
     $OPTS
    ]
  }
}'''


#default_args['align'] = {'input':'s3://ngousso/full/preproc',
#                     'output':'s3://ngousso/full/sam'}

descs['tag']='''use custom perl script to strip-non-aligned seqs out of file,
tag with index for sorting by mapreduce'''
json['tag']='''{
  "Name": "strip non-matched lines & group-tag", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-input",       "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     	"perl splitAlignTags.pl -jobs $NO_SAM_FILES -header ./sam.header",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/stripUnaligned.pl#stripUnaligned.pl",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/splitAlignTags.pl#splitAlignTags.pl",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['splitsam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['splitsam']='''{
  "Name": "field split & remake sam", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=$NO_SAM_TASKS",
      "-D", "map.output.key.value.fields.spec=0:2-",
      "-D", "mapred.text.key.partitioner.options=-k0",
      "-D", "stream.num.map.output.value.fields=2",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl splitAlignTags.pl -jobs $NO_SAM_FILES -header sam.header",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned --cliptag",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/splitAlignTags.pl#splitAlignTags.pl",

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
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/splitByChunk.pl#splitByChunk.pl",

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
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/splitByFeature2.pl#splitByFeature.pl",
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
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/clipChunkTag.pl#clipChunkTag.pl",
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
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
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


