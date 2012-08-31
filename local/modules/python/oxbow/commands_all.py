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
#default_args['preprocess'] = {'input':'ngousso.manifest',
#                              'OUT':'preproc',
#                              'NULL':'null'}


# ALIGN WITH Bowtie (sam)
descs['align'] = '''align seqs in tgz files using align.pl script from crossbow'''
json['align'] = '''{
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
descs['sam'] = '''aligns seqs in tgz files using samtools'''
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
      "-mapper",     	"perl splitAlignTags.pl -jobs 100 -header ./AgamP3.header",
      "-cacheFile",   	"s3n://hts-analyses/resources/AgamP3.header#AgamP3.header",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/stripUnaligned.pl#stripUnaligned.pl",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/splitAlignTags.pl#splitAlignTags.pl"
    ] 
  }
}'''

descs['vcf']='''use samtools (mpileup) to create vcf file; takes output of 'sam' as input'''
json['vcf']='''{
  "Name": "sam to vcf (samtools)", 
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
      "-reducer",     "perl sam_to_vcf.pl -ref ./genomes/AgamP3.I.fa -header ./sam.header",
      "-cacheFile",   "$CROSSBOW_EMR/bowtie64#bowtie",
      "-cacheFile",   "$CROSSBOW_EMR/samtools64#samtools64",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
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

descs['fs_strip']='''strip nulls or unaligned from SAM; takes output of 'sam' as input'''
json['fs_strip']='''{
  "Name": "stripsam fs.reduce 1-2:3-", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-D", "stream.num.map.output.key.fields=3",
      "-D", "reduce.output.key.value.fields.spec=1-2:3-",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl stripSam.pl",
      "-reducer",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['strip_nr']='''strip nulls or unaligned from SAM; takes output of 'sam' as input'''
json['strip_nr']='''{
  "Name": "strip part noFS nored", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=0",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.key.fields=3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl stripSam.pl",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['strip_cat']='''strip nulls or unaligned from SAM; takes output of 'sam' as input'''
json['strip_cat_k12']='''{
  "Name": "strip part cat k1-2", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=50",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.value.fields=3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl stripSam.pl",
      "-reducer", "cat",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_12_1']='''{
  "Name": "fs.map 1,2:1", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:1",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_12_1']='''{
  "Name": "fs.map 1,2:1", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:1",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_12_2']='''{
  "Name": "fs.map 1,2:2", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:2",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_12_3']='''{
  "Name": "fs.map 1,2:3", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_01_2-']='''{
  "Name": "fs.map 0,1:2-", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=0,1:2-",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_01_2-_cat']='''{
  "Name": "fs.map 0,1:2-", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=0,1:2-",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-reducer",     "cat",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''



descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_01_0-']='''{
  "Name": "fs.map 0,1:0-", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=0,1:0-",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['fs']='''strip key vals, return rest of BAM. takes output of 'tag' as input'''
json['fs_01_4-']='''{
  "Name": "fs.map 0,1:4-", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=0,1:4-",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''
