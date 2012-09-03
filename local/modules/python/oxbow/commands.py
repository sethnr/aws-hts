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
descs['samtag'] = '''aligns seqs in tgz files using samtools'''
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
      "-reducer",     "perl splitAlignTags.pl -jobs 100 -header ./sam.header",
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
      "-mapper",     	"perl splitAlignTags.pl -jobs 100 -header ./sam.header",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/stripUnaligned.pl#stripUnaligned.pl",
      "-cacheFile",   	"s3n://hts-analyses/scripts/perl/splitAlignTags.pl#splitAlignTags.pl",
      "-cacheFile",   "$SAM_HEADER#sam.header"
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

descs['splitbams']='''strip nulls or unaligned from SAM; takes output of 'sam' as input'''
json['splitbams']='''{
  "Name": "split bams aligned", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-D", "mapred.text.key.partitioner.options=-k0,1",
      "-D", "stream.num.map.output.key.fields=2",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",      "perl splitAlignTags.pl -jobs 50 -header sam.header",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/splitAlignTags.pl#splitAlignTags.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''



descs['makebams']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['makebams']='''{
  "Name": "field split & remake bam", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=20",
      "-D", "map.output.key.value.fields.spec=0-1:2-",
      "-D", "mapred.text.key.partitioner.options=-k1.1n,2.1n",
      "-D", "stream.num.map.output.value.fields=3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-reducer",     "perl stripSam.pl --header sam.header",
      "-partitioner", "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/samtools#samtools",
      "-cacheFile",   "s3n://hts-analyses/software/samtools-0.1.18/bcftools#bcftools",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['tag_vcf']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tag_vcf']='''{
  "Name": "field split & make vcf", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=20",
      "-D", "map.output.key.value.fields.spec=0-1:2-",
      "-D", "mapred.text.key.partitioner.options=-k0.1n,1.1n",
      "-D", "stream.num.map.output.value.fields=3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-reducer",     "perl sam_to_vcf.pl -ref ./genomes/AgamP3.I.fa -header ./sam.header -cliptag",
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

descs['tag_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tag_sam']='''{
  "Name": "field split & remake bam", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=20",
      "-D", "map.output.key.value.fields.spec=0,1:2-",
      "-D", "mapred.text.key.partitioner.options=-k1.1n,0.1n",
      "-D", "stream.num.map.output.value.fields=3",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "org.apache.hadoop.mapred.lib.FieldSelectionMapReduce",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
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
descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_12']='''{
  "Name": "field split & remake bam im", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:2-",
      "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k1.1n,2.1n",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

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

descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tag_samcat_01']='''{
  "Name": "field split & remake bam 01", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-D", "map.output.key.value.fields.spec=0,1:2-",
      "-D", "mapred.text.key.partitioner.options=-k0,1",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "perl stripSam.pl --header sam.header --aligned",
      "-reducer",  "cat",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tag_samcat_01']='''{
  "Name": "field split & remake bam 01", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "mapred.reduce.tasks=100",
      "-D", "map.output.key.value.fields.spec=0,1:2-",
      "-D", "mapred.text.key.partitioner.options=-k0,1",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "perl stripSam.pl --header sam.header --aligned",
      "-reducer",  "cat",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_12']='''{
  "Name": "field split & remake bam 12", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:2-",
      "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k1.1n,2.1n",
      "-D", "mapred.text.key.partitioner.options=-k1,2",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_12']='''{
  "Name": "field split & remake bam 12", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,2:2-",
       "-D", "mapred.text.key.partitioner.options=-k1,2",
     "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k1.1n,2.1n",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_02']='''{
  "Name": "field split & remake bam 02", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=0,2:2-",
      "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k0.1n,2.1n",
      "-D", "mapred.text.key.partitioner.options=-k0,2",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_21']='''{
  "Name": "field split & remake bam 21", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=2,1:2-",
      "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k2.1n,1.1n",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''


descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_10']='''{
  "Name": "field split & remake bam 10", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=1,0:2-",
      "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k1.1n,0.1n",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''

descs['tagim_sam']='''takes output of tag, splits based on chr & keys, remakes bam file'''
json['tagim_sam_20']='''{
  "Name": "field split & remake bam 20", 
  "ActionOnFailure": "$ACTION_ON_FAILURE", 
  "HadoopJarStep": { 
    "Jar": "$HADOOP_JAR", 
    "Args": [ 
      "-D", "map.output.key.value.fields.spec=2,0:2-",
      "-D", "mapred.output.key.comparator.class=org.apache.hadoop.mapred.lib.KeyFieldBasedComparator",
      "-D", "mapred.text.key.comparator.options=-k2.1n,0.1n",
      "-D", "stream.num.map.output.value.fields=6",
      "-input",	      "$INPUT",
      "-output",      "$OUTPUT",
      "-mapper",     "cat",
      "-reducer",     "perl stripSam.pl --header sam.header --aligned",
      "-partitioner",     "org.apache.hadoop.mapred.lib.KeyFieldBasedPartitioner",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/stripSam.pl#stripSam.pl",
      "-cacheFile",   "s3n://hts-analyses/scripts/perl/sam_to_vcf.pl#sam_to_vcf.pl",

      "-cacheArchive","$GENOME_REF_JAR#genomes",
      "-cacheFile",   "$SAM_HEADER#sam.header"
    ] 
  }
}'''
