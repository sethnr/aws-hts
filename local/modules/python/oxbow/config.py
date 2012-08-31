#!/usr/bin/python

import re
import sys
"""
defaults for all commonly used EMR bits
"""
ACTION_ON_FAILURE = 'TERMINATE_JOB_FLOW'
HADOOP_JAR = '/home/hadoop/contrib/streaming/hadoop-0.20-streaming.jar'
CROSSBOW_EMR = 's3n://crossbow-emr/1.1.2'
GENOME_REF = 'AgamP3.I.fa'
GENOME_REF_JAR = 's3n://hts-analyses/bowtie_jars/AgamP3.jar'
SAM_HEADER = 's3n://hts-analyses/resources/AgamP3.header'
EMR_COMMAND = '/Users/Seth/Work/EC2/emrtools/elastic-mapreduce'

def get_conf_hash():
  mod_dict = sys.modules[__name__].__dict__
  def_dict = {}
  for var in mod_dict:
    if re.match('^[A-Z,_]+$',var):
      #print var+"->"+mod_dict[var]
      def_dict[var] = mod_dict[var]
  return def_dict