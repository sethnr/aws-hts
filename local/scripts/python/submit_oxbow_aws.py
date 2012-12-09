#!/usr/bin/python

import sys
import os
import re
import argparse
import string
import tempfile
import subprocess  # pass subcommands to system
import shlex   # get token, unix-style lexical analysis
               #import commands (deprecated in favour of subprocess)

from oxbow import config 
import oxbow.commands

_args = {}
_steps = []
_step_args = {}
_emr_args = []
_input = 'input'
_output = 'output'
_log_dir = ''
_s3dir = ''


def parse_steps_string(step_string):
  for step in string.split(step_string,':'):
    _steps.append(step)
  
def parse_steps_file(filename):
  print "sorry, not written this thing yet."

  
def get_steps_json_string(input, s3dir):
  '''starting at input directory / file, write json for each step, output to s3dir/stepName
read in next step from s3dir/prevName'''
  jsons = []
  for step in _steps:
    json = oxbow.commands.json[step]
    #    template = template.safe_substitute(oxbow.commands.default_args[commands[0]])
    json = fill_config_defaults(json)
    # S3:// work around for preprocess step
    if step == 'preprocess':
      capS3dir = string.replace(s3dir, 's3://', 'S3://')
      if step in _step_args:
        json = fill_step_defaults(json, input, capS3dir+step, s3dir, _step_args[step])
      else:
        json = fill_step_defaults(json, input, capS3dir+step, s3dir, None)
    else:
      if step in _step_args:
        json = fill_step_defaults(json, input, s3dir+step, s3dir, _step_args[step])
      else:
        json = fill_step_defaults(json, input, s3dir+step, s3dir, None)

    jsons.append(json)
    input = s3dir+step
  return '[\n'+string.join(jsons,',\n')+'\n]'

def fill_config_defaults(json):
  '''from oxbow.config object, make hash of variables and fill into json objects from oxbow.commands'''
  defs = config.get_conf_hash()
  json_t = string.Template(json)
  json = json_t.safe_substitute(defs)
  return json

def fill_step_defaults(json, input, output, dir, step_args):
  '''fill in json template with values from input string, add extra vals to end of json block
nb. this can over-write previous commands further up the same json block'''
  opts = ''
  if step_args:
    for key in step_args:
      opts += ', "-'+key+'",  "'+step_args[key]+'"'  
  defs = {}
  json_t = string.Template(json)
  # print json
  defs['INPUT'] = input
  defs['OUTPUT'] = output
  defs['S3DIR'] = dir
  defs['OPTS'] = opts
  #  print defs
  
  json = json_t.substitute(defs)
  
  return json


def parse_args(args):
  '''parse all arguments from the input string'''
  parser = argparse.ArgumentParser()
  parser.add_argument('-s', '--steps')
  parser.add_argument('-sf', '--stepfile')
  parser.add_argument('-i', '--input')
  parser.add_argument('-d', '--dir', '--s3dir', dest='s3dir')
  parser.add_argument('--dryrun', action='store_true')
  parser.add_argument('--alive', action='store_true')
  parser.add_argument('--bioc','--BioC', action='store_true')
  parser.add_argument('--rdeps','--Rdeps', action='store_true')
  parser.add_argument('--s3','--s3install', action='store_true')
  parser.add_argument('--stream')
  parser.add_argument('--num_instances', '--instances')
  parser.add_argument('--instance_type','--instance-type','--itype')

  global _args 
  _args, rem_args = parser.parse_known_args()
  # check through remaining for step-specific arguments
  # pass rest of commands to EMR
  _emr_args = parse_step_defaults(rem_args)

  # check steps are defined one way or t'other
  if _args.steps:
    parse_steps_string(_args.steps)
  elif _args.stepfile:
    parse_steps_file(_args.stepfile)
  else:
    raise IOError("needs either steps or stepfile to be defined")

  #  if _args.alive or _args.stream:
  if _args.alive:
    config.ACTION_ON_FAILURE = 'CANCEL_AND_WAIT'

  #check output dir / s3dir is defined
  if not (_args.s3dir) or (not re.match('^s3://',_args.s3dir, re.I)):  #absolute output
    raise IOError("needs output / s3dir to be defined with an absolute path")
  elif not re.match('.*\/$',_args.s3dir):
    _args.s3dir += '/'

  #check input is defined, or guess
  if not (_args.input):
    _args.input = _args.s3dir+'input'
    sys.stderr.write('input file not defined, trying:\n  '+_args.input+'\n')
  elif not re.match('^s3://',_args.input, re.I):
    _args.input = _args.s3dir+_args.input
    sys.stderr.write('input is relative, trying:\n  '+_args.input+'\n')

  if not _args.num_instances: _args.num_instances=1
  if not _args.instance_type: _args.instance_type='c1.xlarge'

def parse_step_defaults(args):
  '''for anything with a "--step.key=value" or "--step.key value" format, put into hash with key-val pairs for each step
  append the correctly formatted values to the end of the json file'''
  rem_args = []
  while len(args):
    arg = args.pop(0)
    print arg+"!\n"
    m = re.match('--(\w+)\.(\w+)=*(\S*)', arg)
    if m:
      step = m.group(1)
      key = m.group(2)
      if m.group(3): val = m.group(3)
      else: val = args.pop(0)
      print step+": "+key+" > "+val
      if step not in _step_args:
        _step_args[step] = {}
      _step_args[step][key] = val
      print _step_args[step][key]
    else:
      rem_args.append(arg)  
  return rem_args

def submit_emr_job(json_file):
    #tsocks ~/Work/EC2/emrtools/elastic-mapreduce
  jobflow_name = string.join(_steps,':')
  emr_args = []
  emr_args.append(config.EMR_COMMAND)
  #emr_args.append("echo")
  
  if not _args.stream:
    emr_args.append('--create')

    emr_args.append('--bootstrap-action s3://hts-analyses/scripts/sh/setup.sh --bootstrap-name "Setup paths & libs" ')
    emr_args.append('--bootstrap-action s3://elasticmapreduce/bootstrap-actions/configurations/latest/memory-intensive --bootstrap-name "Set memory-intensive mode" ')
    emr_args.append('--bootstrap-action s3://elasticmapreduce/bootstrap-actions/configure-hadoop --bootstrap-name "Configure Hadoop" --args "-s,mapred.job.reuse.jvm.num.tasks=1,-s,mapred.tasktracker.reduce.tasks.maximum=8,-s,io.sort.mb=100"')
    emr_args.append('--bootstrap-action s3://elasticmapreduce/bootstrap-actions/add-swap --bootstrap-name "Add Swap" --args "1234"')

    # update R & BioConductor if needed
    if _args.rdeps: emr_args.append('--bootstrap-action s3://hts-analyses/scripts/sh/updateRdeps.sh --bootstrap-name "Update R" --args "1234"')

    if _args.bioc: emr_args.append('--bootstrap-action s3://hts-analyses/scripts/sh/installBioC.sh --bootstrap-name "install BioC" --args "1234"')

    if _args.s3: emr_args.append('--bootstrap-action s3://hts-analyses/scripts/sh/installS3cmd.sh --bootstrap-name "install S3" --args "1234"')

    emr_args.append('--enable-debugging')
    emr_args.append('--log-uri '+_args.s3dir+'log')
    emr_args.append('--name '+jobflow_name)
    emr_args.append('--instance-type '+ _args.instance_type)
    emr_args.append('--num-instances ' + str(_args.num_instances))

    if _args.alive: emr_args.append('--alive')



  emr_args.append('--json '+json_file)
  if _args.stream:
    emr_args.append('--stream '+_args.stream)

  for arg in _emr_args:
    emr_args.append(arg)

  emr_args = shlex.split(' '.join(emr_args))

  print ' '.join(emr_args)
  if not _args.dryrun:
    try:
      outstring = subprocess.check_output(emr_args)
      print outstring
      m = re.match('Created job flow (j-\w+)',outstring)
      if m:
        job_id = m.group(1)
        sys.stderr.write('export AWSLAST='+job_id+'\n')
        subprocess.call('export AWSLAST='+job_id, shell=True)
    except subprocess.CalledProcessError as err:
      sys.stderr.write(err.output)
      sys.stderr.write('command returned with status code '+str(err.returncode)+'\n')




def main():
  parse_args(sys.argv)
  if not _steps: raise IOError("no steps found to perform!")
  json = get_steps_json_string(_args.input,_args.s3dir)
  (tmp_handle, tmpname) = tempfile.mkstemp(prefix='oxbow_',suffix='.json',dir='/tmp/',text=True)
  sys.stderr.write('json written to:\n  '+tmpname+'\n')
  tmp_json = open(tmpname, 'w')
  tmp_json.write(json)
  tmp_json.close()
  submit_emr_job(tmpname)
if __name__ == "__main__":
  main()
