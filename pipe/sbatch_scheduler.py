#!/usr/bin/env python3
"""
SLURM scheduler script
"""

import os
import sys
import warnings

from snakemake.utils import read_job_properties

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)
job_resources= job_properties["resources"]

cluster_param={}
cluster_param['name'] = job_properties['rule']
cluster_param["threads"] = job_properties.get("threads",1)
cluster_param['mem'] = int(job_resources.get("mem",5))+ 5 #GB + overhead

eprint("Submit job with parameters:\n"+"\n".join(["\t{} : {}".format(key,cluster_param[key]) for key in cluster_param]))

os.system("sbatch --tasks-per-node {threads} --mem={mem}G --job-name={name} {script}".format(script=jobscript, **cluster_param))
