import os
from fmbiopy.fmtest import *

def test_estimate_seq_error_produces_valid_output(dat, sandbox):
    bam = str(dat['tiny']['bam'][0])
    csv = str(sandbox / "seq_error.csv")
    os.system('samtools index ' + bam)
    os.system(' '.join(["./estimate_seq_error.R", bam, csv]))
    assert os.path.isfile(csv)
