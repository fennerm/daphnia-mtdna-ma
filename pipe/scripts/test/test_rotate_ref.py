from plumbum import local
from pytest import fixture

@fixture(name="rotate_ref")
def gen_rotate_ref():
    return local['./rotate_ref.py']

def test_rotseq1(rotate_ref):
    rotseq = "test/dat/rotseq1.fa"
    rotated = rotate_ref(rotseq).split("\n")
    assert rotated[0].startswith(">4:5_rot_to_start")
    assert rotated[1] == "GGAAA"


def test_rotseq2(rotate_ref):
    rotseq = "test/dat/rotseq2.fa"
    rotated = rotate_ref(rotseq).split("\n")
    assert rotated[0].startswith(">4:6_rot_to_start")
    assert rotated[1] == "GGGAAA"


def test_rotseq3(rotate_ref):
    rotseq = "test/dat/rotseq3.fa"
    rotated = rotate_ref(rotseq).split("\n")
    assert rotated[0].startswith(">5:7_rot_to_start")
    assert rotated[1] == "GGGAAAA"


def test_rotseq4(rotate_ref):
    rotseq = "test/dat/rotseq4.fa"
    rotated = rotate_ref(rotseq).split("\n")
    assert rotated[0].startswith(">5:8_rot_to_start")
    assert rotated[1] == "GGGGAAAA"
