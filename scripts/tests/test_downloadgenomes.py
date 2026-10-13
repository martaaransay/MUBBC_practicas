##### To execute tests:
##### python -m scripts.tests


import pytest, os
from scripts.pipeline.downloadgenomes import datasets_flags, taxon_acc_flag

# ------------------------ Function to get flags for download genomes -> datasets_flags
def test_flags_rerseq_chromosome_complete():
    flags = datasets_flags(assembly_source = "refseq",
                           assembly_level = "chromosome,complete",
                           exclude_atypical = True, 
                           exclude_multi_isolate = True)
    # asserts checks if the results are OK
    # If not, AssertionError

    assert flags == ["--assembly-source", "refseq",
                     "--assembly-level", "chromosome,complete",
                     "--exclude-atypical",
                     "--exclude-multi-isolate"]


def test_flags_genbank_complete_mag_only():
    flags = datasets_flags(assembly_source = "genbank",
                           assembly_level = "complete",
                           mag = "only")
    assert flags == ["--assembly-source", "genbank",
                     "--assembly-level", "complete",
                     "--mag", "only"]


def test_flags_annotated_and_mag_true():
    flags = datasets_flags(annotated = True,
                           exclude_atypical = True, 
                           exclude_multi_isolate = True, 
                           mag = True
                           )
    # Changing assert
    assert "--annotated" in flags
    assert "--exclude-atypical" in flags
    assert "--exclude-multi-isolate" in flags
    assert "--mag" in flags

def test_flags_without_args():
    flags = datasets_flags()
    assert isinstance(flags, list)


#------------------------ Function to define taxon/accesion ->  taxon_acc_flag
def test_source_novalid():
    with pytest.raises(ValueError):
        taxon_acc_flag(source = "aa", 
                       name = "Escherichia coli")

def test_inputfile_nofile():
    with pytest.raises(FileNotFoundError):
        taxon_acc_flag(source = "inputfile", 
                       name = "prueba.txt")

def test_inputfile_existente(tmp_path):
    path = os.path.join(tmp_path, "prueba.txt")
    with open(path, "w") as f:
        f.write("Escherichia coli\n")
    taxon_acc_flag(source = "inputfile", 
                   name = path)  # no error

# To execute the same test function with distinct args

@pytest.mark.parametrize(
    "source, name",  # args of the function
    [ # each tuple is a execution
        ("taxon", "Escherichia coli"),
        ("accession", "GCF_000001405.40"),
        ("taxon", "Escherichia coli,Klebsiella pneumoniae"), 
        ("accession", "GCF_000001405.40,GCF_000001635.26"),   
    ],
)
def test_valids_taxon_acc_flag(source, name):
    result = taxon_acc_flag(source = source, 
                            name = name)
    assert isinstance(result, list)
    for names in name.split(","):
        assert names in result
    assert source in result


