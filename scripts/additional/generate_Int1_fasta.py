from Bio import SeqIO
from Bio.Seq import Seq, MutableSeq
from Bio.SeqRecord import SeqRecord
import re

fasta = "data/Int1/Int1.fa"
record = SeqIO.read(fasta, "fasta")
Int1_sequence = record.seq
rc_Int1_sequence = str(Seq(Int1_sequence).reverse_complement())

# Variants of Pc to search which is in the ref seq
# ONLY Pc (NOT P2 as is in the attI site)
pc_variants_regex = {
    "PcS": ("TTGACA[ACGT]{17}TAAACT"),
    "PcW": ("TGGACA[ACGT]{17}TAAGCT"),
    "PcWTGN-10": ("TGGACA[ACGT]{14}TG[ACGT]TAAGCT"),
    "PcH1": ("TGGACA[ACGT]{17}TAAACT"),
    "PcH1TGN-10": ("TGGACA[ACGT]{14}TG[ACGT]TAAACT"),
    "PcH2": ("TTGACA[ACGT]{17}TAAGCT"),
    "PcH2TGN-10": ("TTGACA[ACGT]{14}TG[ACGT]TAAGCT"),
    "TGGGCA-N14-TGn-TAAGCT": ("TGGGCA[ACGT]{14}TG[ACGT]TAAGCT"),
    "PcSS": ("TTGATA[ACGT]{17}TAAACT"),
    "Pcln42": ("TTGGCA[ACGT]{17}TAAACT"),
    "Pcln116": ("TTGACA[ACGT]{17}TGAACT"),
    "PcPUO": ("TCGACA[ACGT]{17}TAAACT")}

for name_pc_variant, pc_sequence in pc_variants_regex.items():
    # Find the pc sequence into the RC seq of the Integrase
    matches = re.finditer(pc_sequence, rc_Int1_sequence) 
    for match in matches: # For each match, stores the name, the positions and the sequence
        results = {
            "Pc_variant": name_pc_variant,
            "start": int(match.start()),
            "end": int(match.end())
            }

print(f"Variant: {results['Pc_variant']}\n- start pos: {results['start']}\n- end pos: {results['end']}")

# Only 1 hit - PcW

######### -------------------- Sequence with Pc pattern:

# Manually-generated combined pattern of all variants 
combined_pattern = "T[TGC]G[AG][CT]A[ACGT]{17}T[AG]A[AG]CT"
combined_pattern_BLAST = "TBGRYANNNNNNNNNNNNNNNNNTRARCT"

# Change the positions to the pattern
rc_general_Int1_seq = MutableSeq(rc_Int1_sequence)
rc_general_Int1_seq[results["start"]:results["end"]] = combined_pattern_BLAST
# Reverse complement:
general_Int1_seq = Seq(rc_general_Int1_seq).reverse_complement()

# Store results
pattern_record = SeqRecord(seq=general_Int1_seq,
                           id="Int1_Pc_Pattern",
                           description=""
)


SeqIO.write(pattern_record, "data/Int1/Int1_Pc_pattern.fa", "fasta")

######### -------------------- Sequence without Pc region (start of Pc and 100 pb upstream)

rc_100_Int1_sequence = MutableSeq(rc_Int1_sequence)
rc_100_Int1_sequence = rc_100_Int1_sequence[0:results["start"] - 100]
# Reverse complement:
Int1_100_sequence = Seq(rc_100_Int1_sequence).reverse_complement()
# Store results
no_pc_record = SeqRecord(seq=Int1_100_sequence,
                         id="Int1_NO_Pc",
                         description=""
)

SeqIO.write(no_pc_record, "data/Int1/Int1_no_Pc.fa", "fasta")