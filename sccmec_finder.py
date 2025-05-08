import os
import re
import argparse

##############################################################################
# 1. Configuration: Patterns for DR, IR, and a substring for mecA
##############################################################################
# For demonstration, we define simple regexes for DR & IR (with possible '-')
# and a substring for mecA. Adjust these as needed.
DR_REGEX = re.compile(r"TTATGATACGC([ACGTN]+)TCT", re.IGNORECASE)
IR_REGEX = re.compile(r"GC([ACGTN]+)TATC",        re.IGNORECASE)
MEC_A_SUBSEQ = "ATGGTCAAGAAAA"  # partial or full mecA

##############################################################################
# 2. Helper Functions
##############################################################################
def read_single_contig_fasta(fasta_path):
    """Return the entire genome sequence (single-contig) as a string."""
    seq_parts = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_parts.append(line)
    return "".join(seq_parts)

def find_dr_hits(genome_seq):
    """Return a list of (start, end) for all DR matches."""
    hits = []
    for m in DR_REGEX.finditer(genome_seq):
        hits.append((m.start(), m.end()))
    return hits

def find_ir_hits(genome_seq):
    """Return a list of (start, end) for all IR matches."""
    hits = []
    for m in IR_REGEX.finditer(genome_seq):
        hits.append((m.start(), m.end()))
    return hits

def find_mecA_hits(genome_seq):
    """
    Return a list of (start, end) for all (possibly overlapping)
    occurrences of MEC_A_SUBSEQ in a case-insensitive manner.
    """
    seq_upper = genome_seq.upper()
    subseq_upper = MEC_A_SUBSEQ.upper()
    hits = []
    i = 0
    while True:
        i = seq_upper.find(subseq_upper, i)
        if i == -1:
            break
        hits.append((i, i + len(MEC_A_SUBSEQ)))
        i += 1  # allow overlapping
    return hits

##############################################################################
# 3. Core Logic: For Each mecA, Find Closest DR/IR on Each Side
##############################################################################
def get_closest_earlier(hit_list, mecA_start):
    """
    Among hits in hit_list = [(start, end)], find the one whose 'end'
    is <= mecA_start and is closest to mecA_start (maximal end).
    Return (start, end) or None if none qualifies.
    """
    valid = [ (st, en) for (st, en) in hit_list if en <= mecA_start ]
    if not valid:
        return None
    # we want the DR/IR whose 'end' is largest but still <= mecA_start
    # so sort by 'end' descending, pick the first
    valid.sort(key=lambda x: x[1], reverse=True)
    return valid[0]  # (start, end) with max end

def get_closest_later(hit_list, mecA_end):
    """
    Among hits in hit_list = [(start, end)], find the one whose 'start'
    is >= mecA_end and is closest to mecA_end (minimal start).
    Return (start, end) or None if none qualifies.
    """
    valid = [ (st, en) for (st, en) in hit_list if st >= mecA_end ]
    if not valid:
        return None
    # we want the DR/IR whose 'start' is smallest but still >= mecA_end
    valid.sort(key=lambda x: x[0])
    return valid[0]  # (start, end) with min start

def detect_sccmec_closest(genome_seq):
    """
    For each mecA hit, find the closest DR & IR earlier, and DR & IR later.
    If found, define bounding region from min(earlier_DR.start, earlier_IR.start)
    to max(later_DR.end, later_IR.end). Return a list of bounding regions.
    """
    dr_hits = find_dr_hits(genome_seq)  # list of (start, end)
    ir_hits = find_ir_hits(genome_seq)  # list of (start, end)
    mecA_hits = find_mecA_hits(genome_seq)

    # If no hits, no detection
    if not dr_hits or not ir_hits or not mecA_hits:
        return []

    bounding_regions = []

    # Sort DR/IR for potential faster searching (not strictly needed)
    # We'll do simplest approach: for each mecA, find the 4 hits.
    for (mstart, mend) in mecA_hits:
        # 1) closest earlier DR
        earlier_dr = get_closest_earlier(dr_hits, mstart)
        # 2) closest earlier IR
        earlier_ir = get_closest_earlier(ir_hits, mstart)
        # 3) closest later DR
        later_dr = get_closest_later(dr_hits, mend)
        # 4) closest later IR
        later_ir = get_closest_later(ir_hits, mend)

        if (earlier_dr is not None and
            earlier_ir is not None and
            later_dr   is not None and
            later_ir   is not None):
            # earlier_dr: (edr_st, edr_end)
            # earlier_ir: (eir_st, eir_end)
            # later_dr:   (ldr_st, ldr_end)
            # later_ir:   (lir_st, lir_end)
            edr_st, edr_end = earlier_dr
            eir_st, eir_end = earlier_ir
            ldr_st, ldr_end = later_dr
            lir_st, lir_end = later_ir

            # bounding region
            region_start = min(edr_st, eir_st)
            region_end   = max(ldr_end, lir_end)

            # Optionally check region validity: region_end > region_start, etc.
            if region_end > region_start:
                bounding_regions.append((region_start, region_end))

    return bounding_regions

##############################################################################
# 4. Main (example usage)
##############################################################################
def main():
    parser = argparse.ArgumentParser(
        description="Closest DR/IR on each side of mecA detection"
    )
    parser.add_argument(
        "--ref_dir",
        required=True,
        help="Directory of single-contig reference FASTA files."
    )
    args = parser.parse_args()

    ref_dir = args.ref_dir
    if not os.path.isdir(ref_dir):
        print(f"[ERROR] {ref_dir} not found or not a directory.")
        return

    fas_list = [f for f in os.listdir(ref_dir)
                if any(f.endswith(ext) for ext in [".fna", ".fas", ".fasta"])]

    if not fas_list:
        print(f"[WARNING] No .fna/.fas/.fasta found in {ref_dir}.")
        return

    for fname in sorted(fas_list):
        fpath = os.path.join(ref_dir, fname)
        seq = read_single_contig_fasta(fpath)
        if not seq:
            print(f"[WARNING] Skipping empty or invalid FASTA: {fname}")
            continue

        regions = detect_sccmec_closest(seq)
        if regions:
            # if multiple mecA => multiple bounding regions
            for (rstart, rend) in regions:
                print(f"[SCCmec Candidate] {fname}: region={rstart}-{rend}")
        else:
            print(f"[No SCCmec] {fname}")

if __name__ == "__main__":
    main()
