#!/usr/bin/env python3

"""
pattern_sccmec_detector.py

Detect SCCmec regions with:
  1) Regex-based DR/IR detection (all orientations).
  2) BWA MEM alignment for mecA.
  3) (New) BWA MEM alignment for rlmH (for reporting fallback if SCCmec not found).
  4) Require exactly ONE DR and ONE IR on the left boundary (any order) and
     exactly ONE DR and ONE IR on the right boundary (any order), with mecA
     fully contained between those two boundaries.

Outputs:
  - <ref_name>_DR_hits.tsv
  - <ref_name>_IR_hits.tsv
  - <ref_name>_mecA_hits.tsv
  - <ref_name>_SCCmec_candidate.tsv (if a valid region is found)
  - summary_SCCmec_detection.tsv

Usage Example:
  ./pattern_sccmec_detector.py \
    --ref_dir /path/to/reference_fastas \
    --dr_fasta /path/to/DR.fasta \
    --ir_fasta /path/to/IR.fasta \
    --mecA_fasta /path/to/mecA.fasta \
    --rlmH_fasta /path/to/rlmH.fasta \
    --output_dir /path/to/output
"""

import os
import sys
import re
import argparse
import subprocess
from collections import defaultdict

###############################################################################
# 1. Utility Functions
###############################################################################

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(complement)[::-1]

def load_fasta_sequences(fasta_path):
    """
    Load all contigs from a FASTA file.
    Returns a dictionary: {contig_name: sequence}
    """
    sequences = {}
    current_contig = None
    seq_chunks = []
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_contig:
                    sequences[current_contig] = ''.join(seq_chunks).upper()
                current_contig = line[1:].split()[0]  # First word as contig name
                seq_chunks = []
            else:
                seq_chunks.append(line)
        # Add the last contig
        if current_contig:
            sequences[current_contig] = ''.join(seq_chunks).upper()
    return sequences

def find_all_matches(pattern, genome_seq):
    """
    Find all matches of the regex pattern in the genome sequence.
    Returns a list of (start0, end) tuples (0-based, end-exclusive).
    """
    return [(m.start(), m.end()) for m in re.finditer(pattern, genome_seq)]

def write_hits_to_tsv(hits, contig_name, out_tsv, label):
    """
    Append alignment hits to a TSV file with headers.
    Each hit is a tuple: (start0, end)
    """
    with open(out_tsv, 'a') as f:
        for (start0, end) in hits:
            f.write(f"{contig_name}\t{start0}\t{end}\t{label}\n")

def write_summary(summary_path, ref_name, detected, contig="NA", start0="NA", end="NA"):
    """Append a line to the summary TSV."""
    with open(summary_path, 'a') as f:
        f.write(f"{ref_name}\t{detected}\t{contig}\t{start0}\t{end}\n")

def parse_bwa_mem_sam(sam_file):
    """
    Parse SAM file from BWA MEM.
    Returns a list of (contig, start0, end) tuples.
    """
    hits = []
    if not os.path.isfile(sam_file) or os.path.getsize(sam_file) == 0:
        return hits

    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            flag = int(parts[1])
            if (flag & 0x4) != 0:  # Unmapped
                continue
            if (flag & 0x100) != 0:  # Secondary alignment
                continue
            contig = parts[2]
            pos1 = int(parts[3])
            cigar = parts[5]
            aln_length = cigar_ref_length(cigar)
            start0 = pos1 - 1
            end = start0 + aln_length
            hits.append((contig, start0, end))
    return hits

def cigar_ref_length(cigar):
    """
    Calculate the alignment length on the reference from the CIGAR string.
    Sum the lengths of M, D, N, =, X operations.
    """
    length = 0
    matches = re.findall(r'(\d+)([MIDNSHP=X])', cigar)
    for count, op in matches:
        count = int(count)
        if op in ['M', 'D', 'N', '=', 'X']:
            length += count
    return length

def make_bwa_index(ref_fasta):
    """
    Create a BWA index for the reference if not already present.
    Index files: <ref>.amb, <ref>.ann, <ref>.bwt, <ref>.pac, <ref>.sa
    """
    index_exts = [".amb", ".ann", ".bwt", ".pac", ".sa"]
    has_index = all(os.path.exists(ref_fasta + ext) for ext in index_exts)
    if has_index:
        print(f"[BWA Index] Index already exists for {ref_fasta}. Skipping.")
        return
    print(f"[BWA Index] Creating BWA index for {ref_fasta}")
    cmd = ["bwa", "index", ref_fasta]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        print(f"[BWA Index] Successfully created BWA index for {ref_fasta}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] BWA index failed for {ref_fasta}: {e}", file=sys.stderr)
        sys.exit(1)

def run_bwa_mem(ref_fasta, gene_fasta, sam_out):
    """
    Align a gene (e.g., mecA or rlmH) using BWA MEM and write SAM output.
    """
    print(f"[BWA MEM] Aligning {os.path.basename(gene_fasta)} to {os.path.basename(ref_fasta)}")
    cmd = ["bwa", "mem", ref_fasta, gene_fasta]
    try:
        with open(sam_out, 'w') as fout:
            subprocess.run(cmd, check=True, stdout=fout, stderr=subprocess.DEVNULL)
        print(f"[BWA MEM] Alignment complete: {sam_out}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] BWA MEM failed for {ref_fasta}: {e}", file=sys.stderr)
        sys.exit(1)

###############################################################################
# 2. SCCmec Detection Logic
###############################################################################

def detect_minimal_sccmec(dr_hits, ir_hits, mecA_hits):
    """
    Identify a minimal region containing:
      - Exactly 1 DR AND 1 IR on the left boundary (in any order).
      - Exactly 1 DR AND 1 IR on the right boundary (in any order).
      - mecA fully in between these boundaries.

    Pseudocode:
      For each contig:
        for each possible (dr_left, ir_left) combination:
           left_min = min of their coordinates
           left_max = max of their coordinates
           # This pair forms the "left boundary"

           for each possible (dr_right, ir_right) combination:
              right_min = min of their coordinates
              right_max = max of their coordinates
              # This pair forms the "right boundary"

              # Must ensure right_min > left_max (so there's space in between)
              if right_min <= left_max:
                  continue

              # Then check each mecA hit:
              # We want entire mecA to be within (left_max, right_min)
              #    => mecA_start >= left_max AND mecA_end <= right_min

    Return (contig, bounding_min, bounding_max) if found, else None.
    """

    dr_by_contig = defaultdict(list)
    for (ctg, s, e, _orient) in dr_hits:
        dr_by_contig[ctg].append((s, e))

    ir_by_contig = defaultdict(list)
    for (ctg, s, e, _orient) in ir_hits:
        ir_by_contig[ctg].append((s, e))

    mecA_by_contig = defaultdict(list)
    for (ctg, s, e) in mecA_hits:
        mecA_by_contig[ctg].append((s, e))

    # Sort them by start coordinate (for consistency)
    for c in dr_by_contig:
        dr_by_contig[c].sort(key=lambda x: x[0])
    for c in ir_by_contig:
        ir_by_contig[c].sort(key=lambda x: x[0])
    for c in mecA_by_contig:
        mecA_by_contig[c].sort(key=lambda x: x[0])

    candidate_contigs = set(dr_by_contig.keys()) & set(ir_by_contig.keys()) & set(mecA_by_contig.keys())

    minimal_region = None
    minimal_length = None

    for contig in candidate_contigs:
        dr_list = dr_by_contig[contig]
        ir_list = ir_by_contig[contig]
        mecA_list = mecA_by_contig[contig]

        # For each pair DR_left, IR_left => left boundary
        for (dr_s1, dr_e1) in dr_list:
            for (ir_s1, ir_e1) in ir_list:
                left_min = min(dr_s1, dr_e1, ir_s1, ir_e1)
                left_max = max(dr_s1, dr_e1, ir_s1, ir_e1)

                # For each pair DR_right, IR_right => right boundary
                for (dr_s2, dr_e2) in dr_list:
                    for (ir_s2, ir_e2) in ir_list:
                        right_min = min(dr_s2, dr_e2, ir_s2, ir_e2)
                        right_max = max(dr_s2, dr_e2, ir_s2, ir_e2)

                        # Must ensure the right boundary is to the right of the left boundary
                        if right_min <= left_max:
                            continue

                        # Check mecA in between
                        for (ma_s, ma_e) in mecA_list:
                            # We want the entire mecA to be between left_max and right_min
                            if ma_s >= left_max and ma_e <= right_min:
                                # We have a valid region
                                bounding_min = min(left_min, right_min)
                                bounding_max = max(left_max, right_max)
                                length = bounding_max - bounding_min

                                if minimal_length is None or length < minimal_length:
                                    minimal_length = length
                                    minimal_region = (contig, bounding_min, bounding_max)

    return minimal_region

###############################################################################
# 3. Main Script
###############################################################################

def main():
    parser = argparse.ArgumentParser(
        description="Comprehensive Pattern-based SCCmec detection: Regex for DR and IR (all orientations), BWA MEM for mecA, plus fallback to rlmH if SCCmec not found.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example Usage:
  ./pattern_sccmec_detector.py \\
    --ref_dir /path/to/reference_fastas \\
    --dr_fasta /path/to/DR.fasta \\
    --ir_fasta /path/to/IR.fasta \\
    --mecA_fasta /path/to/mecA.fasta \\
    --rlmH_fasta /path/to/rlmH.fasta \\
    --output_dir /path/to/output
"""
    )
    parser.add_argument("--ref_dir", required=True,
                        help="Directory containing reference FASTA files (.fna, .fas, .fasta).")
    parser.add_argument("--dr_fasta", required=True,
                        help="FASTA file containing the DR sequence.")
    parser.add_argument("--ir_fasta", required=True,
                        help="FASTA file containing the IR sequence.")
    parser.add_argument("--mecA_fasta", required=True,
                        help="FASTA file containing the mecA gene sequence.")
    parser.add_argument("--rlmH_fasta", required=True,
                        help="FASTA file containing the rlmH gene sequence.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to store all output files.")

    args = parser.parse_args()

    ref_dir = args.ref_dir
    dr_fasta = args.dr_fasta
    ir_fasta = args.ir_fasta
    mecA_fasta = args.mecA_fasta
    rlmH_fasta = args.rlmH_fasta
    output_dir = args.output_dir

    os.makedirs(output_dir, exist_ok=True)

    # Initialize summary file
    summary_file = os.path.join(output_dir, "summary_SCCmec_detection.tsv")
    with open(summary_file, 'w') as sf:
        sf.write("Reference\tSCCmec_Detected\tContig\tStart0\tEnd\n")

    #===========================#
    # 1) Prepare DR Regex
    #===========================#
    dr_seqs = load_fasta_sequences(dr_fasta)
    if not dr_seqs:
        print(f"[ERROR] No sequences found in DR FASTA: {dr_fasta}", file=sys.stderr)
        sys.exit(1)
    _, dr_seq = next(iter(dr_seqs.items()))

    dr_complements = {
        "DR": dr_seq,
        "DR_reverse": dr_seq[::-1],
        "DR_complement": reverse_complement(dr_seq),
        "DR_reverse_complement": reverse_complement(dr_seq)[::-1]
    }
    dr_patterns = {}
    for key, seq in dr_complements.items():
        pattern = re.escape(seq).replace("N", "[ACGT]")
        dr_patterns[key] = pattern

    #===========================#
    # 2) Prepare IR Regex
    #===========================#
    ir_seqs = load_fasta_sequences(ir_fasta)
    if not ir_seqs:
        print(f"[ERROR] No sequences found in IR FASTA: {ir_fasta}", file=sys.stderr)
        sys.exit(1)
    _, ir_seq = next(iter(ir_seqs.items()))

    ir_complements = {
        "IR": ir_seq,
        "IR_reverse": ir_seq[::-1],
        "IR_complement": reverse_complement(ir_seq),
        "IR_reverse_complement": reverse_complement(ir_seq)[::-1]
    }
    ir_patterns = {}
    for key, seq in ir_complements.items():
        pattern = re.escape(seq).replace("N", "[ACGT]")
        ir_patterns[key] = pattern

    #===========================#
    # 3) Gather Reference Files
    #===========================#
    ref_files = [
        os.path.join(ref_dir, f) for f in os.listdir(ref_dir)
        if any(f.endswith(ext) for ext in [".fna", ".fas", ".fasta"])
    ]
    if not ref_files:
        print(f"[ERROR] No reference FASTA files found in {ref_dir}. Exiting.", file=sys.stderr)
        sys.exit(1)

    #===========================#
    # 4) Process Each Reference
    #===========================#
    for ref_fasta in sorted(ref_files):
        ref_name = os.path.basename(ref_fasta)
        base_name = os.path.splitext(ref_name)[0]
        print(f"\n[PROCESSING] {ref_name}")

        #--------------------#
        # Make BWA Index
        #--------------------#
        make_bwa_index(ref_fasta)

        #--------------------#
        # Load Contigs
        #--------------------#
        contigs = load_fasta_sequences(ref_fasta)
        if not contigs:
            print(f"[WARNING] No contigs found in {ref_fasta}. Skipping.", file=sys.stderr)
            continue

        #--------------------#
        # Regex DR/IR Search
        #--------------------#
        dr_hits = []
        ir_hits = []
        for contig_name, seq in contigs.items():
            # DR
            for orientation, pattern in dr_patterns.items():
                matches = find_all_matches(pattern, seq)
                for (start0, end) in matches:
                    dr_hits.append((contig_name, start0, end, orientation))

            # IR
            for orientation, pattern in ir_patterns.items():
                matches = find_all_matches(pattern, seq)
                for (start0, end) in matches:
                    ir_hits.append((contig_name, start0, end, orientation))

        # Write DR hits
        dr_tsv = os.path.join(output_dir, f"{base_name}_DR_hits.tsv")
        with open(dr_tsv, 'w') as f:
            f.write("contig\tstart0\tend\torientation\n")
        for (contig_name, s0, e0, orient) in dr_hits:
            write_hits_to_tsv([(s0, e0)], contig_name, dr_tsv, orient)

        # Write IR hits
        ir_tsv = os.path.join(output_dir, f"{base_name}_IR_hits.tsv")
        with open(ir_tsv, 'w') as f:
            f.write("contig\tstart0\tend\torientation\n")
        for (contig_name, s0, e0, orient) in ir_hits:
            write_hits_to_tsv([(s0, e0)], contig_name, ir_tsv, orient)

        #--------------------#
        # BWA MEM for mecA
        #--------------------#
        mecA_sam = os.path.join(output_dir, f"{base_name}_mecA.sam")
        run_bwa_mem(ref_fasta, mecA_fasta, mecA_sam)
        mecA_hits = parse_bwa_mem_sam(mecA_sam)

        mecA_tsv = os.path.join(output_dir, f"{base_name}_mecA_hits.tsv")
        with open(mecA_tsv, 'w') as f:
            f.write("contig\tstart0\tend\n")
            for (contig, st, en) in mecA_hits:
                f.write(f"{contig}\t{st}\t{en}\n")

        #--------------------#
        # BWA MEM for rlmH
        #--------------------#
        rlmH_sam = os.path.join(output_dir, f"{base_name}_rlmH.sam")
        run_bwa_mem(ref_fasta, rlmH_fasta, rlmH_sam)
        rlmH_all_hits = parse_bwa_mem_sam(rlmH_sam)
        # We only need one "best" rlmH hit or the first one if multiple
        # (you can implement more sophisticated best-hit selection if needed)
        if rlmH_all_hits:
            # Example: pick the alignment with the largest alignment length
            # or just pick the first.  Here we pick the largest alignment:
            rlmH_hit = max(rlmH_all_hits, key=lambda x: (x[2] - x[1]))
            rlmH_contig, rlmH_start, rlmH_end = rlmH_hit
        else:
            rlmH_contig, rlmH_start, rlmH_end = ("NA", "NA", "NA")

        #--------------------#
        # Write rlmH hits TSV
        #--------------------#
        rlmH_tsv = os.path.join(output_dir, f"{base_name}_rlmH_hits.tsv")
        with open(rlmH_tsv, 'w') as f:
            f.write("contig\tstart0\tend\n")
            for (contig, st, en) in rlmH_all_hits:
                f.write(f"{contig}\t{st}\t{en}\n")

        #--------------------#
        # Detect SCCmec
        #--------------------#
        minimal_region = detect_minimal_sccmec(dr_hits, ir_hits, mecA_hits)

        if minimal_region:
            # SCCmec found
            ctg, min_coord, max_coord = minimal_region
            print(f"[SCCmec Detected] {ref_name}: {ctg}:{min_coord}-{max_coord}")
            sccmec_tsv = os.path.join(output_dir, f"{base_name}_SCCmec_candidate.tsv")
            with open(sccmec_tsv, 'w') as f:
                f.write("contig\tstart0\tend\n")
                f.write(f"{ctg}\t{min_coord}\t{max_coord}\n")

            # Write summary line with SCCmec region
            write_summary(summary_file, ref_name, "Yes", ctg, min_coord, max_coord)

        else:
            # No SCCmec; fallback to reporting rlmH in summary
            print(f"[No SCCmec Detected] {ref_name}")
            if rlmH_contig == "NA":
                # If we didn't find rlmH either, just output "No NA NA NA"
                write_summary(summary_file, ref_name, "No", "NA", "NA", "NA")
            else:
                write_summary(summary_file, ref_name, "No",
                              rlmH_contig, rlmH_start, rlmH_end)

if __name__ == "__main__":
    main()
