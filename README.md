# SCCmecScanner

**Rapid extraction of SCC<sub>mec</sub> boundaries from *Staphylococcus* genome assemblies**


SCCmecScanner is a lightweight Python utility that pinpoints the **start and end coordinates of the staphylococcal cassette chromosome mec (SCC<sub>mec</sub>)** element in draft or complete *Staphylococcus* genome FASTA files.

I am just sharing the code to make people not to suffer the same problem.

> ⚠️ **Scope** This tool **does not classify** SCC<sub>mec</sub> into types or subtypes. If you want, please use the companion repository **[SCCmec\_classifier](https://github.com/jeju2486/sccmec_classifier)**.



## Installation

> Tested on **Python ≥ 3.8** under Linux and macOS. Windows should work but is not routinely tested.

```bash
# 1. Clone the repository
$ git clone https://github.com/jeju2486/SCCmecScanner.git
$ cd SCCmecScanner

# 2. Install dependencies
$ pip install -r requirements.txt
```

Dependencies are minimal—principally **Biopython** and the Python standard library.  See `requirements.txt` for exact versions.


## Quick start

```bash
python sccmec_finder.py \
    --ref_dir /path/to/assemblies \
    --out results.tsv
```

| Argument    | Required | Description                                                               |
| ----------- | -------- | ------------------------------------------------------------------------- |
| `--ref_dir` | ✅        | Folder containing genome assemblies (`*.fa`, `*.fasta`, `*.fna`, `*.fas`) |
| `--out` |          | Path to write tab‑delimited results (default: `sccmec_coordinates.tsv`)   |

The output file lists, for each assembly, the contig ID plus the 1‑based **start** and **end** nucleotide coordinates of the SCC<sub>mec</sub> element.


## How it works — algorithm overview

The implementation reproduces the decision rules described by Farhat *et al.* (2021) \[[PubMed 34370581](https://pubmed.ncbi.nlm.nih.gov/34370581/)]. Briefly:

1. **Detect Direct/Inverted Repeats** Search the direct repeats (DR) and inverted repeats (IR) within the genome 
2. **Save the locational information of DR and IRs** Save the location information
3. **Detect the mecA sequence** Search the mecA sequence
4. **Make the proper combination of DR/IR/mecA** Make sure there is a combination of DR-IR-mecA-IR (but inverted)-DR and save the location info 

For full details, consult the original publication or read the commented source code (`sccmec_finder.py`).


## Limitations 

* Draft assemblies with highly fragmented SCC<sub>mec</sub> may lead to truncated coordinates.
* Non‑canonical SCC elements (e.g., SCC<sup>fus</sup>) are **not** currently detected.

## Contact

I don't think I will update this code. It is just for sharing. But if you have any problem feel free to post it **Issue** tab
