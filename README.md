# ONeSAMP 3.0 ![Latest Release](https://img.shields.io/github/v/release/AaronHong1024/ONeSAMP_3?display_name=tag&sort=semver) ![License](https://img.shields.io/github/license/AaronHong1024/ONeSAMP_3)


<!-- Badges (optional)
[![GitHub release](https://img.shields.io/github/v/release/<owner>/ONeSAMP_3)]()
[![License](https://img.shields.io/github/license/<owner>/ONeSAMP_3)]()
-->

[ONeSAMP 3.0: estimation of effective population size via single nucleotide polymorphism data from one population](https://academic.oup.com/g3journal/article/14/10/jkae153/7712987)

ONeSAMP 3.0 reads a dataset in **GENEPOP** format, computes five summary statistics, and fits a linear regression model to estimate the effective population size. We **strongly recommend** reading the accompanying manuscript before applying ONeSAMP to your data.

---

## Features

- Single-population **Ne** estimation from SNP data  
- Input: **GENEPOP** format  
- Configurable ranges for **Ne**, **θ**, and population size/duration  

---

## Requirements

- **Linux** (tested)
- **Python** ≥ 3.8
- **R** (required by parts of the pipeline)

> Install R through your package manager (e.g., `apt`, `yum`, `brew`) or from r-project.org.

---

## Getting Started

### 1) Clone the repository

```bash
mkdir ONeSAMP
cd ONeSAMP
git clone https://github.com/AaronHong1024/ONeSAMP_3.git
cd ONeSAMP_3
```

### 2) (Optional) Create a Python virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3) Ensure the bundled binary is executable

```bash
chmod u=rwx,go=rx build/OneSamp
```

---

## Input

- **Format:** GENEPOP (standard header, locus names line, followed by per-individual allele codes per locus).
- **Example data:** `exampleData/genePop10Ix30L`

If you’re new to GENEPOP formatting, validate that:
- Locus names are listed once after the header line.
- Each subsequent line contains `SampleID , a1a2 a1a2 ...` (two digits per allele per locus, or standard GENEPOP allele codes).

---

## Usage

### Command

```bash
python main [--s <trials>] [--o <input_file>]
```

- **positional**:
  - `input` — input file name (GENEPOP)

- **options**:

| Flag | Meaning | Type / Range (as implemented) | Notes |
|---|---|---|---|
| `--n` | Treat monomorphic loci as present | flag (default: False) | Include monomorphic loci |
| `--m` | Minimum allele frequency | float in `[0,1]` | Filter rare alleles |
| `--r` | Mutation rate | float in `[0,1]` | Per-locus/per-site depending on model |
| `--lNe` | Lower bound of Ne range | integer `≥ 10` | Search lower limit |
| `--uNe` | Upper bound of Ne range | integer `≤ 500` | Search upper limit |
| `--lT` | Lower bound of θ range | integer `≥ 1` |  |
| `--uT` | Upper bound of θ range | integer `≤ 10` |  |
| `--s` | Number of ONeSAMP trials | integer `1000–50000` | More trials → more stable estimates |
| `--lD` | Lower bound of duration range | integer `≥ 2` | Generations |
| `--uD` | Upper bound of duration range | integer `≤ 8` | Generations |
| `--i` | Missing data (per-individual) | float in `[0,1]` | Proportion allowed |
| `--l` | Missing data (per-locus) | float in `[0,1]` | Proportion allowed |
| `--o` | **Input** file path (GENEPOP) | path | Matches examples below |

> Note: `--o` is used here to provide the input file path (as in the examples below).

### Examples

Run 1,000 trials on the provided example data and save text output:

```bash
python main --s 1000 --o exampleData/genePop10Ix30L > output.txt
```

Specify ranges and filters:

```bash
python main   --o path/to/your_dataset.gen   --s 20000   --lNe 20 --uNe 400   --lT 1  --uT 10   --lD 2  --uD 6   --m 0.02   --i 0.1   --l 0.1
```

---

## Output

- Console (stdout) summary with the **Ne** estimate and the summary statistics used.  
- Redirect to a file if needed, e.g., `> results/your_run.txt`.

---

## Tips & Troubleshooting

- **Permissions**: If you see a “Permission denied” when calling the bundled binary, re-run  
  `chmod u=rwx,go=rx build/OneSamp`.
- **R not found**: Ensure `R` is installed and available in your `PATH` (e.g., `R --version` should work).
- **GENEPOP parsing issues**: Check that locus names and per-individual genotype lines follow the GENEPOP spec and that allele codes are consistent across loci.

---

## Citation

If you use **ONeSAMP 3.0** in your research, please cite:

> Aaron Hong, Rebecca G Cheek, Suhashi Nihara De Silva, Kingshuk Mukherjee, Isha Yooseph, Marco Oliva, Mark Heim, Chris W. Funk, David Tallmon, Christina Boucher, ONeSAMP 3.0: estimation of effective population size via single nucleotide polymorphism data from one population, G3 Genes|Genomes|Genetics, Volume 14, Issue 10, October 2024, jkae153, https://doi.org/10.1093/g3journal/jkae153

---

## License

This project is released under an open-source license. See `LICENSE` in the repository for details.

---

## Acknowledgments

We thank all contributors and users who tested ONeSAMP and provided feedback to improve its robustness and usability.

 
