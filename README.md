# vdjmatch2

`vdjmatch2` is a CLI tool for pairwise repertoire matching.

The tool compares two clonotype tables in TSV format:
- the first table is treated as the query repertoire
- the second table is treated as the target repertoire
- for each query clonotype, all target clonotypes within the specified search radius are reported

The result file contains information about both the query clonotype and the matched target clonotype. The same target clonotype may therefore appear multiple times if it matches multiple queries.

## CLI arguments

| Flag             |        Type | Default                 | Description                                                                                                |
|------------------|------------:|-------------------------|------------------------------------------------------------------------------------------------------------|
| `query_tsv`      |  positional | —                       | First input clonotype table (query repertoire, TSV).                                                       |
| `target_tsv`     |  positional | —                       | Second input clonotype table (target repertoire, TSV).                                                     |
| `--out`          |       `str` | `match_result.tsv`      | Output TSV path.                                                                                           |
| `--max-sub`      |       `int` | `1`                     | Maximum substitutions allowed.                                                                             |
| `--max-ins`      |       `int` | `0`                     | Maximum insertions allowed.                                                                                |
| `--max-del`      |       `int` | `0`                     | Maximum deletions allowed.                                                                                 |
| `--max-edits`    |       `int` | `sum(sub + ins + del)`  | Maximum total number of edit operations allowed.                                                           |
| `--matrix-path`  |       `str` | —                       | Path to a substitution matrix file. If set, matrix search mode is used.                                    |
| `--max-cost`     | `int/float` | `6`                     | Maximum total alignment cost allowed in matrix mode.                                                       |
| `--match-v`      |      `flag` | off                     | Require V gene to match when counting a match.                                                             |
| `--match-j`      |      `flag` | off                     | Require J gene to match when counting a match.                                                             |
| `--gene`         |       `str` | —                       | Chain filter applied while reading input tables.                                                           |
| `--species`      |       `str` | `HomoSapiens`           | Species filter applied while reading input tables.                                                         |
| `--epitope`      |       `str` | —                       | Epitope filter applied while reading input tables. If the epitope column is absent, the filter is ignored. |
| `--threads`      |       `int` | `4`                     | Maximum number of worker threads for search.                                                               |
| `--junction-col` |       `str` | `junction_aa or cdr3`   | Column name for the CDR3 amino-acid sequence. If not specified, the reader tries common alternatives.      |
| `--v-col`        |       `str` | `v_call or v.segm`      | Column name for the V gene call. If not specified, the reader tries common alternatives.                   |
| `--j-col`        |       `str` | `j_call or j.segm`      | Column name for the J gene call. If not specified, the reader tries common alternatives.                   |
| `--epitope-col`  |       `str` | `antigen.epitope`       | Column name for the epitope. This column is optional.                                                      |
| `--species-col`  |       `str` | `species`               | Column name for the species label.                                                                         |
| `--chain-col`    |       `str` | `gene`                  | Column name for the chain label, for example `TRB` or `TRA`.                                               |
| `--align`        |      `flag` | off                     | Add alignment information to the output.                                                                   |

## Search modes

`vdjmatch2` supports two search modes.

### 1. Bounded edit mode

This is the default mode if `--matrix-path` is not provided.

The search is controlled by:
- `--max-sub`
- `--max-ins`
- `--max-del`
- `--max-edits`

Example:

```bash
vdjmatch2 query.tsv target.tsv \
  --out match_result.tsv \
  --max-sub 1 \
  --max-ins 0 \
  --max-del 0 \
  --max-edits 1 \
  --threads 4 \
  --match-v \
  --match-j
```

### 2. Substitution-matrix mode

If `--matrix-path` is provided, the tool switches to matrix-based matching.

The search is controlled by:
- `--matrix-path`
- `--max-cost`

Example:

```bash
vdjmatch2 query.tsv target.tsv \
  --matrix-path blosum62.txt \
  --max-cost 6 \
  --threads 4 \
  --align
```

## Installation

#### Requirements
- Python 3.8+
- CMake
- C++17-compatible compiler

### 1. Install from PyPI

```bash
pip install vdjmatch2
```

### 2. Install directly from GitHub


```bash
pip install "git+https://github.com/antigenomics/vdjmatch2.git"
```

### 3. Clone the repository and install locally

```bash
git clone https://github.com/antigenomics/vdjmatch2.git
cd vdjmatch2
pip install .
```

After installation, the CLI should be available as:

```bash
vdjmatch2 --help
```

If installation fails, check that:
- `cmake` is available in your shell
- a working C++ compiler is installed
- Python and `pip` are available in the environment where you run the installation


## Input and output format

`vdjmatch2` takes **two input TSV files**: a **query repertoire** and a **target repertoire**.  
Each row is a single clonotype. The tool iterates over clonotypes from the first file and searches for matches in the second one.

Both files must be **tab-separated** and contain a **header row**.

### Input columns

Only the CDR3 amino-acid sequence column is required. By default, `vdjmatch2` tries to detect it as:

- `junction_aa`
- `cdr3`

You can also set it explicitly with `--junction-col`.

Other columns are optional and are used only when needed:

- V gene: `--v-col` (*v_call* or *v.segm*)
- J gene: `--j-col` (*j_call* or *j.segm*)
- epitope: `--epitope-col` (*antigen.epitope*)
- species: `--species-col` (*species*)
- chain: `--chain-col` (*gene*)

If a filter is requested, the corresponding column must be present:

- `--match-v` → V column
- `--match-j` → J column
- `--epitope` → epitope column
- `--species` → species column
- `--gene` → chain column

Dataset-level filters (`--gene`, `--species`, `--epitope`) are applied **before** trie construction and **before** matching.

### Output file

The output is a TSV where **each row is one query–target match**.  
If several target clonotypes fall into the allowed radius for one query clonotype, all of them are written as separate rows. The same target clonotype may therefore appear multiple times.

The output contains fields from both the query and target clonotypes, plus match metadata.

Depending on the mode, it also includes:

- **edit mode**: distance, substitutions, insertions, deletions
- **matrix mode**: match cost
- **optional alignment mode**: query alignment, target alignment

So the output is **not** a deduplicated repertoire intersection, but a full list of all matched query–target pairs.
