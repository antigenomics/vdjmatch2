# vdjmatch2

`vdjmatch2` is a CLI tool for pairwise repertoire matching.

The tool compares two clonotype tables in TSV format:
- the first table is treated as the query repertoire
- the second table is treated as the target repertoire
- for each query clonotype, all target clonotypes within the specified search radius are reported

The result file contains information about both the query clonotype and the matched target clonotype. The same target clonotype may therefore appear multiple times if it matches multiple queries.

## CLI arguments

| Flag             |        Type | Default               | Description                                                                                                |
|------------------|------------:|-----------------------|------------------------------------------------------------------------------------------------------------|
| `input_tsv`      | positional  | —                     | First input clonotype table (query repertoire, TSV).                                                       |
| `input_tsv`      |  positional | —                     | Second input clonotype table (target repertoire, TSV).                                                     |
| `--out`          |       `str` | `match_result.tsv`    | Output TSV path.                                                                                           |
| `--max-sub`      |       `int` | `1`                   | Maximum substitutions allowed.                                                                             |
| `--max-ins`      |       `int` | `0`                   | Maximum insertions allowed.                                                                                |
| `--max-del`      |       `int` | `0`                   | Maximum deletions allowed.                                                                                 |
| `--max-edits`    |       `int` | `1`                   | Maximum total number of edit operations allowed.                                                           |
| `--matrix-path`  |       `str` | —                     | Path to a substitution matrix file. If set, matrix search mode is used.                                    |
| `--max-cost`     | `int/float` | `6`                   | Maximum total alignment cost allowed in matrix mode.                                                       |
| `--match-v`      |      `flag` | off                   | Require V gene to match when counting a match.                                                             |
| `--match-j`      |      `flag` | off                   | Require J gene to match when counting a match.                                                             |
| `--gene`         |       `str` | —                     | Chain filter applied while reading input tables.                                                           |
| `--species`      |       `str` | `HomoSapiens`         | Species filter applied while reading input tables.                                                         |
| `--epitope`      |       `str` | —                     | Epitope filter applied while reading input tables. If the epitope column is absent, the filter is ignored. |
| `--threads`      |       `int` | `4`                   | Maximum number of worker threads for search.                                                               |
| `--junction-col` |       `str` | `junction_aa or cdr3` | Column name for the CDR3 amino-acid sequence. If not specified, the reader tries common alternatives.      |
| `--v-col`        |       `str` | `v_call or v.segm`    | Column name for the V gene call. If not specified, the reader tries common alternatives.                   |
| `--j-col`        |       `str` | `j_call or j.segm`    | Column name for the J gene call. If not specified, the reader tries common alternatives.                   |
| `--epitope-col`  |       `str` | `antigen.epitope`     | Column name for the epitope. This column is optional.                                                      |
| `--species-col`  |       `str` | `species`             | Column name for the species label.                                                                         |
| `--chain-col`    |       `str` | `gene`                | Column name for the chain label, for example `TRB` or `TRA`.                                               |
| `--align`        |      `flag` | off                   | Add alignment information to the output.                                                                   |

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
./build/vdjmatch2 query.tsv target.tsv \
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
./build/vdjmatch2 query.tsv target.tsv \
  --matrix-path blosum62.txt \
  --max-cost 6 \
  --threads 4 \
  --align
```

## Output

The output is a TSV file with one row per match.

Each row contains:
- fields from the query repertoire entry
- fields from the matched target repertoire entry
- search metrics such as edit distance or matrix score
- optional alignment strings if `--align` is enabled

Because the tool reports all matches within the search radius, one target clonotype may appear in multiple rows.