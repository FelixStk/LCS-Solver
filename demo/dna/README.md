# Demo: dna
This demo can be used to benchmark the library using the [afproject genetree dataset](https://afproject.org/media/genetree/swisstree/dataset/swisstree.zip)

## Usage
```text
Usage: dna [-p ProblemType] [-t AlgoType] [-f file1.fasta file2.fasta] [-v] [-h]
           [-s result.tsv] [-V] [-g gc.json] [-a alph.json [-l left.json | -r 
           right.json | both]] [-P path/to/swisstree]
Options:
-P path/to/swisstree : Path to SwissTree dataset.
-V : Display version information.
-a alph.json : Enable symbol-dependent constraints. Requires -l/-r.
-f file1.fasta file2.fasta : Input FASTA files.
-g gc.json : JSON file with gap constraints for MC algorithms.
-h : Show this help message.
-l left.json : Left gap constraint tuple (sorted by alph.json).
-p ProblemType : Specifies the problem's type: LCS_Classic LCS2_MC LCS2_MC_INC 
                 LCS2_MC_1C LCS2_MC_O1C_SYNC LCS2_Sigma_R LCS2_Sigma_L LCS2_SIGMA
-r right.json : Right gap constraints tuple (sorted by alph.json).
-s result.tsv : Save results to TSV file.
-t AlgoType : Specifies the algorithm: LLCS2_STD_FL LLCS2_MC LLCS2_MC_INC 
              LLCS2_MC_INC_E LLCS2_MC_1C LLCS2_MC_O1_SYNC LLCS2_SR_MQ 
              LLCS2_SR_RMQ LLCS2_SL_R_LLCS2_SR_MQ LLCS2_SL_R_LLCS2_SR_RMQ 
              LLCS2_SA_MQ LLCS2_SA_RMQ
-v : Enable verbose output (shows a progress bar)
```

## Notes
- After building the project, a copy of the dataset will be placed in the `swisstree` directory next to the `dna` executable. This means the program can be run with the `-P ./swisstree` flag.
- If the `-f` flag is not used, the algorithm specified by the `-t` flag will be executed on all combinations of files in the dataset directory.
- The JSON files must be correctly formatted. Examples:
    - `alph.json` contains the alphabet of the input sequences: `["A", "B", "C"]`
    - `left.json` defines the left mapping based on the order of symbols in `alph.json`: `[[0, 9], [0, 9], [0, 9]]`
    - `right.json` defines the right mapping based on the order of symbols in `alph.json`: `[[0, 9], [0, 9], [0, 9]]`
    - `gc.json` defines a gap tuple. For example: `[[0, 9], [0, 9], [0, 9]]`
- If the number of pairs in a JSON file does not match the length of the alphabet or the length of the shorter input sequence - 1, the constraint will be padded with relaxed gap constraint pairs of the form `(0, length of the longest sequence)`.
- If the `-t` flag is not provided, the program will default to `-t LLCS2_STD_FL`.
- If the `-s result.tsv` flag is not provided, the program will output the LLCS solutions for all file combinations to the standard output.
