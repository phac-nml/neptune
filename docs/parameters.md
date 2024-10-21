# Parameters #

A help message may be viewed by running:


```bash
neptune --help
```

## Mandatory ##

Neptune requires the location of the inclusion, exclusion, and output directories. The remaining parameters will be estimated based on the input sequence or revert to default settings. The following is the minimum number of command line parameters required to run Neptune:

```bash
neptune
    --inclusion /path/to/inclusion/directory/
    --exclusion /path/to/exclusion/directory/
    --output /path/to/output/directory/
```

The following parameters are required by Neptune:

| Option | Alternative | Parameter | Description |
|--------|-------------|-----------|-------------|
| -i | --inclusion | FASTA | A list of inclusion targets in FASTA format. You may list multiple file or directory locations following the parameter. Neptune will automatically include all files within directories. However, Neptune will not recurse into additional directories. |
| -e | --exclusion | FASTA | A list of exclusion targets in FASTA format. You may list multiple file or directory locations following the parameter. Neptune will automatically include all files within directories. However, Neptune will not recurse into additional directories. |
| -o | --output | directory | The location of the output directory. If this directory exists, any files produced with existing names will be overwritten. If this directory does not exist, then it will be created. |

## Optional ##

The optional parameters will either be automatically calculated or be assigned default values.

### *k*-mer ###

The following parameters relate to *k*-mer generation and aggregation:

| Option | Alternative | Parameter | Description |
|--------|-------------|-----------|-------------|
| -k | --kmer | integer | The size of the *k*-mers. This must be a positive integer and should be large enough such that random intra-genome *k*-mer matches, within the largest genome, are unexpected. The size of *k*-mers cannot be larger than the smallest sequence record. This will be automatically calculated if not specified. |
| | --organization | integer | The degree of organization of *k*-mer counting and aggregation. This parameter determines the number nucleotide bases used in parallelized *k*-mer counting and, in turn, the number of parallel instances of *k*-mer aggregation. The number of parallel instances is determined by 4^n, where n is the specified organization parameter. This value must be a non-negative integer smaller than *k*. If the parameter is not specified, then n = 0 and there will be no parallel *k*-mer aggregation. This will likely require a much longer computation time to complete *k*-mer aggregation. |

### Filtering ###

The following command-line parameters relate to signature filtering:

| Option | Alternative | Parameter | Description |
|--------|-------------|-----------|-------------|
| | --filter-length | float | The minimum percent length of a signature candidate against a exclusion target required to filter out the candidate. This value is a percentage expressed as a floating point number [0.0, 1.0]. If the any exclusion hit exceeds the percent length **and** percent identity of any candidate, the candidate is removed. The default value is 0.5. |
| | --filter-percent | float | The minimum percent identity of a signature candidate against a exclusion target required to filter out the candidate. The percent identity is calculated as identities divided by the alignment length. This value is a percentage expressed as a floating point number [0.0, 1.0]. If the any exclusion hit exceeds the percent length **and** percent identity of any candidate, the candidate is removed. The default value is 0.5. |
| | --seed-size | integer | The seed size used for alignments. This value must be no smaller than 4. The default value is 11. |
  
### Extraction ###

The following command-line parameters relate to signature extraction:

| Option | Alternative | Parameter | Description |
|--------|-------------|-----------|-------------|
| -r | --reference | FASTA | A list of references from which to extract signatures. If this parameter is not specified, signatures will be extracted from **all** inclusion targets. You may list multiple file locations following the --reference parameter. |
| | --rate | float | The probability (0.0, 1.0) that any two homologous bases are different from each other. This should incorporate mutation rates, sequencing error rates, and assembly error rates. The rate is used to calculate the maximum allowable gap size in a signature and the minimum expected number of exact *k*-mer matches in a signature. If this value is not specified, the rate is assumed to be 0.01. |
| | --gc-content | float | The expected GC-content of the environment. The GC-content is used to calculate the maximum allowable gap size in a signature and the minimum expected number of exact *k*-mer matches in a signature. If this value is not specified, it is calculated by observing the GC-content of each target during signature extraction. The value must be between (0.0, 1.0). |
| | --confidence | float | The statistical confidence of decision making in the software. The confidence affects the automatic calculation of both the maximum gap size and minimum number of inclusion hits. If this value is not specified, a default of 0.95 is used. The value must be between (0.0, 1.0). |
| | --inhits | integer | The minimum number of inclusion hits required to start and continue signature extraction. If this value is not specified, it will be automatically calculated using the number of inclusion targets, the GC-content, the rate, and the *k*-mer size. The calculation can be found in the *Mathematics* documentation. This value must be a positive integer. |
| | --exhits | integer | The minimum number of exclusion hits necessary to stop extraction of a signature. If this value is not specified, it is assumed to be 1. This value must be a positive integer. |
| | --gap | int | The maximum allowable number of base positions shifted before seeing an exact *k*-mer match. If this value is not specified, it will be automatically calculated using the rate, GC-content, and the *k*-mer size. The calculation can be found in the *Mathematics* documentation. This value must be a positive integer. |
| | --size | int | The minimum size for a signature. Signatures which are shorter than this length will not be reported. If this value is not specified, the minimum signature size will be four times the length of the *k*-mer size. It is not recommended to locate signatures smaller than this size, unless application-specific. This value must be a positive integer. |
  
### Parallelization ###

The following parameters relate to the parallelization of Neptune:

| Option | Alternative | Parameter | Description |
|--------|-------------|-----------|-------------|
| -p | --parallelization | integer | The number of parallel working processes to create. This parameter will directly increase the speed of many stages of the software, provided there are sufficient resources available to run the worker process simultaneously. This value must be a positive integer. The default value is 8. |
