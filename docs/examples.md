# Examples #

## Basic ##

The following basic example will report all of the signatures that are sufficiently shared by the (FASTA) sequences in the inclusion directory and sufficiently absent from the (FASTA) sequences in the exclusion directory. Neptune will automatically calculate many of the parameters used in this execution.

```bash
neptune
    --inclusion inclusion_directory/
    --exclusion exclusion_directory/
    --output output_directory/
```

The output of immediate interest will be located in the follow file:

    output_directory/consolidated/consolidated.fasta

This file will contain a consolidated list of signatures, sorted by their Neptune score, which is a combined estimate of sensitivity and specificity. The signatures with higher scores, near the top of the file, are considered the most discriminatory signatures.

## Specifying File Locations ##

You may wish to specify particular files used in signature discovery. This may be important when specifying references for signature extraction:

```bash
neptune
    --inclusion inclusion_dir/ in1.fasta in2.fasta
    --exclusion exclusion_dir/ ex1.fasta ex2.fasta
    --reference in1.fasta in2.fasta
    --output output/
```

## DRMAA Parameters ##

It may be necessary to specify DRMAA native specification parameters to accommodate Neptune job scheduling. This example specifies the resources required by all jobs (--default-specification) and further specifies that *k*-mer aggregation jobs (--aggregate-specification) will require more memory. The remaining Neptune parameters are automatically calculated.

```bash
neptune
    --inclusion inclusion/
    --exclusion exclusion/
    --output output/
    --default-specification "-l h_vmem=6G -pe smp 4"
    --aggregate-specification "-l h_vmem=10G -pe smp 4"
```
