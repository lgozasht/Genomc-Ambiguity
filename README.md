# Genomic-Ambiguity

Ambiguous characters in genomic assemblies can render various challenges with regard to phylogenetics and comparative genomics since they can often be interperetted as either a reference or alternate allele. Here, I provide a script that summarizes the number of ambiguous characters and "N"s in each genome within a phylogeny and provides the location of ambiguous characters at sites above a user specified minor allele frequency.

### Options

**-msa**: Path to input multiple sequence alignment file

**-min_MAF**: Minimum minor allele frequency for output sites

### Usage


```
python3 ambiguityAnalysis.py -msa [path to msa file] -min_MAF [(float) minimum minor allele frequency]
```
