# Pitaya transcriptomics scripts

This is a collection of scripts associated with the re-analysis of a pitaya RNA-Seq dataset. Results of this analysis can be found [here](https://doi.org/10.4119/unibi/2946374). These scripts were developed for use with Python 2.7.


### Calculation of transcriptome assembly statistics 

```
python contig_stats.py
--input <FILENAME>
				
optional:
--min_contig_len <INTEGER> [500]
--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
--exp <EXPRESSION_FILE(normalized)>
```        

`--input` specifies a (multiple) FASTA file which contains the contigs of an assembly or any other collection of sequences. Statistics will be calculated for these sequences under the assumption that it is an assembly.

`--min_contig_len` specifies a cutoff to exclude shorter contigs. A new FASTA file will be generated which contains only sequences passing this filter. Default value is 500bp.

`--out` specifies an output folder where all the generated files will be placed.

`--exp` specifies a text file which contains normalized expression data (e.g. TPMs, RPKMs, FPKMs). These values are used to calculate ExNx statistics like E90N50. Therefore, it is crucial that the sequence names in this file match the sequence names in the provided FASTA file.

This version is based on previously developed scripts ([1](https://doi.org/10.1371/journal.pone.0164321), [2](https://doi.org/10.3389/fmolb.2018.00062)).


### RPKM calculation

```
python calculate_RPKMs.py
--counts <COUNT_TABLE>
--rpkm <RPKM_FILE>
--assembly <ASSEMBLY_FILE>
```

`--counts` this text file contains the raw counts in a matrix with samples in columns and sequences (transcripts/genes) in rows.

`--rpkm` this text file is the output. The layout matches the input file, but all values are normalized based on sequencing depth and sequence length (RPKMs).

`--assembly` this FASTA file contains the sequences mentioned in the expression data file. Lengths of these transcript sequences will be used to calculate RPKM values.


### Plotting transcript abundances

The scripts `pitaya_MYB_exp_plots_summary.py` and `pitaya_exp_plots_summary.py` contain dataset specific elements. These scripts are included for documentation purposes and should not be used on other datasets without adjustments.

```
python pitaya_exp_plots_summary.py
--genes <GENES_FILE>
--exp <EXPRESSION_FILE>
--out <FIGURE_OUTPUT_FILE>
```
`--genes` this text file has two TAB-separated columns. The first column specifies a display name for the transcript/gene. The second column is a single sequence ID or a comma-separated list of sequence IDs which are used to connect transcript abundance values to the display name.

`--exp` this file contains the transcript abundance values as described above. Normalized data (TPMs, RPKMs) should be used for the construction of figures.

`--out` specifies the figure output file. The file extension specifies the file type (e.g. PNG, JPG, PDF, SVG).




## References

Pucker B. and Brockington S. (2020). Pitaya transcriptome assemblies and investigation of transcript abundances. Bielefeld University. doi: [10.4119/unibi/2946374](https://doi.org/10.4119/unibi/2946374).

Haak, M., Vinke, S., Keller, W., Droste, J., Rückert, C., Kalinowski, J., & Pucker, B. (2018). High Quality de novo Transcriptome Assembly of Croton tiglium. Frontiers in Molecular Biosciences, 5. doi:[10.3389/fmolb.2018.00062](https://doi.org/10.3389/fmolb.2018.00062).

Pucker, B., Holtgräwe, D., Rosleff Sörensen, T., Stracke, R., Viehöver, P., and Weisshaar, B. (2016). A de novo Genome Sequence Assembly of the Arabidopsis thaliana Accession Niederzenz-1 Displays Presence/Absence Variation and Strong Synteny. PloS-ONE 11:e0164321. doi:[10.1371/journal.pone.0164321](https://doi.org/10.1371/journal.pone.0164321).

