# pitaya

This is a collection of script associated with the re-analysis of a pitaya RNA-Seq dataset.


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





python calculate_RPKMs.py
					--counts <COUNT_TABLE>
					--rpkm <RPKM_FILE>
					--assembly <ASSEMBLY_FILE>


The script pitaya_MYB_exp_plots_summary.py and pitaya_exp_plots_summary.py contain dataset specific elements. These scripts are included for documentation purposes, but should not be used on other datasets without adjustments.

## References

