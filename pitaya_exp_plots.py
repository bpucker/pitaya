### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python pitaya_exp_plots.py
					--genes <GENES_FILE>
					--exp <EXPRESSION_FILE>
					--out <PREFIX_OF_OUTPUT_FILES>
					"""

import os, sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# --- end of imports --- #


def load_genes( gene_file ):
	"""! @brief load genes """
	
	genes = {}
	gene_order = []
	with open( gene_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			if "," in parts[1]:
				genes.update( { parts[0]: parts[1].split(',') } )
			else:
				genes.update( { parts[0]: [ parts[1] ] } )
			gene_order.append( parts[0] )
			line = f.readline()
	return genes, gene_order


def load_exp( gene_file ):
	"""! @brief load expression values """
	
	exp = {}
	with open( gene_file, "r" ) as f:
		headers = f.readline().strip().split('\t')[1:]
		for header in headers:
			exp.update( { header: {} } )
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			values = map( float, parts[1:] )
			for idx, val in enumerate( values ):
				exp[ headers[ idx ] ].update( { parts[0]: val } )
			line = f.readline()
	return exp


def generate_gene_exp_figure( plot_prefix, genes, exp, gene_order, sample_groups, sample_order ):
	"""! @brief generate figure """
	
	cultivar_colors = { "BR":"blue", "FR":"darkorchid", "DH":"deeppink", "BSJ":"grey" }
	counter = 1
	for gene in gene_order:	#ordered by position in pathway
		# --- collection of data for plot --- #
		data_to_plot = []
		labels = []
		positions = []
		for k, genotype in enumerate( sample_order ):	#run over all 3/4 genotypes
			vals_per_sample = []
			for sample in sample_groups[ genotype ]:	#run over replicates
				for g in genes[ gene ]:	#run over all paralogous genes
					vals_per_sample.append( exp[ sample ][ g ] )
			data_to_plot.append( vals_per_sample )
			labels.append( genotype  )
			positions.append( counter )
			counter+= 1
	
		# --- generate plot --- #
		plot_file = plot_prefix + gene + ".pdf"
		fig, ax = plt.subplots()			
		
		for idx, data in enumerate( data_to_plot ):
			cultivar = sample_order[ idx ]
			ax.boxplot( data, positions=[ idx ], widths=0.75, patch_artist=True, showmeans=True, #notch=True,
								boxprops=dict( facecolor=cultivar_colors[ cultivar ], color=cultivar_colors[ cultivar ]),
								capprops=dict(color=cultivar_colors[ cultivar ]),
								whiskerprops=dict(color=cultivar_colors[ cultivar ]),
								flierprops=dict(color=cultivar_colors[ cultivar ], markeredgecolor=cultivar_colors[ cultivar ]),
								medianprops=dict(color="black"),
								meanprops=dict(color="black")
								 )
		
		
		ax.xaxis.set_ticks( range( 3 ) )
		ax.set_xticklabels( labels, fontsize=12 )
		
		ax.text( 1, 1.05*max( [ x for sublist in data_to_plot for x in sublist ] ), gene, ha="center", fontsize=12 )
		
		ax.set_ylabel( "transcript abundance [TPMs]", fontsize=12 )
		ax.set_xlim( -0.5, 3.5 )
		ax.set_ylim( 0, 1.1*max( [ x for sublist in data_to_plot for x in sublist ] ) )
		
		#ax.xaxis.set_tick_params( labelsize=12 )
		
		my_legend = [	mpatches.Patch(color=cultivar_colors[ "BR" ], label='BR'),
									mpatches.Patch(color=cultivar_colors[ "FR" ], label='FR'),
									mpatches.Patch(color=cultivar_colors[ "DH" ], label='DH')#,
									#mpatches.Patch(color=cultivar_colors[ "BSJ" ], label='BSJ')
								]
		ax.legend( handles=my_legend, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.15), fontsize=12 )
		
		fig.savefig( plot_file, dpi=300 )
		plt.close( "all" )


def main( arguments ):
	"""! @brief run generation of plots """
	
	gene_file = arguments[ arguments.index('--genes')+1 ]
	exp_file = arguments[ arguments.index('--exp')+1 ]
	plot_prefix = arguments[ arguments.index('--out')+1 ]


	sample_groups = { 'BR': [ "SRR11190794", "SRR11190793", "SRR11190792" ],
									'DH': [ "SRR11190802", "SRR11190801", "SRR11190798" ],
									'FR': [ "SRR11190797", "SRR11190796", "SRR11190795" ],
									'BSJ': [ "SRR11190791", "SRR11190800", "SRR11190799" ]
								}
	#DH = red; FR=pink, BR=white; BSJ=white
	sample_order = [ "DH", "FR", "BR" ]	#, "BSJ"

	genes, gene_order = load_genes( gene_file )

	exp = load_exp( exp_file )

	generate_gene_exp_figure( plot_prefix, genes, exp, gene_order, sample_groups, sample_order )


if '--genes' in sys.argv and '--exp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
