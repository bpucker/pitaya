### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

__usage__ = """
					python pitaya_exp_plots_summary.py
					--genes <GENES_FILE>
					--exp <EXPRESSION_FILE>
					--out <FIGURE_OUTPUT_FILE>
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


def generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_groups, genotype_order, upper_border, lower_border ):
	"""! @brief generate figure """
	
	cultivar_colors = { "BR":"blue", "FR":"darkorchid", "DH":"deeppink", "BSJ":"grey" }
	
	data_to_plot = {}
	positions = {}
	label_positions = []
	gene_lables = []
	for cultivar in genotype_order:
		data_to_plot.update( { cultivar: [] } )
		positions.update( { cultivar: [] } )
	counter = 1
	for idx, gene in enumerate( gene_order ):	#ordered by position in pathway
		# --- collection of data for plot --- #
		for k, genotype in enumerate( genotype_order ):	#run over all 3/4 genotypes
			vals_per_sample = []
			new_line = []
			for sample in sample_groups[ genotype ]:	#run over replicates
				for g in genes[ gene ]:	#run over all paralogous genes
					try:
						vals_per_sample.append( exp[ sample ][ g ] )
					except KeyError:
						vals_per_sample.append( 0 )
			
			data_to_plot[ genotype ].append( vals_per_sample )
			positions[ genotype ].append( counter )
			counter+= 1
		counter += 1
		label_positions.append( idx*( len( genotype_order )+1 ) + 0.5*len( genotype_order ) + 0.5 )
		if "-" in gene:
			gene_lables.append( "$\it{" + gene.split('-')[0].replace( "_", "'" ) + "}$" + "-" + gene.split('-')[1] )
		else:
			gene_lables.append( "$\it{" + gene.replace( "_", "'" ) + "}$" )
	
	# --- generate plot --- #
	fig, ( ax, ax2 ) = plt.subplots( 2, 1, figsize=( 7, 5 ) )
	
	for idx, cultivar in enumerate( genotype_order ):
		ax.boxplot( data_to_plot[ cultivar ], positions= positions[ cultivar ], widths=0.75, patch_artist=True, #showmeans=True, #notch=True,
							boxprops=dict( facecolor=cultivar_colors[ cultivar ], color=cultivar_colors[ cultivar ]),
							capprops=dict(color=cultivar_colors[ cultivar ]),
							whiskerprops=dict(color=cultivar_colors[ cultivar ]),
							flierprops=dict(color=cultivar_colors[ cultivar ], markeredgecolor=cultivar_colors[ cultivar ]),
							medianprops=dict(color="black"),
							meanprops=dict(color="black"),
							zorder=0
							 )
		ax2.boxplot( data_to_plot[ cultivar ], positions= positions[ cultivar ], widths=0.75, patch_artist=True, #showmeans=True, #notch=True,
							boxprops=dict( facecolor=cultivar_colors[ cultivar ], color=cultivar_colors[ cultivar ]),
							capprops=dict(color=cultivar_colors[ cultivar ]),
							whiskerprops=dict(color=cultivar_colors[ cultivar ]),
							flierprops=dict(color=cultivar_colors[ cultivar ], markeredgecolor=cultivar_colors[ cultivar ]),
							medianprops=dict(color="black"),
							meanprops=dict(color="black"),
							zorder=0
							 )
	
	ax2.xaxis.set_ticks( label_positions )
	ax2.set_xticklabels( gene_lables, fontsize=10, rotation=90 )
	
	#ax.set_ylabel( "transcript abundance [RPKMs]", fontsize=10 )
	ax2.text( -5, upper_border, "transcript abundance [RPKMs]", fontsize=10, rotation=90, ha="center", va="center" )
	ax.set_xlim( -0.5, max( [ x for sublist in positions.values() for x in sublist ] )+1 )
	ax2.set_xlim( -0.5, max( [ x for sublist in positions.values() for x in sublist ] )+1 )
	
	simplified_list =  [ x for sublist in data_to_plot.values() for x in sublist ]
	ax.set_ylim( lower_border+1, 1.01*max( [ x for sublist in simplified_list for x in sublist ] ) )
	ax2.set_ylim( 0, upper_border )
	
	ax.tick_params(axis='y', which='major', labelsize=10)
	ax.tick_params(axis='y', which='minor', labelsize=10)
	ax2.tick_params(axis='y', which='major', labelsize=10)
	ax2.tick_params(axis='y', which='minor', labelsize=10)
	
	ax.spines['bottom'].set_visible(False)

	ax2.spines['top'].set_visible(False)
	ax.set_xticks( [] )
	ax.tick_params(labeltop='off')
	ax2.xaxis.tick_bottom()
	
	d = .005  # how big to make the diagonal lines in axes coordinates
	# arguments to pass to plot, just so we don't keep repeating them
	kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
	ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
	ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

	kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
	ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
	ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
	
	my_legend = [	mpatches.Patch(color=cultivar_colors[ "DH" ], label='$\it{Hylocereus}$ $\it{polyrhizus}$ Da Hong (DH, red)'),
								mpatches.Patch(color=cultivar_colors[ "FR" ], label='$\it{Hylocereus}$ $\it{polyrhizus}$ x $\it{undatus}$ Fen Rou (FR, pink)'),
								mpatches.Patch(color=cultivar_colors[ "BR" ], label='$\it{Hylocereus}$ $\it{undatus}$ Bai Rou (BR, white)')
								#mpatches.Patch(color=cultivar_colors[ "BSJ" ], label='BSJ')
							]
	#ax.legend( handles=my_legend, loc="upper left", ncol=1, bbox_to_anchor=(0.425, 0.95), fontsize=8.5 )
	ax.legend( handles=my_legend, loc="upper right", ncol=1, bbox_to_anchor=(0.999, 0.95), fontsize=8.5 )
	
	plt.subplots_adjust( left=0.1, right=0.999, top=0.99, bottom=0.13, hspace=0.03  )
	
	fig.savefig( figfile, dpi=300 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief run generation of plots """
	
	gene_file = arguments[ arguments.index('--genes')+1 ]
	exp_file = arguments[ arguments.index('--exp')+1 ]
	figfile = arguments[ arguments.index('--out')+1 ]

	sample_groups = { 'BR': [ "SRR11190794", "SRR11190793", "SRR11190792" ],
									'DH': [ "SRR11190802", "SRR11190801", "SRR11190798" ],
									'FR': [ "SRR11190797", "SRR11190796", "SRR11190795" ],
									'BSJ': [ "SRR11190791", "SRR11190800", "SRR11190799" ]
								}
	#DH = red; FR=pink, BR=white; BSJ=white
	genotype_order = [ "DH", "FR", "BR" ]	#, "BSJ"

	genes, gene_order = load_genes( gene_file )

	exp = load_exp( exp_file )

	generate_gene_exp_figure( figfile, genes, exp, gene_order, sample_groups, genotype_order, 180, 280 )


if '--genes' in sys.argv and '--exp' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
