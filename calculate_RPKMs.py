### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python calculate_RPKMs.py
					--counts <COUNT_TABLE>
					--rpkm <RPKM_FILE>
					--assembly <ASSEMBLY_FILE>
					"""

import os, sys

# --- end of imports --- #


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
	return exp, headers


def load_sequences( fasta_file ):
	"""! @brief load candidate gene IDs from file """
	
	sequences = {}
	with open( fasta_file ) as f:
		header = f.readline()[1:].strip()
		if " " in header:
			header = header.split(' ')[0]
		seq = []
		line = f.readline()
		while line:
			if line[0] == '>':
					sequences.update( { header: "".join( seq ) } )
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
					seq = []
			else:
				seq.append( line.strip() )
			line = f.readline()
		sequences.update( { header: "".join( seq ) } )	
	return sequences


def calculate_normed_values( raw_counts, samples, assembly ):
	"""! @brief calculate normalized values """
	
	RPKMs = {}
	for sample in samples:
		n = sum( raw_counts[ sample ].values() )
		sample_size_factor = n / 1000000.0
		RPKMs.update( { sample: {} } )
		genes = raw_counts[ sample ].keys()
		for gene in genes:
			val = raw_counts[ sample ][ gene ]
			length = len( assembly[ gene ] ) / 1000.0
			RPKMs[ sample ].update( { gene: val /  ( length * sample_size_factor ) } )
	return RPKMs, genes


def main( arguments ):
	"""! @brief run everything """
	
	raw_counts_file = arguments[ arguments.index('--counts')+1 ]
	rpkm_file = arguments[ arguments.index('--rpkm')+1 ]
	assembly_file = arguments[ arguments.index('--assembly')+1 ]
	
	assembly = load_sequences( assembly_file )
	
	raw_counts, samples = load_exp( raw_counts_file )
	
	RPKMs, genes = calculate_normed_values( raw_counts, samples, assembly )
	
	with open( rpkm_file, "w" ) as out:
		out.write( "genes\t" + "\t".join( samples ) + "\n" )
		for gene in sorted( genes ):
			new_line = [ gene ]
			for sample in samples:
				new_line.append( RPKMs[ sample ][ gene ] )
			out.write( "\t".join( map( str, new_line ) ) + "\n" )


if '--counts' in sys.argv and '--rpkm' in sys.argv and '--assembly' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
