##script to generate false enrich 0.2 read aligner output fiel
## takes a directory of fq file of trimmed and filtered barcodes, DNA and PRO read aligner files from subassembly as input
## converts the fq file into a false read_aligner output file. 

import sys, os
from optparse import OptionParser

def makeAlignedDict(d, p, DNAdict, PROdict):
	for line in d:
		if line.startswith('readID') : continue
		
		else:
			readID, sequence, match_count, mutation_count, mutation_location, mutation_identity, max_mutation_run = line.strip().split()
			barcode = readID[1:]
			
			
			DNAdict[barcode] = line.strip().split()
	
	for line in p:
		if line.startswith('readID') : continue
		else:
			readID, sequence, match_count, mutation_count, mutation_location, mutation_identity, max_mutation_run = line.strip().split()
			barcode = readID[1:]
			
			PROdict[barcode] = line.strip().split()
	return (DNAdict, PROdict)

def buildAlignFiles(DNAdict, PROdict, barcode_file, DNAout, PROout, stats, name):
	ctr =	0
	noAssembly = []
	match = []
	
	for line in barcode_file:
		ctr += 1
		barcode = barcode_file.next().strip()
		barcode_file.next().strip()
		barcode_file.next().strip()
		
		if barcode in DNAdict:
			DNAout.write("\t".join( DNAdict[barcode]) + '\n')
			PROout.write("\t".join( PROdict[barcode]) + '\n')
			match.append(barcode)
		else:
			noAssembly.append(barcode) 
			continue
		
		ctrStatusPrinter( ctr, 1000000, "Hashed %i barcodes...\n", name )
		
	stats.write(name + '\t' 
				 + str( ctr ) + '\t' 
				 + str(len(noAssembly)) + '\t'
				 + str(len(set(noAssembly))) + '\t' 
				 + str(len(match)) + '\t'
				 + str(len(set(match))) +'\n')
	

def ctrStatusPrinter( ctr, interval, forprint, name ):
	if ctr % interval == 0:
		sys.stderr.write(	 name + '\t'	+ forprint % ctr )
		sys.stderr.flush()
		
if __name__ == '__main__':

	parser = OptionParser()
	parser.add_option('--DNA', action = 'store', type = 'string', dest = 'dna', help = "name of DNA subAssReadAlignerFile")
	parser.add_option('--PRO', action = 'store', type = 'string', dest = 'pro', help = "name of PRO subAssReadAlignerFile")
	parser.add_option('-p','--paths_to_barcodes', action = 'store', type = 'string', dest = 'path', help = "directory of filtered trimmed barcode fq files")
	(option, args) = parser.parse_args()


	DNAdict = {}
	PROdict = {}	
	
	dirlist = os.listdir(option.path)
	print dirlist
	
	if 'barcode_linked_aligned' not in dirlist:
		os.mkdir(option.path + "barcode_linked_aligned" )
				
	with open(option.path +  "barcode_linked_aligned/stats.txt", 'w') as stats:
		stats.write("sample" + '\t' 
				 + "total barcodes" + '\t' 
				 + "total_NotFoundInAssembly" + '\t'
				 + "unique_NotFoundInAssembly" + '\t' 
				 + "total_AssemblyMatch" + '\t'
				 + "unique_AssemblyMatch" +'\n')	
				 

		with open( option.dna, 'r' ) as DNAin:
			DNAheader = DNAin.readline()
			with open( option.pro, 'r' ) as PROin:	
				PROheader = PROin.readline()
				
				(DNAdict, PROdict) = makeAlignedDict(DNAin, PROin, DNAdict, PROdict)

		for item in dirlist:
			if '.fq' in item:
				name = '_'.join(item.strip().split("_")[0:2])
				print name
				with open( option.path + item ) as barcodef:
					with open( option.path + "barcode_linked_aligned/" + name + "_R1_DNA_qc", 'w+b' ) as DNAout:
						DNAout.write(DNAheader )
						with open( option.path + "barcode_linked_aligned/" + name + "_R1_PRO_qc", 'w+b' ) as PROout:
							PROout.write(PROheader )
							buildAlignFiles(DNAdict, PROdict, barcodef, DNAout, PROout, stats, name)


