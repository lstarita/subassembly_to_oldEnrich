##base trimmer.	 Trims reads down to barcodes and/or useful sequence

from Bio.Seq import Seq
from optparse import OptionParser
import sys
import re
import os


	
def Trim_Filter_RC( reads, barcodesout, bp, runlen, commonS, bp1, taglowqual, scaling, rc, debug, name ):
	#regular expression to find Ns and runs of homopolymers
	re1 = re.compile( "[Nn]" )
	re2 = re.compile( "[Aa]{%d,}|[Gg]{%d,}|[Cc]{%d,}|[Tt]{%d,}" % ( runlen, runlen, runlen, runlen ), flags=re.IGNORECASE )
	char = chr(scaling+2)
	discarded = [ 0 for _ in xrange(5) ]
	readctr = 0
	
	
	
	for line in reads:
		readctr+=1
		header = line.strip()
		read = reads.next().strip()
		barcode = read[0:bp]				
		reads.next()
		qualscore = reads.next().strip()[0:bp]
		i = [header, barcode, qualscore]
	
		localdiscard = [ "0" for _ in xrange(3) ]
		shouldprint = 1
		if re1.search( i[1] ):
			shouldprint = 0
			discarded[0]+=1
		if re2.search( i[1] ):	  
			shouldprint = 0
			discarded[1]+=1
			localdiscard[0]="1"
		if screenStringForChar( i[2], char ) > taglowqual:
			shouldprint = 0
			discarded[2]+=1
			localdiscard[1]="1"
		
		if hammingDistance(read[bp1: (bp1 + len(commonS))], commonS) > 3:
			shouldprint = 0
			discarded[3] += 1
			localdiscard[2]="1"
		
		if debug:
			sys.stderr.write( "\t".join( localdiscard ) + "\n" )
			sys.stderr.flush()
		
		if shouldprint:
			if rc == "y":
				barcode_seq = Seq(barcode)
				barcode_rc = barcode_seq.reverse_complement()
				barcode_rc_str = barcode_rc.tostring()
				barcodesout.write( "\n".join([header, barcode_rc_str, "+", qualscore[::-1]]) + "\n")
			else:
				barcodesout.write( "\n".join([header, barcode, "+", qualscore]) + "\n")
		else:
			discarded[4]+=1
		
		ctrStatusPrinter( readctr, 1000000, "Processed %d reads...\n", name)
		
	strings = ( "N in tag:", "homopolymer tag:", "low quality tag:", "doesn't match common seq:", "total actually discarded:" )
	for i,v in enumerate( discarded ):
		sys.stdout.write( "%s\t%d\n" % ( strings[i], v ) )
		sys.stdout.flush()
	
def screenStringForChar( string, char ):
	ctr=0
	for x in string:
		if x==char:
			ctr+=1
	return ctr

def ctrStatusPrinter( ctr, interval, forprint, name ):
	if ctr % interval == 0:
		sys.stderr.write( name + '\t' + forprint % ctr )
		sys.stderr.flush()

def hammingDistance(seq1, seq2):
	return len(seq1) - sum([ seq1[x] == seq2[x] for x in range(len(seq1))])



if __name__ == '__main__':	
	parser = OptionParser()
	parser.add_option('--rc',  action = 'store', type = 'string', dest = 'rc', help = "reverse complement, y or n")
	parser.add_option('--bp', action = 'store', type = 'int', dest = 'bp', help = "number of basepairs to trim to")
	##parser.add_option('-b', '--barcodes', action = 'store', type = 'string', dest = 'bc_reads', help = "name of barcode fq file")
	parser.add_option('-p', '--path', action = 'store', type = 'string', dest = 'path', help = "path to fq barcode files")
	##parser.add_option('-o', '--outfile', action = 'store', type = 'string', dest = 'out', help = "name for out fq file")
	parser.add_option("--runlen",
						  action="store", type='int', dest="runlen",
						  help="minimum homopolymer length to discard a tag, recommend 7 ( default = %default )", default=7)
	parser.add_option("--commonSeq",
						  action="store", type='string', dest="commonSeq",
						  help="common sequence in vector to match")
	parser.add_option("--bp1",
						  action="store", type='int', dest="bp1",
						  help="index of 1st base of common seq") 

	parser.add_option("--taglowqual",
						  action="store", type='int', dest="taglowqual",
						  help="maximum number of Q2 bases in tag to allow, recommend 4 ( default = %default )", default=4)
	
	parser.add_option("--scaling",
						  action="store", type='int', dest="scaling",
						  help="phred scaling ( default = %default )", default=33)
	parser.add_option("--debug",
					  action="store_true", dest="debug",
					  help="print whether each read passes various filters ( default = %default )", default=False)
	(option, args) = parser.parse_args()
						  
	dirlist = os.listdir(option.path)
	print dirlist
	
	if 'barcodeTrimFilterRC_output' not in dirlist:
		os.mkdir(option.path + "barcodeTrimFilterRC_output")
				

	for item in dirlist:
		if '.fq' in item:
			name = item.strip().split(".")[0]
			print name
			with open( option.path + item ) as I:
				with open( option.path + "barcodeTrimFilterRC_output/" + name + "_barcodeTrimFilter.fq", 'w+b' ) as Iout:
					Trim_Filter_RC( I, Iout, option.bp, option.runlen, option.commonSeq, option.bp1, option.taglowqual, option.scaling, option.rc, option.debug, name )
					
						