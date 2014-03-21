 ##script to reconnect barcode with sequence after subassembly.  Will to link barcode to sequence to run Nrich
 
 ##need .counts text file that was used by readNamer.py and subassembly out file
 
 ## .counts file should be sorted and the rank number corresponds to the 1st number of the SA out file
 
 
import sys 
from optparse import OptionParser
import cPickle

def readAndHashBarcodes( barcodes, dict ):
    with open( barcodes ) as f:
        ctr = 1
        for line in f:
            line = line.rstrip( '\n' )
            linesplit = line.split ( '\t' )
            dict[ ctr ] = ( linesplit[ 0 ], int( linesplit[ 1 ] ) )
            ctr += 1
            ctrStatusPrinter( ctr, 10000, "Hashed %i barcodes...\n", sys.stderr )
            
def ctrStatusPrinter( ctr, interval, forprint, fileh ):
    if ctr % interval == 0:
        sys.stderr.write( forprint % ctr )
        sys.stderr.flush()


def fakeFastq(SA_out_file, barcodeDict, ofile, assemblyDict):
    with open(SA_out_file) as sa_out:
        with open(ofile, 'w') as writefile:
            ctr = 0
            for line in sa_out:
                ctr += 1
                items = line.strip().split()
                bc_id = int(items[0].split('_')[1])
                seq = items[1]
                qualscore = items[2]
            
        
                barcode = barcodeDict[bc_id][0]
                assert  int(items[0].split('_')[2]) == barcodeDict[bc_id][1]
                assemblyDict[barcode] = (seq, qualscore)
                ctrStatusPrinter( ctr, 10000, "linked %i barcodes...\n", sys.stderr )
                writefile.write( "@" + barcode + '\n'  + seq + "\n" + "+" + "\n" + qualscore + "\n")
                

if __name__ == '__main__':

    parser = OptionParser()
    parser.add_option('-o', '--out_file', action = 'store', type = 'string', dest = 'out', help = "name of out file")
    parser.add_option('--counts', action = 'store', type = 'string', dest = 'counts', help = "barcode counts file used for naming reads")
    parser.add_option('-f','--SA_file', action = 'store', type = 'string', dest = 'SA_out', help = "name of SA out file")
    (option, args) = parser.parse_args()


    bdict = {}
    readAndHashBarcodes( option.counts, bdict )

    adict = {}
    fakeFastq(option.SA_out, bdict, option.out, adict)
    

        
