'''
@author: Marta Luksza, mluksza@ias.edu
'''

import sys
from Aligner import Aligner

from Neoantigen import Neoantigen

def readNeoantigens(neofilename):
    neoantigens={}
    f=open(neofilename)
    header=f.readline()
    htab=header.strip().split("\t")
    hdict={}
    for i in range(0,len(htab)):
        hdict[htab[i]]=i

    line=f.readline()
    while line:   
        line=line.strip()
        nparams=line.split("\t")
        if nparams[7]=="NA":
            line=f.readline()
            continue
        neoantigen=Neoantigen(nparams)
        neoantigens[neoantigen.id]=neoantigen       
        neoantigen.setA()
        line=f.readline()
    f.close()
    samples=set(map(lambda neo: neo.getSampleName(),neoantigens.values()))
    return [neoantigens,samples]

def main(argv):
    
    '''
    command line parameters:
    neofile - text file with neoantigen data (supplementary data)
    alignmentDirectory - folder with precomputed alignments
    a - midpoint parameter of the logistic function, alignment score threshold
    k - slope parameter of the logistic function
    outfile - path to a file where to output neoantigen fitness computation
    '''
        
    neofile=argv[1]
    alignmentDirectory=argv[2]
    a=float(argv[3])
    k=float(argv[4])
    outfile=sys.argv[5]
    
    [neoantigens,samples]=readNeoantigens(neofile)    
    #Compute TCR-recognition probabilities for all neoantigens
    aligner=Aligner()    
    for sample in samples:
        xmlpath=alignmentDirectory+"/neoantigens_"+sample+"_iedb.xml"
        aligner.readAllBlastAlignments(xmlpath)    
    aligner.computeR(a, k)    
    
    #Write neoantigen recognition potential
    of=open(outfile,'w')
    header=["NeoantigenID","Mutation","Sample","MutatedPeptide","ResidueChangeClass","MutantPeptide","WildtypePeptide","A","R","Excluded","NeoantigenRecognitionPotential"]
    header="\t".join(header)
    of.write(header+"\n")
    for i in neoantigens:
        neoantigen=neoantigens[i]
        w=neoantigen.getWeight() #excludes neoantigens that mutated from a nonhydrophobic residue on position 2 or 9
        A=neoantigens[i].getA() #MHC amplitude A
        mtpeptide=neoantigens[i].mtPeptide #mutant peptide
        wtpeptide=neoantigens[i].wtPeptide
        R=aligner.getR(i)        
        
        # Residue change:
        # HH: from hydrophobic to hydrophobic, 
        # NN: from non-hydrophobic to non-hydrophobic
        # HN: from hydrophobic to non-hydrophobic, 
        # NH: from non-hydrophobic to hydrophobic
        # other (WW, WH, HW, NW, WN) which include aminoacids without a clear classification
        residueChange=neoantigen.residueChange 
        
        fitnessCost=A*R*w
        
        l=[i, neoantigen.mid, neoantigen.sample, neoantigen.position, residueChange, mtpeptide, wtpeptide, A,R, 1-w, fitnessCost]#, neoAlignment, epitopeAlignment, score, species]
        l="\t".join(map(lambda s: str(s),l))
        of.write(l+"\n")

if __name__ == '__main__':
    main(sys.argv)
    
