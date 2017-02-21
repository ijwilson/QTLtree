#!/usr/bin/python

""" runs sima, collects the results and converts to a useful format"""

import sys
import subprocess
import random
import bisect
import string

def read_tnt(inf):
    """ read tnt matrix output from a stream"""
    line = inf.readline()
    b=line.strip().split()
    if len(b)==2:
        rows,cols=int(b[0]),int(b[1])
        m=list()
        for r in range(rows):
            m.append([int(x) for x in inf.readline().strip().split()])
        return m
        
    elif len(b)==1:
        res = [int(x) for x in inf.readline().strip().split()]
        if len(res)!=int(b[0]):
            print("Error, line is of wrong length in read_tnt, read",line)
            sys.exit(-1)
        return res        
    else:
        print("Error, should read a line with one or two integers in  read_tnt_matrix, read",line)
        sys.exit(-1)


def readsima(in_file,pops=True):
    """ parses the results of a sima run"""
    b = in_file.readline()
    m = read_tnt(in_file)
    pos = read_tnt(in_file)
    if pops: pops=read_tnt(in_file)
    else: pops=[-1]
    return [m,pos,pops]
    

def sima(**kw):
    """ read sima with the list of options given by command line arguments
    and then read the results into a list which is returned"""
    args = ["sima"] + ["--%s=%s" % (k,v) for k,v in kw.iteritems()]
    pr = subprocess.Popen(args,stdout=subprocess.PIPE)
    b=readsima(pr.stdout)
    return(b)


def WeightedChoice(n,p):
    """ get a sample of n integers from 0 to len(p)-1 weighted by
    the values in p (that need not add to 1"""
    choices=range(len(p))
    cp = list(accumulate(p))
    return [choices[bisect.bisect(cp,random.random()*cp[-1])] for i in range(n)]
    

def accumulate(xx):
    """ get a vector with  cumulative sums of the first i values of xx"""
    yy=xx
    for i in range(1,len(xx)): yy[i] = yy[i-1]+xx[i]
    return yy


def writeGenotypes(m,filestem="genotypes"):
    
    SNPtypes = ["CT","TC","AG","GA","AC","CA","CG","GC","AT","TA","GT","TG"]
    SNPfreqtype = [0.23,0.23,0.23,0.23,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]
    
    SNPvar = []
   
    for chr,xx in enumerate(m):
        nSNP=len(xx[1]) 
        out = file(filestem+str(chr+1)+".csv","w")
        out.write("ID,"+",".join("SNP%s%d" % (string.ascii_lowercase[chr],i+1) for i in range(nSNP))+"\n")
        SNP = WeightedChoice(nSNP,SNPfreqtype)
        SNPvar.append([SNPtypes[i] for i in SNP])
    
        for i in range(0,len(xx[0]),2):
            out.write("sample%d" % ((i+2)/2,))
            for j in range(nSNP):
                out.write(","+SNPtypes[SNP[j]][xx[0][i][j]] + "/" +SNPtypes[SNP[j]][xx[0][i+1][j]])
            out.write("\n")
    
        out.close()
    return SNPvar


def writeQTL(q,filename="qtl.csv"):
    out = file(filename,"w")
    out.write("ID,QTL1\n")
    for i in range(0,len(q),2):     out.write("sample%d,%g\n" % ((i+2)/2,q[i]+q[i+1]))
    out.close()


def writeMarkers(xx,SNPs,filename="markers.csv"):
    """ takes a list of lists xx [[pos11,pos12,pos13,...,pos1k],[pos21,pos22,pos32...] ...]
    and outputs a list with snp_name marker position gene snp_type variant and maf"""
    out = file(filename,"w")
    out.write("snp_name,chromosome,position,gene,snp_type,variant,maf\n")

    for chr,yy in enumerate(xx):
        SNPstem = "SNP"+string.ascii_lowercase[chr]
        for index,pos in enumerate(yy[1]):
            out.write("%s%d,%d,%d,None,None,%s,0.2\n" % (SNPstem,index+1,chr+1,pos,SNPs[chr][index][1]))
    
    out.close()
    
    
def writeTargets(chromosomes,chrlength,tarlength=1000,overlap=800,filename="targets.csv"):
    out = open(filename,"w")
    out.write("gene,chromosome,position,start,end,length,orientation,geneid\n")
    count=1
    for i in range(chromosomes):
        for j in range(1,chrlength-1,tarlength-overlap):
            out.write("gene"+str(i+1)+"_"+str(count)+","+str(i+1)+",%d,%d,%d,%d,+,%d\n" % (j+tarlength/2,j,j+tarlength-1,tarlength,count)) 
            count+=1
            
    out.close()
            
def writeQTLregions(reg,filename="targets.csv"):
    out = open(filename,"w")
    out.write("gene,chromosome,start,end\n")
    for i,yy in enumerate(reg):
        out.write("gene"+str(i+1)+","+str(yy[0]+1)+","+str(yy[1])+","+str(yy[2])+","+str(yy[3])+"\n")
    out.close()
    
    
def writePED(xx,filename="pops.ped"):
    """ writes a ped file in the same style as the GAW 17 data"""
    out = file(filename,"w")
    out.write("ID,SEX,AGE,Population\n")
    for i in range(0,len(xx),2):
        out.write("sample%d,%d,%d,pop%d\n" % ((i+2)/2,random.randint(1,2),random.randint(40,80),xx[i]+1))
    out.close()
    

def main():
    chromosomes=22    ## Number of replicate chromosomes
    sites=50000       ## Sites per chromosome
    theta=100         ## Value of theta
    migmodel='Island(4,4)'
    growthmodel="exponential(10)"
    sampleSize=200
    QTLregionsize=1000 ## length in baces of each QTL region
    nQTL=5             ## number of QTL regions
    QTLincrease=0.5    ## part of a standard deviation each mutation in a 
                       ## QTL gene increases the QTL by
    
    data=list()

    for i in range(chromosomes):
        seed = random.randint(1,10000)
        sys.stderr.write("chromosome "+str(i+1)+"\n")
        a = sima(ss=sampleSize,sites=sites,theta=theta,migmatrix=migmodel,growthmodel=growthmodel,seed=seed)
        data.append(a)
        
    targetm=[[] for i in range(len(data[0][0]))]
    QTLregions=list()
    ## Now get the QTLregions
    for ii in range(nQTL):
        chr=random.randrange(chromosomes)
        regionStart=random.randint(QTLregionsize,sites-2*QTLregionsize)
        QTLregions.append([chr,regionStart,regionStart+QTLregionsize])
       
    for region in QTLregions:
        dat=data[region[0]]
        useIndex = [index for index,value in enumerate(dat[1]) if value >= region[1] and value < region[2]]  ## which rows to use
        region.append(len(useIndex))
        for i,row in enumerate(dat[0]):
            targetm[i].extend(row[j] for j in useIndex)
           
    indval=[sum(xx) for xx in targetm]
    QTLval = [random.gauss(0,1)+QTLincrease*y for y in indval]
    
    SNPs = writeGenotypes(data,"genotypes")
    writePED(data[0][2],"pops.ped");                          ## should all be the same
    writeMarkers(data,SNPs,"markers.csv")
    writeQTL(QTLval,"qtl.csv")    
    writeQTL(indval,"muts.csv")

    writeQTLregions(QTLregions,"QTLregions.csv")

    writeTargets(chromosomes,sites,1000,800,"targets.csv")


if __name__=="__main__": main()
