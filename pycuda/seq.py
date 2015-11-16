#this is sequential version of gwas

from sys import argv
import numpy as np
import scipy.stats
import pandas
import time

start = time.clock()

#open list of files
txt = open(argv[1])
filepath = argv[1].split("/")
filepath = filepath[0]
#open and read disease status Y
y = open(argv[2])
y_line = y.read()
Y = np.array(list(y_line)).astype(int) - 1

num_f = 0

list_files = txt.readline()

while list_files:
    files = list_files.replace("\n","")
    files = "%s/%s" % (filepath,files)
    txt_geno = open(files, 'r')
    print files
    #write to txt
    f = open('output/out_seq_%s.txt' % num_f,'w')

    line = txt_geno.readline()
    #while there are more rows
    while line:
        #split key and value
        isi = line.split(' ', 1)
        key = isi[0]
        value = isi[1]

        #parse each SNP
        gt = np.array(value.split(' ')).astype(int)   
        #mapping disease status with each SNP
        t = zip(Y,gt)
        #get count frequency array
        ps = pandas.Series([tuple(i) for i in t])
        counts = ps.value_counts()

        try:
            one_zero = counts[(1,0)]
        except KeyError:
            one_zero = 0
        try:
            one_one = counts[(1,1)]
        except KeyError:
            one_one = 0
        try:
            one_two = counts[(1,2)]
        except KeyError:
            one_two = 0
        try:
            zero_one = counts[(0,1)]
        except KeyError:
            zero_one = 0
        try:
            zero_zero = counts[(0,0)]
        except KeyError:
            zero_zero = 0
        try:
            zero_two = counts[(0,2)]
        except KeyError:
            zero_two = 0

        N = len(t)
        R = one_zero + one_one + one_two
        r1 = one_one
        r2 = one_two
        n1 = zero_one + one_one
        n2 = zero_two + one_two
        
        #N = len(t) #sample size
        #R = counts[(1,0)] + counts[(1,1)] + counts[(1,2)] #cases
        #r1 = counts[(1,1)] #cases with one copy of the disease allele
        #r2 = counts[(1,2)] #cases with two copies of the disease allele
        #n1 = counts[(0,1)] + counts[(1,1)] #number samples with one copy of the disease allele
        #n2 = counts[(0,2)] + counts[(1,2)] #number samples with two copies of the disease allele
        
        chi1 = float(N*(N*(r1+2*r2)-R*(n1+2*n2)) ** 2)    
        try:
            chi2 = chi1/((N-R)*R*(N*(n1+4*n2)-(n1+2*n2) ** 2))
        except ZeroDivisionError:
            chi2 = 0

        pval = 1 - scipy.stats.chi2.cdf(chi2,1)

        #print N, R, r1, r2, n1, n2, chi1, chi2, pval
        #print counts
        #print key, " ", pval
        wr = "%s %s \n" % (key, pval)
        f.write(wr)

        line = txt_geno.readline()
    f.close()
    txt_geno.close()
    num_f = num_f + 1
    list_files = txt.readline()
    

end = time.clock()
print "Time execution: ", end - start, " s"
