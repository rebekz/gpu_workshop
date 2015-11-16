#this is parallel version of gwas written with pyCuda

from sys import argv
import numpy
import scipy.stats
import time
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

start = time.clock()

txt = open(argv[1], 'r')
txt2 = open(argv[2], 'r')
y = txt2.read()
filepath = argv[1].split("/",1)
filepath = filepath[0]
#calculate number of sample

Y = numpy.array(list(y)).astype(int) - 1
num_sample = len(Y)

# find the position of control
i = 0
while Y[i] != 0:
    i = i + 1

idx_ctrl = i


#########################################
num_f = 0

list_files = txt.readline()

while list_files:
    files = list_files.replace("\n","")
    files = "%s/%s" % (filepath, files)
    print files
    txt_geno = open(files, 'r')
    #create array of geno
    geno = txt_geno.read()
    genos = numpy.array(geno.split("\n"))

    gt = []
    gt_name = []
    for i in genos:
        isi =  i.split(' ', 1)
        gt_name.append(isi[0])
        gt.append(isi[1].replace(" ",""))

    num_threads = len(gt)
    flag = False
    num_threads_2 = 0
    st = 0
    if(num_threads > 1000):
        num_blocks = round(num_threads / 1000, 1)
        num_t = num_threads % 1000
        if(num_t != 0):
           flag = True
           num_threads_2 = num_t
           st = num_blocks * 1000
        num_threads = 1000
    else:
        num_blocks = 1

    g = numpy.array(gt)

    load_up_file = time.clock()

    mod = SourceModule("""

    #include <stdio.h>
    #include <string.h>
    #include <math.h>

    __global__ void map(char *a, float *chi2s, int numCols, int ctrl)
    {

     int threadsPerblock = blockDim.x * blockDim.y;
     int blockId = blockIdx.x + (blockIdx.y * gridDim.x);
     int threadId = threadIdx.x + (threadIdx.y * blockDim.x);
     int idx = (blockId * threadsPerblock) + threadId; 
     //int idx = blockIdx.x * blockDim.x + threadIdx.x;

     int start = idx * numCols;
     int stop = start + numCols-1;
     int ctrl_start = idx * numCols + ctrl;
     int zero_case = 0;
     int zero_ctrl = 0;
     int one_case = 0;
     int one_ctrl = 0;
     int two_case = 0;
     int two_ctrl = 0;
     float N = numCols;
     float R = ctrl;
     float r1 = 0;
     float r2 = 0;
     float n1 = 0;
     float n2 = 0;
     float chi1 = 0;
     float chi2 = 0;

     //printf("start: %d ",start);

     for(int i = start; i <= stop; i++) 
     {
      int temp = a[i] - '0';
      //printf("%d == %d ", i, ctrl_start); 
      if(i == ctrl_start) {
       ctrl_start = ctrl_start + 1;
        switch(temp)
        {
         case 0:
          zero_ctrl = zero_ctrl + 1;
          break;
         case 1:
          one_ctrl = one_ctrl + 1;
          break;
         case 2:
          two_ctrl = two_ctrl + 1;
          break;
        }
      } else {
        switch(temp)
        {
         case 0:
          zero_case = zero_case + 1;
          break;
         case 1:
          one_case = one_case + 1;
          break;
         case 2:
          two_case = two_case + 1;
          break;
        }
      }
     }

     
     r1 = one_case;
     r2 = two_case;
     n1 = one_ctrl + one_case;
     n2 = two_ctrl + two_case;
     
     chi1 = N*(r1+2*r2)-R*(n1+2*n2);
     chi1 = powf(chi1,2);
     chi1 = N*chi1;
     
     float chi2_left = (N-R)*R;
     float chi2_right = N*(n1+4*n2);
     float chi2_right_2 = (n1+2*n2);
     chi2_right_2 = powf(chi2_right_2,2);
     chi2_right = chi2_right - chi2_right_2;
     chi2 = chi2_left * chi2_right;

     chi2 = chi1/chi2;
     chi2s[idx] = chi2;
     
    }

    __global__ void map_2(char *a, float *chi2s, int numCols, int ctrl, int st)
    {

     int threadsPerblock = blockDim.x * blockDim.y;
     int blockId = blockIdx.x + (blockIdx.y * gridDim.x);
     int threadId = threadIdx.x + (threadIdx.y * blockDim.x);
     int idx = (blockId * threadsPerblock) + threadId; 
     //int idx = blockIdx.x * blockDim.x + threadIdx.x;

     int start = (numCols * idx) + st * numCols;
     int stop = start + numCols-1;
     int ctrl_start = (numCols * idx) + st * numCols + ctrl;
     int zero_case = 0;
     int zero_ctrl = 0;
     int one_case = 0;
     int one_ctrl = 0;
     int two_case = 0;
     int two_ctrl = 0;
     float N = numCols;
     float R = ctrl;
     float r1 = 0;
     float r2 = 0;
     float n1 = 0;
     float n2 = 0;
     float chi1 = 0;
     float chi2 = 0;

     //printf("start %d ", start);
     for(int i = start; i <= stop; i++) 
     {
      int temp = a[i] - '0';
      if(i == ctrl_start) {
       ctrl_start = ctrl_start + 1;
        switch(temp)
        {
         case 0:
          zero_ctrl = zero_ctrl + 1;
          break;
         case 1:
          one_ctrl = one_ctrl + 1;
          break;
         case 2:
          two_ctrl = two_ctrl + 1;
          break;
        }
      } else {
        switch(temp)
        {
         case 0:
          zero_case = zero_case + 1;
          break;
         case 1:
          one_case = one_case + 1;
          break;
         case 2:
          two_case = two_case + 1;
          break;
        }
      }
     }

     
     r1 = one_case;
     r2 = two_case;
     n1 = one_ctrl + one_case;
     n2 = two_ctrl + two_case;
     
     chi1 = N*(r1+2*r2)-R*(n1+2*n2);
     chi1 = powf(chi1,2);
     chi1 = N*chi1;
     
     float chi2_left = (N-R)*R;
     float chi2_right = N*(n1+4*n2);
     float chi2_right_2 = (n1+2*n2);
     chi2_right_2 = powf(chi2_right_2,2);
     chi2_right = chi2_right - chi2_right_2;
     chi2 = chi2_left * chi2_right;

     chi2 = chi1/chi2;
     int new_idx = st + idx;
     chi2s[new_idx] = chi2;
     
    }


    """)

    map = mod.get_function("map")
    ns_gpu = numpy.int32(num_sample)
    idx_ctrl_gpu = numpy.int32(idx_ctrl)
    chi2s = numpy.zeros(g.size).astype(numpy.float32)

    map(cuda.In(g), cuda.Out(chi2s), ns_gpu, idx_ctrl_gpu, block=(num_threads,1,1), grid=(int(num_blocks),1))
    cuda.Context.synchronize() 
    if(flag == True):
       st_gpu = numpy.int32(st)
       map_2 = mod.get_function("map_2")
       map_2(cuda.In(g), cuda.Out(chi2s), ns_gpu, idx_ctrl_gpu, st_gpu, block=(num_threads_2,1,1), grid=(1,1))
       cuda.Context.synchronize()
    gpu_execution = time.clock()
    f = open('output/out_par_%s.txt' % num_f,'w')

    for idx, val in enumerate(chi2s.tolist()):
        pval = 1 - scipy.stats.chi2.cdf(val,1)
        #print gt_name[idx], " ", pval
        wr = "%s %s \n" % (gt_name[idx], pval)
        f.write(wr)
    f.close()
    txt_geno.close()
    num_f = num_f + 1
    list_files = txt.readline()

end = time.clock()

print "Time execution : ", end - start, " s"
