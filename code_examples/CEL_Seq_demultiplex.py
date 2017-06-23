#!/share/apps/6.0/python/Python-2.7.3_gnu446_x86_64/bin/python2.7

## Josh Campbell
## 9/10/2014
## This script takes in a SampleSheet.csv file, parses it, then parses
## the fastq files in order to split them up into the proper samples

## Example of the comma separated file with header:

# Sample_ID,Plate_Barcode,Cell_Barcode,Lane,FASTQ_R1_Path,FASTQ_R2_Path
# HK293_1,AGTC,CATCACGC,1,/path/to/fastq1,/path/to/fastq2
# H1299_1,AGTC,GTCGTCGC,1,/path/to/fastq1,/path/to/fastq2
# HK293_2,AGTC,CATCAATC,1,/path/to/fastq1,/path/to/fastq2

import numpy as np
import time
import argparse
import gzip
import distance

## Functions

## Determines the closest matching 3' and 5' barcodes within the max mismatch limit
## If more than 1 barcode matches with <= max mismatch, then no barcodes are chosen
def best_match(my_bc5, my_bc3, bc5_list, bc3_list, bc5_max_mismatch, bc3_max_mismatch):
  
  ## Find the edit distance between the given 3' barcode and all barcodes within the list
  edit_dist3 = [];
  for i in bc3_list:
    edit_dist3.append(distance.hamming(my_bc3, i))
  edit_dist3_min = min(edit_dist3)
  
  ## Find the edit distance between the given 5' barcode and all barcodes within the list
  edit_dist5 = [];
  for i in bc5_list:
    edit_dist5.append(distance.hamming(my_bc5, i))
  edit_dist5_min = min(edit_dist5)
  
  ## If the number of mismatches is lower than the given threshold AND only one barcode has the lowest
  ## number of mismatches, then the index is assigned
  if edit_dist3_min <= bc3_max_mismatch and edit_dist3.count(edit_dist3_min) == 1 and edit_dist5_min <= bc5_max_mismatch and edit_dist5.count(edit_dist5_min) == 1:
    tmp_bc3 = bc3_list[edit_dist3.index(edit_dist3_min)]
    tmp_bc5 = bc5_list[edit_dist5.index(edit_dist5_min)]
  else: 
    tmp_bc3 = ""
    tmp_bc5 = ""
  
  return (tmp_bc5, tmp_bc3)

################################
## Main
################################

## Import options from command line
parser = argparse.ArgumentParser(description='Demultiplexes MARS-Seq FASTQs.')
parser.add_argument('sample_sheet_name', metavar='SampleSheet.csv', type=str, nargs=1, help='Comma seperated with Sample Name and Barcodes')
parser.add_argument('-m', "--max_mismatchBC", action='store', default=2, type=int,
                   help='The maximum number of mismatches to allow in the 3-prime barcode [Default=2]')
parser.add_argument('-startBC', "--start_bc", action='store', default=6, type=int,
                   help='The start position of the cell barcode in read R1 [Default=1]')                   
parser.add_argument('-endBC', "--end_bc", action='store', default=11, type=int,
                   help='The end position of the cell barcode in read R1 [Default=1]')                   
parser.add_argument('-startUMI', "--start_UMI", action='store', default=1, type=int,
                   help='The start position of the UMI barcode in read R1 [Default=1]')                   
parser.add_argument('-endUMI', "--end_UMI", action='store', default=5, type=int,
                   help='The end position of the UMI barcode in read R1 [Default=1]')                   

args = parser.parse_args()
sample_sheet_name = args.sample_sheet_name[0]
max_mismatch_5bc = args.max_mismatch_5bc
max_mismatch_3bc = args.max_mismatch_3bc

print "\nSample sheet used:", sample_sheet_name
print "Number of mismatches allowed in cell barcode:", m
print "Cell barcode location in R1:", startBC, "-", endBC
print "UMI barcode location in R1:", startUMI, "-", endUMI

## Read in the sample sheet
sample_sheet = np.genfromtxt(sample_sheet_name, delimiter=',', dtype="string", skip_header=1)


## Parsing sample sheet
fastq_pairs = set()
index_pairs = set()
index_5BC = set()
index_3BC = set()
sample = dict()
sample_check = set()
lanes = set()
for x in range(0, len(sample_sheet)):

  ## Devise a sample name
  sample_name = sample_sheet[x,0] + "_L" + sample_sheet[x,3] + "_5" + sample_sheet[x,1] + "_3" + sample_sheet[x,2]  
  index = (sample_sheet[x,1], sample_sheet[x,2])

  ## Check for unique sample names
  if sample_sheet[x,0] in sample_check:
    raise Exception("Sample ids need to be unique. Offender: %s" % sample_sheet[x,0])
  else:
    sample[index] = sample_name
    sample_check.add(sample_sheet[x,0])
  
  ## Check for unique combination of barcodes
  if index in index_pairs:
    raise Exception("The combination of 5' and 3' indices need to be unique in each row. Offender: %s" % str((sample_sheet[x,1], sample_sheet[x,2])))
  else:
    sample[index] = sample_name
  
  ## Adds relevant info to sets
  fastq_pairs.add((sample_sheet[x,4], sample_sheet[x,5], sample_sheet[x,3], sample_sheet[x,1]))
  index_pairs.add(index)
  index_5BC.add(sample_sheet[x,1])
  index_3BC.add(sample_sheet[x,2])
  lanes.add((sample_sheet[x,3], sample_sheet[x,1]))

## Set up file handles for each sample
fastq1_handles = dict()
fastq2_handles = dict()

for i in index_pairs:
  fastq1_handles[i] = gzip.open(sample[i] + "_R1.fastq.gz", "w")
  fastq2_handles[i] = gzip.open(sample[i] + "_R2.fastq.gz", "w")

for i in lanes:
  fastq1_handles["Unknown_L" + i[0] + "_" i[1]] = gzip.open("Unknown_L" + i[0] + "_" i[1] + "_R1.fastq.gz", "w") 
  fastq2_handles["Unknown_L" + i[0] + "_" i[1]] = gzip.open("Unknown_L" + i[0] + "_" i[1] + "_R2.fastq.gz", "w")


bc3_miss = dict()
bc5_miss = dict()
total_hit = 0
total_miss = 0

## Go through each pair of FASTQ files and demultiplex into proper files
for s in fastq_pairs:
  fastq1 = s[0]
  fastq2 = s[1]
  lane = s[2]
  plateBC = s[3]
  print "\nStarting to demultiplex FASTQ files:", "\n", fastq1, "\n", fastq2, "\n", "on", time.strftime("%c"), "..."

  ## Open a filehandle for the first fastq file
  if fastq1[-3:] == '.gz':
    handle1 = gzip.open(fastq1, "rb")
  else:
    handle1 = open(fastq1, "rU")

  ## Count number of lines in first fastq file
  fastq1_line_count = 0
  for record in handle1:
      fastq1_line_count += 1

  ## Open a filehandle for the second fastq file
  if fastq2[-3:] == '.gz':
    handle2 = gzip.open(fastq2, "rb")
  else:
    handle2 = open(fastq2, "rU")

  ## Count number of lines in second fastq file
  fastq2_line_count = 0
  for record in handle2:
      fastq2_line_count += 1

  handle1.close()
  handle2.close()

  ## Check for equal number of reads between fastq files
  if(fastq1_line_count % 4 != 0):
    raise Exception("The number of lines in %s is not a multiple of 4. Do you have a truncated file?" % fastq1)
  if(fastq2_line_count % 4 != 0):
    raise Exception("The number of lines in %s is not a multiple of 4. Do you have a truncated file?" % fastq2)
  if(fastq1_line_count != fastq2_line_count):
    raise Exception("The number of reads in %(f1)s is %(d1)s while the number of reads in %(f2)s is %(d2)s. Do you have the correct fastq files?" % {"f1": fastq1, "d1": str(fastq1_line_count/4), "f2": fastq2, "d2": str(fastq2_line_count/4)})
  

  print "Number of reads in", fastq1, ":", fastq1_line_count/4
  print "Number of reads in", fastq2, ":", fastq2_line_count/4

  ## Reopen fastq files for splitting
  if fastq1[-3:] == '.gz':
    handle1 = gzip.open(fastq1, "rb")
  else:
    handle1 = open(fastq1, "rU")
  if fastq2[-3:] == '.gz':
    handle2 = gzip.open(fastq2, "rb")
  else:
    handle2 = open(fastq2, "rU")
  
  for id1 in handle1:
    ## Read in next record from both FASTQ files
    id1 = id1.rstrip('\r\n')
    seq1 = next(handle1).rstrip('\r\n')
    qual_id1 = next(handle1).rstrip('\r\n')
    qual1 = next(handle1).rstrip('\r\n')

    id2 = next(handle2).rstrip('\r\n')
    seq2 = next(handle2).rstrip('\r\n')
    qual_id2 = next(handle2).rstrip('\r\n')
    qual2 = next(handle2).rstrip('\r\n')
 
    ## Get Barcodes and UMIs out of reads
    bc3 = seq2[(startBC-1):(endBC-1)]
    umi3 = seq2[(startUMI-1):(endUMI-1)]
    current_index = (bc5, bc3)
    
    ## Check to see if barcodes are an exact match
    ## If not, then the next best match is attempted
    if current_index in index_pairs:
      read_index = current_index
    else:
      read_index = best_match(bc5, bc3, list(index_5BC), list(index_3BC), max_mismatch_5bc, max_mismatch_3bc)

    ## Add Barcde and Unique molecular tag information to sequence and quality headers
    new_id1 = id1 + "_BC5-" + bc5 + ":BC3-" + bc3 + ":UMI3-" + umi3 + "\n"
    new_id2 = id2 + "_BC5-" + bc5 + ":BC3-" + bc3 + ":UMI3-" + umi3 + "\n"
    
    new_id1 = new_id1.replace(" ", "_")
    new_id2 = new_id2.replace(" ", "_")
    
    ## Chop barcodes/UMIs off of the reads and add newline
    seq1 = seq1[7:] + "\n"
    seq2 = seq2[16:] + "\n"
    qual1 = qual1[7:] + "\n"
    qual2 = qual2[16:] + "\n"

    ## If both 3' and 5' barcodes match what is in the SampleSheet
    ## Then the reads are output to the respective file
    if(read_index in index_pairs):
      total_hit += 1

      fastq1_handles[read_index].write(new_id1)
      fastq1_handles[read_index].write(seq1)
      fastq1_handles[read_index].write("+\n")
      fastq1_handles[read_index].write(qual1) 

      fastq2_handles[read_index].write(new_id2)
      fastq2_handles[read_index].write(seq2)
      fastq2_handles[read_index].write("+\n")
      fastq2_handles[read_index].write(qual2)
    else:
      total_miss += 1
      fastq1_handles["Unknown_L" + lane + "_" + plateBC].write(new_id1)
      fastq1_handles["Unknown_L" + lane + "_" + plateBC].write(seq1)
      fastq1_handles["Unknown_L" + lane + "_" + plateBC].write("+\n")
      fastq1_handles["Unknown_L" + lane + "_" + plateBC].write(qual1)

      fastq2_handles["Unknown_L" + lane + "_" + plateBC].write(new_id2)
      fastq2_handles["Unknown_L" + lane + "_" + plateBC].write(seq2)
      fastq2_handles["Unknown_L" + lane + "_" + plateBC].write("+\n")
      fastq2_handles["Unknown_L" + lane + "_" + plateBC].write(qual2)

      ## Count the barcode misses
      if bc3 not in index_3BC:
        if bc3 in bc3_miss:
          bc3_miss[bc3] += 1
        else:
          bc3_miss[bc3] = 1 

  read_count = fastq1_line_count/4
  print "Total number of reads with matched barcode:", total_hit, "(", round(float(total_hit)/float(read_count), 3)*100,  "%)"
  print "Total number of reads without a matched barcode:", total_miss,  "(", round(float(total_miss)/float(read_count), 3)*100,  "%)"
  print "Finished on", time.strftime("%c"), "\n"
  handle1.close()
  handle2.close()

  ## Print out 3' barcode misses for lane
  f = open("Unknown_L" + lane + "_" + plateBC + "_Cell_Barcode_counts.txt", "w")
  for w in sorted(bc3_miss, key=bc3_miss.get, reverse=True):
    f.write(w + "\t" + str(bc3_miss[w]) + "\n")
  f.close()
  
print "\nFinished demultiplexing all files on", time.strftime("%c"), "\n"


for i in fastq1_handles.values():
  i.close()
for i in fastq2_handles.values():
  i.close()



