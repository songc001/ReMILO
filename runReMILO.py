#!/usr/bin/python

import os
import argparse
import shutil
import datetime
global  miserr
parser = argparse.ArgumentParser(description='runReMILO.py')
parser.add_argument('contig_path', metavar='contig.fa', help="The path to contig.fa")
parser.add_argument('reference_path', metavar='reference.fa', help="The path to  reference.fa")
parser.add_argument('short_read1_path',metavar='short_read1.fasta',help="The path to  short_read1_path")
parser.add_argument('short_read2_path',metavar='short_read2.fasta',help="The path to  short_read2_path")
parser.add_argument('-i', "--insert",dest='insert',help="Insert length of shortreads")
parser.add_argument('-k','--kmer', help="Kmer length for LoRDEC refinement.(19)", default=19, type=int)
parser.add_argument('-d','--distance',help="minimum alignment distance  of adjacent  contig  positions  to detect a misassembly error(with either reference  genome or longreads(85)",default=85,type=int)
parser.add_argument('-l', '--longread', help="corresponding  long reads(yes)" )
parser.add_argument('-c','--coverage',help="Minimum number of longreads required to detect  a  misassembly error")
parser.add_argument('-o',"--outdir",action="store",dest="outdir",help="output directory")
args = parser.parse_args()

# Default Parameters#####################################
temp_dir = './temp'
output_dir = './output'
prefix = 'ReMILO'

# Parameters Analyzing###################################
reference_path = args.reference_path
contig_path = args.contig_path
short_read1_path=args.short_read1_path
short_read2_path=args.short_read2_path

if args.insert:
        if args.insert > 128 or args.insert < 1:
        	print 'ERROR: argument -i/insert should be between 1 and 128'
        	exit(-1)
if args.kmer > 127 or args.kmer < 4:
	print 'ERROR: argument -k/--kmer should be between 4 and 127'
	exit(-1)
if args.distance:
        if args.distance>65535 or args.distance < 1:
        	print 'ERROR: argument -d/--distance should be between 1 and 65535'
        	exit(-1)
if args.coverage:
	if args.coverage > 65535 or args.coverage < 1:
		print 'ERROR: argument -c/--coverage should be between 1 and 65535'
		exit(-1)
if args.longread:
	longread_path = args.longread
	

#start_time = datetime.datetime.now()
#print start_time

# Step 1 detect misassembled errors by reference  ####################
print'''/////STEP1 STARTED //////////////////////////////////////////////////////////////////////////////////////////////////'''

if os.path.exists(output_dir + '/ref'):
        shutil.rmtree("./output")

if not os.path.exists(output_dir):
        os.makedirs(output_dir)

if not os.path.exists(temp_dir):
	os.makedirs(temp_dir)


if os.path.exists(output_dir + '/ref'):
	print 'ERROR: ' + output_dir + '/ref' + ' already exist, please delete it before running ref'
	exit(-1)
else:
	os.mkdir(output_dir + '/ref')
        os.mkdir(output_dir+'/ref/aln')
str=reference_path
list=str.split('/')
reference_name=list[len(list)-1]
shutil.copyfile(reference_path,output_dir+'/ref/aln/'+reference_name)
bwa_index_command = 'bwa index '  + output_dir +'/ref/aln/' + reference_name + ' 2>' + output_dir+'/ref/aln/bwa.txt' 

print 'Running command:' + bwa_index_command
err = os.system(bwa_index_command)
if err != 0:
        print 'ERROR: ' + 'Failed to run bwa index:' + os.strerror(err)
        exit(-1)
list1 = contig_path.split('/')
contig_name = list1[len(list1)-1]
tempstr = contig_name.split('.')
sam_name = tempstr[0] + '.sam'
bwa_mem_command = 'bwa mem -a ' +  output_dir +'/ref/aln/' + reference_name + ' ' + contig_path + ' > ' + output_dir + '/ref/aln/' + sam_name 
bwa_mem_command += ' 2>'+output_dir+'/ref/aln/bwa.txt'
 
print 'Running command:' + bwa_mem_command
err = os.system(bwa_mem_command)
if err !=0:
        print 'ERROR: ' + 'Failed to run bwa mem : ' + os.strerror(err)
        exit(-1)
#print '\n'
#print 'Detect misassemble errors by reference genome ' 

REF_command = 'processRefGenome ' + output_dir + '/ref/aln/' + reference_name  + ' ' +  output_dir + '/ref/aln/' + sam_name  + ' ' + output_dir + '/ref/locations.txt' 
REF_command+= ' ' + output_dir + '/ref/subContig.fa'
print 'Running command:' + REF_command
err = os.system(REF_command)
if err !=0:
       print 'ERROR: ' + 'Failed ro run ref ' + os.strerror(err)
       exit(-1)
#print '\n'

if not os.path.exists(output_dir):
       os.makedirs(output_dir)

if not os.path.exists(temp_dir):
       os.makedirs(temp_dir)

if os.path.exists(output_dir + '/short'):
       print 'ERROR: ' + output_dir + '/ref' + ' already exist, please delete it before running ref '
       exit(-1)
else:
       os.mkdir(output_dir + '/short')
       os.mkdir(output_dir+'/short/aln')
shutil.move(output_dir+'/ref/subContig.fa',output_dir+ '/short/aln/')
shutil.move(output_dir + '/ref/locations.txt',output_dir + '/short/aln/')
bowtie2_index_command= 'bowtie2-build  ' + output_dir + '/short/aln/subContig.fa ' + output_dir + '/short/aln/subContig' + ' 1> '+output_dir + '/short/aln/bowtie2.txt'
bowtie2_index_command+=' 2> ' +output_dir + '/short/aln/bowtie2.err'
print 'Running command:' + bowtie2_index_command
err = os.system(bowtie2_index_command)
if err!=0:
      print 'ERROR: ' + 'Failed to run bowtie2 build : ' +os.strerror(err)
      exit(-1)
bowtie2_mapread1_command = 'bowtie2 -x ' + output_dir + '/short/aln/subContig -U ' + short_read1_path + ' -S ' + output_dir 
bowtie2_mapread1_command+= '/short/aln/read1.sam' +' 1>> ' + output_dir +'/short/aln/bowtie2.txt'  +' 2>> '+output_dir+'/short/aln/bowtie2.err'
print 'Running command:' + bowtie2_mapread1_command
err = os.system(bowtie2_mapread1_command)
if err !=0:
      print 'ERROR: ' + ' Failed to run bowtie2   read1 : ' + os.strerror(err)
      exit(-1)

bowtie2_mapread2_command = 'bowtie2 -x ' + output_dir + '/short/aln/subContig -U ' + short_read2_path 
bowtie2_mapread2_command+= ' -S ' + output_dir + '/short/aln/read2.sam ' +' 1>> ' +output_dir + '/short/aln/bowtie2.txt' +' 2>> '+output_dir+'/short/aln/bowtie2.err'
print 'Running command:' + bowtie2_mapread1_command
err = os.system(bowtie2_mapread2_command)
if err !=0:
      print 'ERROR: ' + ' Failed to run bowtie2   read2 : ' + os.strerror(err)
      exit(-1)
#print '\n'
miserr=3
print'''/////STEP1 DONE//////////////////////////////////////////////////////////////////////////////////////////////////'''

print'''/////STEP2 STARTED//////////////////////////////////////////////////////////////////////////////////////////////////'''
short_command = 'processShortReads ' + output_dir + '/short/aln/read1.sam'  + ' ' + output_dir + '/short/aln/read2.sam' + ' ' 
short_command += output_dir + '/short/aln/subContig.fa '
short_command += output_dir + '/short/aln/locations.txt '
short_command+= ' ' + output_dir + '/short/misLocations.txt ' + output_dir + '/short/splitContigs.fa ' +output_dir + '/short/unsplitContigs.fa'
print 'Running command:' + short_command
err = os.system(short_command)
if err !=0:
       print 'ERROR: ' + 'Failed ro run shortread ' + os.strerror(err)
       exit(-1)
print'''/////STEP2 DONE //////////////////////////////////////////////////////////////////////////////////////////////////'''
shutil.copy2(output_dir+'/short/misLocations.txt','./output/')
if args.longread:
      print'''/////STEP3 STARTED//////////////////////////////////////////////////////////////////////////////////////////////////'''
      
      if not os.path.exists(output_dir):
             os.makedirs(output_dir)

      if not os.path.exists(temp_dir):
             os.makedirs(temp_dir)

      if os.path.exists(output_dir + '/long'):
             print 'ERROR: ' + output_dir + '/long' + ' already exist, please delete it before running ref'
             exit(-1)
      else:
             os.mkdir(output_dir + '/long')
             os.mkdir(output_dir+'/long/aln')
      long_list=longread_path.split('/')
      longread_name=long_list[len(long_list)-1]
      shutil.copyfile(longread_path,output_dir + '/long/aln/' +  longread_name)
      longread_bwa_index_command='bwa index ' + output_dir + '/long/aln/' + longread_name + ' 2> ' +output_dir+'/long/aln/bwa.txt'
      print 'Running command:' + longread_bwa_index_command
      err = os.system(longread_bwa_index_command)
      if err != 0:
             print 'ERROR: ' + 'Failed to run bwa index : ' + os.strerror(err)
             exit(-1)

      long_sam_name='long_' + sam_name
      longread_bwa_mem_command ='bwa mem -a ' + output_dir + '/long/aln/' + longread_name
      longread_bwa_mem_command +=' ' + contig_path + '> '  + output_dir + '/long/aln/' + long_sam_name  +' 2>> ' +output_dir + '/long/aln/bwa.txt'
      print 'Running command:' + longread_bwa_mem_command
      err = os.system(longread_bwa_mem_command)
      if err != 0:
             print 'ERROR: ' + 'Failed to run bwa mem : ' + os.etrerror(err)
             exit(-1)
#      print '\n'
      print'''/////STEP3 DONE//////////////////////////////////////////////////////////////////////////////////////////////////'''

      print'''/////STEP4 STARTED //////////////////////////////////////////////////////////////////////////////////////////////////'''
      long_command = 'processLongReads ' + output_dir + '/long/aln/' + long_sam_name  + ' ' + output_dir + '/long/aln/' 
      long_command += longread_name  + ' ' + output_dir + '/long/misLocations.txt' 
      long_command += ' ' + output_dir + '/long/splitContigs.fa '
      long_command += output_dir + '/long/unsplitContigs.fa'
      print 'Running command:' + long_command
      err = os.system(long_command)
      if err !=0:
             print 'ERROR: ' + 'Failed ro run  longread ' + os.strerror(err)
             exit(-1)
    #  print '\n'
      print'''/////STEP4 DONE //////////////////////////////////////////////////////////////////////////////////////////////////'''
      print '\n'

