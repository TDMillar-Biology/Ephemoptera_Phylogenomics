#!/usr/bin/env python


# Script was written by Jesse Breinholt (jessebreinholt@gmail.com)
# assume linux command line tools, usearch, ncbi blast, SPADES assembler in path


# UF cluster---- module load gcc/5.2.0 spades python ncbi_blast usearch



import os, sys
import gzip
from Bio import SeqIO

mtemp="mkdir %s"


error_c="spades.py --only-error-correction -o %s -1 %s -2 %s -t %s"


Ucat="zcat %s > %s"
Us="usearch -ublast %s -db %s -evalue 0.01 -id 0.65 -maxaccepts 1 -strand both -blast6out %s"
Us2="usearch -ublast %s -db %s -evalue 0.0001 -id 0.95 -maxaccepts 1 -strand both -blast6out %s"

Ufilter="cut -f1 -d \" \" %s | sort -u | awk \'{print \"@\"$1}\' > %s"
Ufilter2="cat %s %s | cut -f1 -d \" \" | sort -u | awk \'{print \"@\"$1}\' > %s"

ucollapse="usearch -cluster_fast %s -id 0.98 -query_cov 0.98 -strand both -sort length -centroids %s"

rlist="grep -v \"^@\" %s | cut -f1 | sort -u > %s"

getread="zgrep --no-group-separator -A3 -F -f %s %s | gzip > %s"
ass="spades.py --only-assembler -o %s -1 %s -2 %s -s %s -t %s"
ass2="spades.py --only-assembler --careful -o %s -1 %s -2 %s -s %s -t %s"
blastdb="makeblastdb -dbtype \"nucl\" -in %s -out %s"
realblast= "blastn -query %s -db %s -out %s -outfmt 6 -evalue 0.0001 -qcov_hsp_perc 80 -num_threads %s"
fblast= "blastn -query %s -db %s -out %s -outfmt 6 -evalue 0.001 -num_threads %s"

rblast= "tblastx -query %s -db %s -out %s -outfmt 6 -evalue 0.01 -num_threads %s"
rblast2= "tblastx -query %s -db %s -out %s -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe\" -evalue 0.0001 -qcov_hsp_perc 0.50 -num_threads %s"


rusearch="usearch -ublast %s -db %s -evalue 0.001 -id 0.70 -mincols 100 -blast6out %s  -strand both"
fusearch="usearch -ublast %s -db %s -evalue 0.001 -id 0.70 -mincols 100 -blast6out %s  -strand both"

filter1="sort -k12,12rn %s | sort -u -k1,2 > %s"
filter2="sort -k12,12rn %s | sort -u -k1,1 > %s"
joincat="cat %s %s > %s"
sc1="cut -f1 %s | sort -u > %s"
fasta_filter="grep --no-group-separator -A1 -F -f %s %s > %s"
key_filter="grep --no-group-separator -F -f %s %s > %s"

kkill="rm -r %s"

remove = set()
blastdict= {}
locidict={}
seq = []
lista = []
arguments = sys.argv

if len(arguments) == 1:
	sys.exit("try python genome_getprobe_BLAST.py -h for help")
try:
	hflag = arguments.index("-h")
except:
	hflag = None

if hflag:
	sys.exit("\n\n#################################################################\n\n./BASS.py -inp REF_Probe.fasta -ing genome.fasta -threads ## -tname taxaname -flanksize ## -1 read1 -2 read2  -rname Reference name in target file \n\n\n\t\t-inp: fasta file of probe sequences\n\t\t-ing: genome or transcriptome file with sequnces data for each sequnces on a single line\n\t\t-threads: number of threads to use in blast\n\t\t-tname: taxa name you want to put on the output sequnces\n\t\t-flanksize: number of base pairs you want inclulded before or after the probe hit\n\n#################################################################\n\n")

try:
	probeset = arguments[arguments.index("-inp")+1]
except:
	sys.exit("Error: need probe fasta")

try:
	genome = arguments[arguments.index("-ing")+1]
except:
	sys.exit("Error: need genome fasta")

try:
	taxaname = arguments[arguments.index("-tname")+1]
except:
	sys.exit("Error: need taxa name")


try:
	R1 = arguments[arguments.index("-1")+1]
except:
	sys.exit("Error: need read 1")
try:
	R2 = arguments[arguments.index("-2")+1]
except:
	sys.exit("Error: need read 2")

try:
	threads = int(arguments[arguments.index("-threads")+1])
except:
	threads = 1

try:
	flank = int(arguments[arguments.index("-flanksize")+1])
except:
	flank = 0

try:
	probeset = arguments[arguments.index("-inp")+1]
except:
	sys.exit("Error: need probe fasta")
try:
	refname= arguments[arguments.index("-rname")+1]
except:
	sys.exit("Error: need ref name")

try:
	fullprobeset = arguments[arguments.index("-inpf")+1]
except:
	fullprobeset = arguments[arguments.index("-inp")+1]

#S1 strip file name from path
if "/" in R1:
	r1part=R1.split("/") ## split file path at /
	rout1=r1part[-1] ## grab the last part (base name)
elif "/" not in R1:
	rout1=R1
	
if "/" in R2:
	r2part=R2.split("/")
	rout2=r2part[-1]
elif "/" not in R2:
	rout2=R2

prout1=rout1.strip(".gz") ## remove .gz from basename
prout2=rout2.strip(".gz")


#S2 make taxa dir
print("\n\n making tmp dir: "+ taxaname+"_tmp\n\n")
os.system(mtemp % (taxaname+"_tmp")) ## Make tmp directory to write to 
print("Done with making the taxa directory")

#S3 error correction
print("\n\nIllumina Error correction using Spades\n\n")
os.system( error_c % (taxaname+"_tmp/corrected/",R1,R2,threads)) ## error correction (find illumina seq errors) this uses BayesHmmer
os.system( Ucat % (taxaname+"_tmp/corrected/corrected/*.gz",taxaname+"_tmp/BOTH.fq" )) ## Concatinate forward and reverse reads
os.system('>&2 echo "S3"')




#S4 usearch mapping 1 fullprobeset
print("\n\nUsearch mapping......\n\n")
os.system( Us % (taxaname+"_tmp/BOTH.fq", fullprobeset, taxaname+"_tmp/MAP1")) ## Usearch searches the full probeset for sequence allignments and write that to /MAP1
os.system( Ufilter % (taxaname+"_tmp/MAP1",taxaname+"_tmp/mreads")) ## cut field 1 using delimeter ""
os.system('>&2 echo "S4"')

#S5 get ontarget reads
print("\n\nCollecting mapped reads\n\n")
os.system(getread % (taxaname+"_tmp/mreads", taxaname+"_tmp/corrected/corrected/"+prout1+".00.0_0.cor.fastq.gz", taxaname+"_tmp/sub_"+rout1)) ## I fixed up prout1 to be only basename
os.system(getread % (taxaname+"_tmp/mreads", taxaname+"_tmp/corrected/corrected/"+prout2+".00.0_0.cor.fastq.gz", taxaname+"_tmp/sub_"+ rout2))
os.system(getread % (taxaname+"_tmp/mreads", taxaname+"_tmp/corrected/corrected/*unpaired*.gz", taxaname+"_tmp/sub_unpaired_"+ rout2))
os.system('>&2 echo "S5"')


#S6 spades assembly 1
print("\n\n1st ASSEMBLY\n\n")
os.system(ass % (taxaname+"_tmp/ASS1",taxaname+"_tmp/sub_"+rout1,taxaname+"_tmp/sub_"+rout2,taxaname+"_tmp/sub_unpaired_"+ rout2, threads ))
os.system('>&2 echo "S6"')


#S7 ASSEMBLY on Target screen to REF
print("\n\ASSEMBLY on Target screen\n\n")
os.system(rusearch % (taxaname+"_tmp/ASS1/scaffolds.fasta", probeset, taxaname+"_tmp/screen_out"))
os.system(sc1 % (taxaname+"_tmp/screen_out", taxaname+"_tmp/pass_out"))
os.system('>&2 echo "S7"')


#try tblastx screen
os.system(blastdb % (probeset, taxaname+"_tmp/"+probeset+"_db"))
os.system(rblast % (taxaname+"_tmp/ASS1/scaffolds.fasta", taxaname+"_tmp/"+probeset+"_db", taxaname+"_tmp/screen_out", threads ))
os.system(sc1 % (taxaname+"_tmp/screen_out", taxaname+"_tmp/pass_out"))
os.system('>&2 echo "S8a"')

##S8 write asssebly 1 to single line
mfasta=open(taxaname+"_tmp/ASS1/sl_scaffolds.fasta", "w")
lc=0
with open(taxaname+"_tmp/ASS1/scaffolds.fasta") as f:
	for line in f:
		if line.startswith(">"):
			if lc==0:
				mfasta.write(line)
				lc+=1
			elif lc !=0:
				mfasta.write("\n"+line)
		else:	
			mfasta.write(line.strip("\n"))
mfasta.close()
os.system('>&2 echo "S8b"')

#S9 pull seqs that pass filter
os.system(fasta_filter % (taxaname+"_tmp/pass_out", taxaname+"_tmp/ASS1/sl_scaffolds.fasta",taxaname+"_tmp/ASS1/clean_scaffolds.fasta"))
os.system('>&2 echo "S9"')

#S10 mapping raw read to on target
print("MAPPING 2......")
os.system( Us2 % (taxaname+"_tmp/BOTH.fq", taxaname+"_tmp/ASS1/clean_scaffolds.fasta" , taxaname+"_tmp/MAP2"))
os.system( Ufilter2 % (taxaname+"_tmp/MAP2",taxaname+"_tmp/MAP1",taxaname+"_tmp/mreads2"))
os.system('>&2 echo "S10"')


#S11 Collecting mapped reads for assembly2
print("\n\nCollecting mapped reads\n\n")
os.system(getread % (taxaname+"_tmp/mreads2", taxaname+"_tmp/corrected/corrected/"+prout1+".00.0_0.cor.fastq.gz", taxaname+"_tmp/sub2_"+rout1))
os.system(getread % (taxaname+"_tmp/mreads2", taxaname+"_tmp/corrected/corrected/"+prout2+".00.0_0.cor.fastq.gz", taxaname+"_tmp/sub2_"+ rout2))
os.system(getread % (taxaname+"_tmp/mreads2", taxaname+"_tmp/corrected/corrected/*unpaired*.gz", taxaname+"_tmp/sub2_unpaired_"+ rout2))
os.system('>&2 echo "S11"')

#S12 assembly2
print("\n\nASSEMBLY\n\n")
os.system(ass2 % (taxaname+"_tmp/ASS2",taxaname+"_tmp/corrected/corrected/"+prout1+".00.0_0.cor.fastq.gz",taxaname+"_tmp/corrected/corrected/"+prout2+".00.0_0.cor.fastq.gz",taxaname+"_tmp/corrected/corrected/*unpaired*.gz", threads ))
os.system('>&2 echo "S12"')

#S14 Collapsing ASSEMBLY [id:0.98, query_cov: 0.98]
print("\n\nCollapsing ASSEMBLY [id:0.98, query_cov: 0.98]\n\n")
os.system( ucollapse % (taxaname+"_tmp/ASS2/scaffolds.fasta",taxaname+"_tmp/ASS2/rscaffolds.fasta"))


#S14 making blast databases
print("making blast databases")
os.system(blastdb % (genome, taxaname+"_tmp/"+genome+"_db"))
os.system(blastdb % (probeset, taxaname+"_tmp/"+probeset+"_db"))
os.system(blastdb % (taxaname+"_tmp/ASS2/rscaffolds.fasta", taxaname+"_tmp/ASS2/rscaffolds.fasta_db"))

#S15 blast and find probe location and filter for best location
os.system(realblast % (probeset,taxaname+"_tmp/"+genome+"_db",taxaname+"_tmp/prelocal_out", threads))
os.system(filter2 % (taxaname+"_tmp/prelocal_out", taxaname+"_tmp/local_out"))





#S16 ASSEMBLY on Target screen
print("\n\ASSEMBLY on Target screen\n\n")
os.system(fusearch % (probeset,taxaname+"_tmp/ASS2/rscaffolds.fasta", taxaname+"_tmp/screen2_out"))
os.system(filter1 % (taxaname+"_tmp/screen2_out", taxaname+"_tmp/clean_screen2_out"))

os.system(rblast2 % (probeset,taxaname+"_tmp/ASS2/rscaffolds.fasta_db",taxaname+"_tmp/screen2_out",threads))
os.system(filter1 % (taxaname+"_tmp/screen2_out", taxaname+"_tmp/clean_screen2_out"))


#S17 make single line fasta
mfasta2=open(taxaname+"_tmp/ASS2/sl_rscaffolds.fasta", "w")
lc=0
with open(taxaname+"_tmp/ASS2/rscaffolds.fasta") as f:
	for line in f:
		if line.startswith(">"):
			if lc==0:
				mfasta2.write(line)
				lc+=1
			elif lc !=0:
				mfasta2.write("\n"+line)
		else:	
			mfasta2.write(line.strip("\n"))
mfasta2.close()




#S18 Identifying Loci and trimming to target region
print("\n\nIdentifying Loci and trimming to target region\n\n")
blastfile=open(taxaname+"_tmp/clean_screen2_out")



print("\n\nMaking hit dictionary\n\n")
#qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qframe

for line in blastfile:
	line = line.strip().split()
	if line[1] in blastdict.keys():
		lista = blastdict[line[1]]
		locilist = line[0].split("_")
		lista.append(locilist[0]+"_"+line[8]+"_"+line[9]+"_"+line[11]+"_"+line[12])
		blastdict[line[1]] = lista
	else:
		lista = []
		locilist = line[0].split("_")
		lista.append(locilist[0]+"_"+line[8]+"_"+line[9]+"_"+line[11]+"_"+line[12])
		blastdict[line[1]] = lista
		
print("\n\nGetting sequences\n\n")

tb=open(taxaname + "_table.Key", "w")
RCLIST=open(taxaname + "_RC.list", "w")
FL=open(taxaname + "_targetsFULL.fasta", "w")

with open(taxaname + "_targets.fasta", "w") as out:
	with open(taxaname+"_tmp/ASS2/sl_rscaffolds.fasta") as GG:
		line=GG.readline()
		while line:
			if line[0] == ">": ## If the first char in line is > (its the title)
				try:
					line = line.split() ## split str into list
					id = line[0].lstrip(">") ## get rid of '>' on the line and assign that to id
				except:
					id = line.strip().lstrip(">")
			if id in blastdict.keys(): ## if the line is in the blast dictionary
				line=GG.readline()
				for info in blastdict[id]:
					#a loci b d (bit score) e(qframe)
					a, b, c, d, e= info.split("_")
					if a in locidict:
						locidict[a] += 1
					else:
						locidict[a] = 1
					b = int(b)
					c = int(c)
					e = int(e)
					
					if e < 0:
						RCLIST.write(a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"|"+id+"\n")
					
					if c > b:
						if b > flank:
							start_rec = b - 1 - flank
						else:
							start_rec = 0
						end_rec = c + flank
						out.write(">"+a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"\n")
						out.write(line[start_rec:end_rec].strip()+"\n")
						tb.write(id+"\t"+a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"\n")
						FL.write(">"+a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"|"+id+"\n"+line)
					
					elif b > c:
					#backwards coordinated given
						if c > flank:
							start_rec = c - 1 - flank
						else:
							start_rec = 0
						end_rec = b + flank
						out.write(">"+a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"\n")
#						out.write(line[start_rec:end_rec].strip()[::-1]+"\n")
						out.write(line[start_rec:end_rec].strip()+"\n")
						tb.write(id+"\t"+a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"\n")
						FL.write(">"+a+"_"+taxaname+"_comp"+str(locidict[a])+"_"+d+"|"+id+"\n"+line)
						
				line = GG.readline()
			else:
				line=GG.readline()
				line=GG.readline()	
tb.close()
FL.close()

#S19 NCBI TBLASTX for Ortholog filter

print("\n\nNCBI TBLASTX for Ortholog filter\n\n")

os.system(rblast % (taxaname + "_targets.fasta",taxaname+"_tmp/"+genome+"_db",taxaname+"_tmp/preFass_out", threads))
os.system(filter2 % (taxaname+"_tmp/preFass_out", taxaname+"_tmp/Fass_out"))
os.system(joincat % (taxaname+"_tmp/local_out",taxaname+"_tmp/Fass_out", taxaname+"_tmp/Ortho.in"))

#S20 Ortholog filter


print("\n\nNCBI TBLASTX for Ortholog filter Breinholt et al 2018\n\n")
outfile1 =open(taxaname + "_del_list.txt", "w")
outfile2 =open(taxaname + "_keep_list", "w")

hit=[]
loci=set([])
count=0
seqset=set([])
#open table and parse for scaford hit and coordinates of reference
with open(taxaname+"_tmp/Ortho.in", "r") as table:
	makein=table.readlines()
	for i in makein:
		loci.add(i.split("_",-1)[0])
		ALL_loci=list(loci)
table.close()


for x in ALL_loci:
	print("Processing " + x + " .............\n")
	with open(taxaname+"_tmp/Ortho.in", "r") as table2:
		makein2=table2.readlines()
		for i in makein2:
			taxa, scaf, id, length, mismatch, gaps, qs, qend, ts, tend, evalue, bit=i.split()
			if taxa.startswith(str(x) + "_"):
				if taxa == str(x) + "__" + refname + "_R":
					hit.append(scaf)
					print(taxa + " scaffold : " + scaf)
					leftcoord=int(ts)
					rightcoord=int(tend)
					if int(ts) < int(tend):
						direction= int(1)
					if int(ts) > int(tend):
						direction = int(0)
	table2.close()		


# open again as diffrent names to start at the top check to see if it hit scaf and coordinates
	with open(taxaname+"_tmp/Ortho.in", "r") as table3:
		makein3=table3.readlines()
		for i in makein3:
			taxa3, scaf3, id3, length3, mismatch3, gaps3, qs3, qend3, ts3, tend3, evalue3, bit3=i.split()
			seqset.add(taxa3)
			if int(ts3) < int(tend3):
				seqdirection= int(1)
			if int(ts3) > int(tend3):
				seqdirection = int(0)
#			print("seq direction(fwd=1, rev=0) : " + str(seqdirection))
			if taxa3.startswith(str(x) + "_"):
				if scaf3 not in hit:
					print("diffent scaffold " + taxa3 + " scaffold : " + scaf3)
					outfile1.write(taxa3 + "\n")
					count +=1
				if scaf3 in hit: 
					if direction is 1 and seqdirection is 1:
						if int(ts3) < rightcoord and int(tend3) > leftcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print("Same scaffold diffrent location Direction (ref fwd : seq fwd) " + str(direction) + ":"+ str(seqdirection))
							print(str(leftcoord) + " " + str(rightcoord) + "|" + i)
							count +=1
					
					if direction is 1 and seqdirection is 0:
						if int(tend3) < rightcoord and int(ts3) > leftcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print("Same scaffold diffrent location Direction(ref fwd: seq rev) " + str(direction) + ":"+ str(seqdirection))
							print(str(leftcoord) + " " + str(rightcoord) + "|" + i)
							count +=1					
					if direction is 0 and seqdirection is 0:
						if int(tend3) < leftcoord and int(ts3) > rightcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print("Same scaffold diffrent location Direction(ref rev: seq rev) " + str(direction) + ":"+ str(seqdirection))
							print(str(leftcoord) + " " + str(rightcoord) + "|" + i)
							count +=1
					if direction is 0 and seqdirection is 1:
						if int(ts3) < leftcoord and int(tend3) > rightcoord:
							outfile2.write(taxa3 + "\n")
						else:
							outfile1.write(taxa3 + "\n")
							print("Same scaffold diffrent location Direction(ref rev: seq fwds) " + str(direction) + ":"+ str(seqdirection))
							print(str(leftcoord) + " " + str(rightcoord) + "|" + i)
							count +=1


print("Ortholog filter results: "+str(count) + "/" + str(len(seqset)) + " (delete/total)")

table3.close()
outfile1.close()
outfile2.close()

#S21 Generating ortholog fasta files

print("\n\nGenerating ortholog fasta files\n\n")

os.system(fasta_filter % (taxaname+"_keep_list",taxaname+"_targets.fasta",taxaname+"_targets_ORTHO.fasta"))
os.system(fasta_filter % (taxaname+"_keep_list",taxaname+"_targetsFULL.fasta",taxaname+"_targetsFULL_ORTHO.fasta"))
os.system(key_filter % (taxaname+"_keep_list",taxaname+"_table.Key",taxaname+"_keep_list_ORTHO"))



#S22 CLean
print("\n\nclearing tmp files\n\n")
os.system( kkill % (taxaname+"_tmp"))

print("\n\n!!!!!Complete!!!!!\n\n")


sys.exit()