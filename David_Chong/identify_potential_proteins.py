#!/usr/bin/env python
'''
Reads the mass spec fragments from percolator and mascot and filters them against the mouse proteome to identify a set of fragments that do not map to known proteins.
Genomic coordinates of these fragments are then identified using blast and some info about position relative to exons/other genes are noted.
An IGV script is created for each fragment to allow automated visualisation of the alignment in IGV.
'''
import sys, subprocess, os, requests
import pandas as pd

server = "https://rest.ensembl.org"
data_folder = "/media/hdd/Prabakaran_Lab/data"
executable_folder = "/media/hdd/Prabakaran_Lab/bin"
Tcellnucdb = data_folder + "/Tcell_nucleotide_database_ce_29-04-2018.fa"
Bcellnucdb = data_folder + "/Bcell_nucleotide_database_ce_29-04-2018.fa"
Tcellbam = data_folder + "/T_cell_reads.bam"
Bcellbam = data_folder + "/B_cell_reads.bam"

def get_overlapping_features(species, chrom, start, end, strand = None, featurelist = ["gene"]):
	ext = "/overlap/region/" + species + "/" + chrom + ":" + start + "-" + end
	if strand != None:
		ext += ":" + strand
	ext += "?"
	for f in featurelist:
		ext += "feature=" + f + ";"
	ext = ext[:-1]
	r = requests.get(server + ext, headers = { "Content-Type" : "application/json" })
	return(r)
	
def classify_overlap(test_start,test_end,reference_start,reference_end,strand):
	if strand == 1:
		if test_end < reference_start:
			return("upstream")
		elif test_end == reference_end and test_start < reference_start:
			return("upstream_extension")
		elif test_end == reference_end and test_start > reference_start:
			return("upstream_truncation")
		elif test_end >= reference_start and test_end < reference_end and test_start < reference_start:
			return("upstream_overlap")
		elif test_start > reference_end:
			return("downstream")
		elif test_start == reference_start and test_end > reference_end:
			return("downstream_extension")
		elif test_start == reference_start and test_end < reference_end:
			return("downstream_truncation")
		elif test_start > reference_start and test_start <= reference_end and test_end > reference_end:
			return("downstream_overlap")
		elif test_start == reference_start and test_end == reference_end:
			return("exact-match")
		elif test_start < reference_start and test_end > reference_end:
			return("encompass")
		elif test_start > reference_start and test_end < reference_end:
			return("within")
		
	elif strand == -1:
		if test_start > reference_end:
			return("upstream")
		elif test_start == reference_start and test_end > reference_end:
			return("upstream_extension")
		elif test_start == reference_start and test_end < reference_end:
			return("upstream_truncation")
		elif test_start <= reference_end and test_start > reference_start and test_end > reference_end:
			return("upstream_overlap")
		elif test_end < reference_start:
			return("downstream")
		elif test_end == reference_end and test_start < reference_start:
			return("downstream_extension")
		elif test_end == reference_end and test_start > reference_start:
			return("downstream_truncation")
		elif test_end >= reference_start and test_end < reference_end and test_start < reference_start:
			return("downstream_overlap")
		elif test_start == reference_start and test_end == reference_end:
			return("exact-match")
		elif test_start < reference_start and test_end > reference_end:
			return("encompass")
		elif test_start > reference_start and test_end < reference_end:
			return("within")

#read in list of IDs
if len(sys.argv) > 2:
	f = open(sys.argv[1])
	data = f.read()
	f.close()
	#prepare output file format
	f = open(sys.argv[1].split(".")[0] + "_potential_peptides.csv","w")
	f.write("prot_seq,chrom,prot_start,prot_end,evidence_start,evidence_end,e_val\n")
	data = data.split("\n")[:-1]
	#ensure correct database is set and directory structure
	if not os.path.exists(data_folder + "/Potential_Novel_Peptide_ORFs"):
		subprocess.run(["mkdir", data_folder + "/Potential_Novel_Peptide_ORFs"])
	if sys.argv[1].find("T") == 0:
		nucdb = Tcellnucdb
		bamdb = Tcellbam
		celltype = "T"
	elif sys.argv[1].find("B") == 0:
		nucdb = Bcellnucdb
		bamdb = Bcellbam
		celltype = "B"
	else:
		exit()
	print("Using databases: " + nucdb + " , " + bamdb)
	#Prepare output file format
	if os.path.exists(data_folder + "/novel_peptide_region_features.csv"):
		genomicfeatures = open(data_folder + "/novel_peptide_region_features.csv","a")
	else:
		genomicfeatures = open(data_folder + "/novel_peptide_region_features.csv","w")
		genomicfeatures.write("region,fragment,parent_id,feature_type,feature_id,feature_name,strand_relationship,coordinate_relationship,biotype_relationship,translation_relationship,feature_start,feature_end\n")
	if os.path.exists(data_folder + "/novel_peptide_frame_evidence.csv"):
		frameevidence = open(data_folder + "/novel_peptide_frame_evidence.csv","a")
	else:
		frameevidence = open(data_folder + "/novel_peptide_frame_evidence.csv","w")
		frameevidence.write("region,source,frame,fragment,evalue,alignment_start,alignment_end,classification,predicted_aa_seq,predicted_aa_start,predicted_aa_end\n")
	if not os.path.exists(data_folder + "/Novel_Peptide_ORFs_Fragment_Alignment_Info.csv"):
		ORF_Fragment_Info = open(data_folder + "/Novel_Peptide_ORFs_Fragment_Alignment_Info.csv","w")
		ORF_Fragment_Info.write("ORF_id,fragment,alignment_start,alignment_end,evalue\n")
	else:
		ORF_Fragment_Info = open(data_folder + "/Novel_Peptide_ORFs_Fragment_Alignment_Info.csv","a")
	if not os.path.exists(data_folder + "/Novel_Peptide_ORFs_Info.csv"):
		ORF_Info = open(data_folder + "/Novel_Peptide_ORFs_Info.csv","w")
		ORF_Info.write("ORF_region,ORF_chrom,ORF_strand,ORF_start,ORF_end,predicted_ORF_aa_seq,peptide_coverage,number_of_peptides,read_coverage,Cell_Type,regionpath,ORFpath,Classification\n")
	else:
		ORF_Info = open(data_folder + "/Novel_Peptide_ORFs_Info.csv","a")
	if os.path.exists(executable_folder + "/novel_peptide_regions_igv.bat"):
		igvbat = open(executable_folder + "/novel_peptide_regions_igv.bat","a")
	else:
		igvbat = open(executable_folder + "/novel_peptide_regions_igv.bat","w")
		igvbat.write("new\ngenome mm10\nload " + Tcellbam + "\nload " + Bcellbam + "\nsnapshotDirectory /media/hdd/Prabakaran_Lab/igv_snapshots\n")
	#begin analysing sequences
	for d in data:
		print("Processing " + d)
		#extract region info from header
		temp = d.split(":")
		pattern = temp[1].split("-")[2].split("_")[0]
		chrom = temp[0][temp[0].find("chr") + 3:]
		if temp[2] == "+":
			strand = 1
		else:
			strand = -1
		temp2 = temp[1].split("-")
		if d.find("_") != -1:
			checkcoords = subprocess.run(" ".join(["blastn", "-db", data_folder + "/blastdb/mm10.fa", "-query", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fna", "-perc_identity", "100.0", "-qcov_hsp_perc", "100.0", "-outfmt", "\"6", "qseqid", "sseqid", "qlen", "length", "nident", "gaps", "bitscore", "evalue", "qstart", "qend", "sstart", "send\""]), shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[:-1]
			print(checkcoords)
			start = checkcoords[0].split("\t")[10] if checkcoords[0].split("\t")[10] < checkcoords[0].split("\t")[11] else checkcoords[0].split("\t")[11] 
			end = checkcoords[0].split("\t")[11] if checkcoords[0].split("\t")[10] < checkcoords[0].split("\t")[11] else checkcoords[0].split("\t")[10] 
		else:
			start = temp2[0]
			end = temp2[1]
		region_length = int(end) - int(start)
		protid = temp2[2]
		#Extract nucleotide sequence from database based on ID
		print("Extracting nucleotide sequence")
		if not os.path.exists(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fna"):
			seq = open(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fna", "w")
			if temp[2] == "+":
				seq.write(subprocess.run(["pcregrep", "--buffer-size", "320000", "-M", ">" + temp[0] + ":" + temp[1] + ":\\" + temp[2] + "\\n[ACTG]+", nucdb],stdout=subprocess.PIPE).stdout.decode('utf-8'))
			else:
				seq.write(subprocess.run(["pcregrep", "--buffer-size", "320000", "-M", ">" + temp[0] + ":" + temp[1] + ":" + temp[2] + "\\n[ACTG]+", nucdb],stdout=subprocess.PIPE).stdout.decode('utf-8'))
			seq.close()
		#Finding ORFs
		
		#if not os.path.exists(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".ORFfinder.faa"):
		#	subprocess.run([executable_folder + "/ORFfinder", "-s", "2", "-in", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fna", "-out", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".ORFfinder.faa"])
		#	subprocess.run(["makeblastdb", "-in", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".ORFfinder.faa", "-dbtype", "prot"])

		if not os.path.exists(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".transeq.faa"):
			subprocess.run(["transeq", "-frame", "6", "-sequence", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fna", "-outseq", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".transeq.faa"])
			subprocess.run(["makeblastdb", "-in", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".transeq.faa", "-dbtype", "prot"])

		if not os.path.exists(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa"):
			peptide_fragments = list(set(subprocess.run(["grep", d, sys.argv[2]], stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[:-1]))
			frags = open(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa","w")
			counter = 1
			uniquefrags = {}
			for frag in peptide_fragments:
				if frag.split("\t")[14] not in uniquefrags.keys():
					uniquefrags[frag.split("\t")[14]] = [frag.split("\t")[0]]
				else:
					if frag.split("\t")[0] not in uniquefrags[frag.split("\t")[14]]:
						uniquefrags[frag.split("\t")[14]].append(frag.split("\t")[0])
			for frag, sources in uniquefrags.items():
				frags.write(">fragment_" + str(counter) + ":" + "-".join(sources) + "\n" + frag + "\n")
				counter += 1
			frags.close()
			checkknownproteins = subprocess.run(" ".join(["blastp", "-db", data_folder + "/mus_musculus_proteome/uniprot-proteome%3AUP000000589.fasta", "-query", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa", "-outfmt", "\"6", "qseqid", "sseqid", "qlen", "length", "nident", "gaps", "bitscore", "evalue", "qstart", "qend", "sstart", "send\""]), shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[:-1]
			sigalignments = []
			for a in checkknownproteins:
				print(a)
				a_field = a.split("\t")
				if abs(a_field[2] - a_field[3]) <= 1 and abs(int(a_field[2]) - int(a_field[4])) <= 1 and a_field[0] not in sigalignments and int(a_field[5]) == 0:
					sigalignments.append(a_field[0])
			frags = open(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa")
			fragdata = frags.read().split("\n")[:-1]
			frags.close()
			print(sigalignments)
			if len(sigalignments) == len(fragdata) / 2:
				pass
			else:
				frags = open(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa","w")
				for fragline in range(0,len(fragdata),2):
					if fragdata[fragline][1:] not in sigalignments:
						frags.write(fragdata[fragline] + "\n" + fragdata[fragline+1] + "\n")
				frags.close()

		#Identify surrounding genomic environment
		#Using ORFfinder
		'''
		alignments = subprocess.run(["blastp", "-db", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".ORFfinder.faa", "-query", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa", "-outfmt", "6"], stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[:-1]
		orf_regions = {}
		for a in alignments:
			temp3 = a.split("\t")
			if temp3[1] not in orf_regions.keys():
				orf_regions[temp3[1]] = [a]
			else:
				orf_regions[temp3[1]].append(a)
		for k,v in orf_regions.items():
			significant_alignment = False
			for a in v:
				temp3 = a.split("\t")
				if float(temp3[10]) <= 0.01:
					significant_alignment = True
					break
			if significant_alignment:
				ORFid = k.split("|")[1]
				print(ORFid)
				ORF_start = int(ORFid.split(":")[3])
				ORF_end = int(ORFid.split(":")[4])
				aaseq = "".join(subprocess.run(["pcregrep", "-M", ORFid.replace("+","\+") + " unnamed protein product(, partial)*\\n([A-Z]+\\n)+", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".ORFfinder.faa"], stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[1:-1])
				if  ORF_start < ORF_end:
					ORF_Info.write(ORFid + "," + chrom + ",+," + str(int(start) + ORF_start) + "," + str(int(start) + ORF_end) + "," + aaseq + ",")
					ORF_len = ORF_end-ORF_start
					ORF_strand = 1
				else:
					ORF_Info.write(ORFid + "," + chrom + ",-," + str(int(start) + ORF_end) + "," + str(int(start) + ORF_start) + "," + aaseq + ",")
					ORF_len = ORF_start-ORF_end
					ORF_strand = -1
				startstoppairs = []
				for a in v:
					temp3 = a.split("\t")
					ORF_Fragment_Info.write(ORFid + "," + temp3[0] + "," + temp3[8] + "," + temp3[9] + "," + temp3[10] + "\n")
					startstoppairs.append([int(temp3[8]),int(temp3[9])])
				contigs = []
				for s in startstoppairs:
					merged = False
					for c in contigs:
						description = classify_overlap(c[0],c[1],s[0],s[1],1)
						if description.find("overlap") != -1 or description.find("within") != -1 or description.find("encompass") != -1:
							if description.find("upstream") != -1:
								c[1] = s[1]
							elif description.find("downstream") != -1:
								c[0] = s[0]
							elif description == "within":
								c[0] = s[0]
								c[1] = s[1]
							merged = True
						else:
							pass
					if not merged:
						contigs.append(s)
				coveredaa = 0
				for c in contigs:
					coveredaa += c[1]-c[0]
				ORF_Info.write(str(coveredaa/len(aaseq)) + "," + str(len(v)) + ",")
				
				#Enter Bat script
				if ORF_start < ORF_end:
					readcov = int(subprocess.run(["samtools", "view", "-c", bamdb, "chr" + chrom + ":" + str(int(start) + ORF_start) + "-" + str(int(start) + ORF_end)], stdout=subprocess.PIPE).stdout.decode('utf-8')) / ORF_len / 6
					ORF_Info.write(str(readcov) + "," + celltype + ",")
					igvbat.write("goto chr" + chrom + ":" + start + "-" + end + "\n")
					igvbat.write("region chr" + chrom + " " + str(int(start) + ORF_start) + " " + str(int(start) + ORF_end) + "\ncollapse\n")
					igvbat.write("snapshot " + ORFid + "-broad-region.png\n")
					igvbat.write("goto chr" + chrom + ":" + str(int(start) + ORF_start) + "-" + str(int(start) + ORF_end) + "\ncollapse\n")
					igvbat.write("snapshot " + ORFid + "-ORF.png\n")
					r = get_overlapping_features("mus_musculus", chrom, str(int(start) + ORF_start), str(int(start) + ORF_end), str(strand)).json()
				else:
					readcov = int(subprocess.run(["samtools", "view", "-c", bamdb, "chr" + chrom + ":" + str(int(start) + ORF_end) + "-" + str(int(start) + ORF_start)], stdout=subprocess.PIPE).stdout.decode('utf-8')) / ORF_len / 6
					ORF_Info.write(str(readcov) + "," + celltype + ",")
					igvbat.write("goto chr" + chrom + ":" + start + "-" + end + "\n")
					igvbat.write("region chr" + chrom + " " + str(int(start) + ORF_end) + " " + str(int(start) + ORF_start) + "\ncollapse\n")
					igvbat.write("snapshot " + ORFid + "-broad-region.png\n")
					igvbat.write("goto chr" + chrom + ":" + str(int(start) + ORF_end) + "-" + str(int(start) + ORF_start) + "\ncollapse\n")
					igvbat.write("snapshot " + ORFid + "-ORF.png\n")
					r = get_overlapping_features("mus_musculus", chrom, str(int(start) + ORF_end), str(int(start) + ORF_start), str(strand)).json()
				ORF_Info.write(ORFid + "-broad-region.png," + ORFid + "-ORF.png,")

				#Identify genomic features
				if len(r) == 0:
					ORF_Info.write("intergenic")
				for feature in r:
					print(feature)
					classification = feature["id"] + "-" + feature["external_name"] + "-"
					if ORF_strand == feature["strand"]:
						classificationprefix = "sense"
					else:
						classificationprefix = "antisense"
					if ORF_start < ORF_end:
						r2 = get_overlapping_features("mus_musculus", chrom, str(int(start) + ORF_start), str(int(start) + ORF_end), str(strand),["exon"]).json()
						overlapclass = classify_overlap(int(feature["start"]), int(feature["end"]), int(start) + ORF_start, int(start) + ORF_end, ORF_strand)
					else:
						r2 = get_overlapping_features("mus_musculus", chrom, str(int(start) + ORF_end), str(int(start) + ORF_start), str(strand),["exon"]).json()
						overlapclass = classify_overlap(int(feature["start"]), int(feature["end"]), int(start) + ORF_end, int(start) + ORF_start, ORF_strand)
					classification += classificationprefix + "_" + overlapclass
					if overlapclass != "upstream" and overlapclass != "downstream":
						exonic = False
						for feature2 in r2:
							if ORF_start < ORF_end:
								exonoverlapclass = classify_overlap(int(feature2["start"]), int(feature2["end"]), int(start) + ORF_start, int(start) + ORF_end, ORF_strand)
							else:
								exonoverlapclass = classify_overlap(int(feature2["start"]), int(feature2["end"]), int(start) + ORF_end, int(start) + ORF_start, ORF_strand)
							if exonoverlapclass == "encompass":
								classification += "-exonic"
								exonic = True
							elif exonoverlapclass.find("overlap") != -1:
								classification += "-exon_overlap"
								exonic = True
							elif exonoverlapclass == "within":
								classification += "-encompass-exon"
								exonic = True
							if exonic:
								break
						if not exonic:
							classification += "-intronic"
					classification += ";"
					ORF_Info.write(classification)
				ORF_Info.write("\n")
		'''
		#Using Transeq
		alignments = subprocess.run(" ".join(["blastp", "-db", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".transeq.faa", "-query", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".fragments.faa", "-outfmt", "\"6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "nident\""]), stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8').split("\n")[:-1]
		exoninfo = subprocess.run(["grep", pattern, data_folder + "/reference_merged_assembled_transcripts.gtf"], stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[:-1]
		exonboundaries = []
		exonigvfile = open(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".igv.bat","w")
		exonigvfile.write("new\ngenome mm10\nload " + Tcellbam + "\nload " + Bcellbam + "\nsnapshotDirectory /media/hdd/Prabakaran_Lab/igv_snapshots\ngoto chr" + chrom + ":" + start + "-" + end + "\n")
		"""
		transcriptseq = ""
		for exon in exoninfo:
			exonterms = exon.split("\t")
			if exonterms[2] == "exon":
				print(exon)
				exonigvfile.write("region chr" + chrom + " " + str(exonterms[3]) + " " + str(exonterms[4]) + "\n")
				if exonterms[6] == "+":
					requesturl = server + "/sequence/region/mus_musculus/" + chrom + ":" + str(exonterms[3]) + ".." + str(exonterms[4]) + ":1?"
					r = requests.get(requesturl, headers={ "Content-Type" : "text/plain"})
					transcriptseq += r.text
				else:
					requesturl = server + "/sequence/region/mus_musculus/" + chrom + ":" + str(exonterms[3]) + ".." + str(exonterms[4]) + ":-1?"
					r = requests.get(requesturl, headers={ "Content-Type" : "text/plain"})
					transcriptseq += r.text
				exonboundaries.append((int(exonterms[3]),int(exonterms[4])))
		exonigvfile.close()
		transcriptfile = open(data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + "-transcript.fna","w")
		transcriptfile.write(">" + d + "\n" + transcriptseq)
		transcriptfile.close()
		"""
		translationfrombothstrands = ""
		forward = False
		backward = False
		for a in alignments:
			alignment_fields = a.split("\t")
			if alignment_fields[3] == alignment_fields[12] and int(alignment_fields[5]) == 0 and alignment_fields[12] == alignment_fields[13]:
				aaseq = "".join(subprocess.run(["pcregrep", "--buffer-size", "320000", "-M", ">[+|-]_" + alignment_fields[1][2] + "\\n([A-Z*]+\\n)+", data_folder + "/Potential_Novel_Peptide_ORFs/" + temp[0] + "." + temp[1].replace(".","_") + "." + temp[2] + ".transeq.faa"],stdout=subprocess.PIPE).stdout.decode('utf-8').split("\n")[1:-1])
				predictstart = int(alignment_fields[8]) - 1
				predictend = int(alignment_fields[9]) - 1
				print(len(aaseq))
				while predictstart > 0 and aaseq[predictstart] != "M" and aaseq[predictstart] != "*":
					predictstart -= 1
				while predictend < len(aaseq) and aaseq[predictend] != "*":
					predictend += 1
				if aaseq[predictstart] == "M":
					predictaa = aaseq[predictstart:predictend]
				else:
					predictaa = aaseq[predictstart+1:predictend]
				print(predictaa)
				fragnum = alignment_fields[0].split(":")[0]
				fragsource = ""
				if alignment_fields[0].split(":")[1].find("mascot") != -1:
					fragsource += "mascot"
				if alignment_fields[0].split(":")[1].find("perc") != -1:
					fragsource += "perc"
				if strand == 1:
					frameevidence.write(d + "," + fragsource + "," + alignment_fields[1][2] + "," + fragnum + "," + alignment_fields[10] + ",")
					if alignment_fields[1][2] == "4" or alignment_fields[1][2] == "5" or alignment_fields[1][2] == "6":
						alignmentstart = str(int(end) - int(alignment_fields[1][-1])-1 - int(alignment_fields[9])*3)
						alignmentend = str(int(end) - int(alignment_fields[1][-1])-1 - int(alignment_fields[8])*3)
						predictalignmentstart = str(int(end) - int(alignment_fields[1][-1])-1 - int(predictend)*3)
						predictalignmentend = str(int(end) - int(alignment_fields[1][-1])-1 - int(predictstart)*3)
						frameevidence.write(alignmentstart + "," + alignmentend + ",")
						frameevidence.write("antisense," + predictaa + "," + predictalignmentstart + "," + predictalignmentend + "\n")
						if float(alignment_fields[10]) < 0.01:
							backward = True
					else:
						alignmentstart = str(int(start) + int(alignment_fields[1][-1])-1 + int(alignment_fields[8])*3)
						alignmentend = str(int(start) + int(alignment_fields[1][-1])-1 + int(alignment_fields[9])*3)
						predictalignmentstart = str(int(start) + int(alignment_fields[1][-1])-1 + int(predictstart)*3)
						predictalignmentend = str(int(start) + int(alignment_fields[1][-1])-1 + int(predictend)*3)
						frameevidence.write(alignmentstart + "," + alignmentend + ",")
						if float(alignment_fields[10]) < 0.01:
							forward = True
						exonic = False
						for e in exonboundaries:
							if classify_overlap(int(alignmentstart),int(alignmentend),int(e[0]),int(e[1]),1) == "within":
								frameevidence.write("exonic,")
								exonic = True
								break
							elif classify_overlap(int(alignmentstart),int(alignmentend),int(e[0]),int(e[1]),1).find("overlap") != -1 or classify_overlap(int(alignmentstart),int(alignmentend),int(e[0]),int(e[1]),1) == "encompass":
								frameevidence.write("exon_overlap,")
								exonic = True
								break
						if not exonic:
							frameevidence.write("intronic,")
						frameevidence.write(predictaa + "," + predictalignmentstart + "," + predictalignmentend + "\n")
				else:
					alt = {"1" : "4", "2" : "5", "3" : "6", "4" : "1", "5" : "2", "6" : "3"}
					frameevidence.write(d + "," + fragsource + "," + alt[alignment_fields[1][2]] + "," + fragnum + "," + alignment_fields[10] + ",")
					if alignment_fields[1][2] == "4" or alignment_fields[1][2] == "5" or alignment_fields[1][2] == "6":
						alignmentstart = str(int(start) + int(alt[alignment_fields[1][-1]])-1 + int(alignment_fields[8])*3)
						alignmentend = str(int(start) + int(alt[alignment_fields[1][-1]])-1 + int(alignment_fields[9])*3)
						predictalignmentstart = str(int(start) + int(alt[alignment_fields[1][-1]])-1 + int(predictstart)*3)
						predictalignmentend = str(int(start) + int(alt[alignment_fields[1][-1]])-1 + int(predictend)*3)
						frameevidence.write(alignmentstart + "," + alignmentend + ",")
						frameevidence.write("antisense," + predictaa + "," + predictalignmentstart + "," + predictalignmentend + "\n")
						if float(alignment_fields[10]) < 0.01:
							backward = True
					else:
						alignmentstart = str(int(end) - int(alignment_fields[1][-1])-1 - int(alignment_fields[9])*3)
						alignmentend = str(int(end) - int(alignment_fields[1][-1])-1 - int(alignment_fields[8])*3)
						predictalignmentstart = str(int(end) - int(alignment_fields[1][-1])-1 - int(predictend)*3)
						predictalignmentend = str(int(end) - int(alignment_fields[1][-1])-1 - int(predictstart)*3)
						frameevidence.write(alignmentstart + "," + alignmentend + ",")
						if float(alignment_fields[10]) < 0.01:
							forward = True
						exonic = False
						for e in exonboundaries:
							if classify_overlap(int(alignmentstart),int(alignmentend),int(e[0]),int(e[1]),1) == "within":
								frameevidence.write("exonic,")
								exonic = True
								break
							elif classify_overlap(int(alignmentstart),int(alignmentend),int(e[0]),int(e[1]),1).find("overlap") != -1 or classify_overlap(int(alignmentstart),int(alignmentend),int(e[0]),int(e[1]),1) == "encompass":
								frameevidence.write("exon_overlap,")
								exonic = True
								break
						if not exonic:
							frameevidence.write("intronic,")
						frameevidence.write(predictaa + "," + predictalignmentstart + "," + predictalignmentend + "\n")
			if forward and backward:
				translationfrombothstrands = "bidirectional"
			else:
				translationfrombothstrands = "single_strand"

		"""
		r = get_overlapping_features("mus_musculus", chrom, str(int(start)-2000), str(int(end)+2000), str(strand)).json()
		for feature in r:
			print(feature)
			print(d)
			print(feature["gene_id"])
			print(feature["feature_type"])
			print(feature["id"])
			print(feature["external_name"])
			print(feature["start"])
			print(feature["end"])
			print(feature["biotype"])
			print(start)
			print(end)
			if feature['strand'] == strand:
				if feature['strand'] == 1:
					print(",sense_" + classify_overlap(int(start), int(end), feature['start'], feature['end'], 1) + "_biotype:" + feature["biotype"] + "-" + translationfrombothstrands)
					genomicfeatures.write(d + ",0," + feature["gene_id"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",sense," + classify_overlap(int(start), int(end), feature['start'], feature['end'], 1) + "," + feature["biotype"] + "," + translationfrombothstrands + "," + str(feature['start']) + "," + str(feature['end']) + "\n")
				else:
					print(",sense_" + classify_overlap(int(start), int(end), feature['start'], feature['end'], -1) + "_biotype:" + feature["biotype"] + "-" + translationfrombothstrands)
					genomicfeatures.write(d + ",0," + feature["gene_id"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",sense," + classify_overlap(int(start), int(end), feature['start'], feature['end'], -1) + "," + feature["biotype"] + "," + translationfrombothstrands + "," + str(feature['start']) + "," + str(feature['end']) + "\n")
			else:
				if feature['strand'] == 1:
					print(",anti-sense_" + classify_overlap(int(start), int(end), feature['start'], feature['end'], 1) + "_biotype:" + feature["biotype"] + "-" + translationfrombothstrands)
					genomicfeatures.write(d + ",0," + feature["gene_id"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",anti-sense," + classify_overlap(int(start), int(end), feature['start'], feature['end'], 1) + "," + feature["biotype"] + "," + translationfrombothstrands + "," + str(feature['start']) + "," + str(feature['end']) + "\n")
				else:
					print(",anti-sense_" + classify_overlap(int(start), int(end), feature['start'], feature['end'], -1) + "_biotype:" + feature["biotype"] + "-" + translationfrombothstrands)
					genomicfeatures.write(d + ",0," + feature["gene_id"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",anti-sense," + classify_overlap(int(start), int(end), feature['start'], feature['end'], -1) + "," + feature["biotype"] + "," + translationfrombothstrands + "," + str(feature['start']) + "," + str(feature['end']) + "\n")
		"""

		
				
		'''for a in alignments:
			print(a)
			temp3 = a.split("\t")
			frameevidence.write(d + "," + temp3[1][-1] + "," + temp3[10] + ",")
			if int(temp3[1][-1]) <= 3:
				frameevidence.write(str(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3) + "," + str(int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3) + "\n")
				r = get_overlapping_features("mus_musculus", chrom, str(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3 - 500), str(int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3 + 500), "1", ["transcript","cds","exon","regulatory"]).json()
				for feature in r:
					print(feature)
					if feature["feature_type"] == "transcript":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",anti-sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						prefix = "transcript_"
					elif feature["feature_type"] == "cds":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["protein_id"] + ",none,sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["protein_id"] + ",none,anti-sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						prefix = "protein_"
					elif feature["feature_type"] == "exon":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + ",none,sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + ",none,anti-sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
					elif feature["feature_type"] == "regulatory":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + ",none," + feature["feature_type"] + "," + feature["id"] + "," + feature["description"] + ",sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + ",none," + feature["feature_type"] + "," + feature["id"] + "," + feature["description"] + ",anti-sense_" + classify_overlap(int(start) + int(temp3[1][-1])-1 + int(temp3[8])*3, int(start) + int(temp3[1][-1])-1 + int(temp3[9])*3, feature['start'], feature['end'], "+") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")	
					
			else:
				offset = region_length % 3 + int(temp3[1][-1]) - 4
				frameevidence.write(str(int(end) - offset - int(temp3[9])*3) + "," + str(int(end) - offset - int(temp3[8])*3) + "\n")
				r = get_overlapping_features("mus_musculus", chrom, str(int(end) - offset - int(temp3[9])*3 - 500), str(int(end) - offset - int(temp3[8])*3 + 500), "-1", ["transcript","cds","exon","regulatory"]).json()
				for feature in r:
					print(feature)
					if feature["feature_type"] == "transcript":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + "," + feature["external_name"] + ",anti-sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
					elif feature["feature_type"] == "cds":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["protein_id"] + ",none,sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["protein_id"] + ",none,anti-sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
					elif feature["feature_type"] == "exon":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + ",none,sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + "," + feature["Parent"] + "," + feature["feature_type"] + "," + feature["id"] + ",none,anti-sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
					elif feature["feature_type"] == "regulatory":
						if feature["strand"] == strand:
							genomicfeatures.write(d + "," + temp3[0] + ",none," + feature["feature_type"] + "," + feature["id"] + "," + feature["description"] + ",sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
						else:
							genomicfeatures.write(d + "," + temp3[0] + ",none," + feature["feature_type"] + "," + feature["id"] + "," + feature["description"] + ",anti-sense_" + classify_overlap(int(end) - offset - int(temp3[9])*3, int(end) - offset - int(temp3[8])*3, feature['start'], feature['end'], "-") + "," + str(feature['start']) + "," + str(feature['end']) + "," + temp3[10] + "\n")
					
		
		#for i in range(1,counter + 1):
		for a in alignments:
			temp = a.split("\t")
			alignedseq = subprocess.run(["pcregrep", "--buffer-size", "320000", "-M", ">\\" + temp[1] + "\\n([A-Z\\*]+\\n)+"], stdout=subprocess.PIPE).stdout.decode('utf-8')
			alignedseq = alignedseq.split("\n")[1:-1]'''
	ORF_Info.close()
	ORF_Fragment_Info.close()
	igvbat.close()
	genomicfeatures.close()
	frameevidence.close()
