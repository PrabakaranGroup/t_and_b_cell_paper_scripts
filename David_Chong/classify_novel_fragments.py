#!/usr/bin/env python
import sys, subprocess, os, requests, time

server = "https://rest.ensembl.org"
data_folder = "/media/hdd/Prabakaran_Lab/data"

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

f = open(data_folder + "/novel_transcript_fragment_alignment_list_filtered.txt")
data = f.read().split("\n")[:-1]
f.close()


feature_info = {}
output = open(data_folder + "/novel_fragmennt_alignment_classification.txt","w")
output.write("region\tgene_id\tgene_start\tgene_end\tstrand_relationship\tcoordinate_relationship\tbiotype\ttranscript_relationship\texonic_relationship\tcds_relationship\n")

for line in data:
	fragment_fields = line.split(":")
	gene = fragment_fields[0]
	chrom = fragment_fields[1][3:]
	start = fragment_fields[2]
	end = fragment_fields[3]
	strand = fragment_fields[4]
	feature_info[line] = []
	features = get_overlapping_features("mus_musculus",chrom, str(int(start)-500), str(int(end)+500)).json()
	if len(features) > 0:
		genelist = {}
		for f in features:
			if str(f["strand"]) == strand:
				classification = classify_overlap(int(start),int(end),f["start"],f["end"],f["strand"])
				genebiotype = f["biotype"]
				geneid = f["id"]
				genelist[geneid] = []
				transcriptids = {}
				transcriptinfo = get_overlapping_features("mus_musculus",chrom, str(int(start)-500), str(int(end)+500),featurelist=["transcript"]).json()
				#time.sleep(0.1)
				for t in transcriptinfo:
					if t["Parent"] == geneid:
						genelist[geneid] = t["id"]
						transcriptids[t["id"]] = [classify_overlap(int(start),int(end),t["start"],t["end"],t["strand"]),t["start"],t["end"]]
				exoninfo = get_overlapping_features("mus_musculus",chrom, str(int(start)-500), str(int(end)+500),featurelist=["exon"]).json()
				exon_relation = {}
				for e in exoninfo:
					if e["Parent"] in transcriptids.keys():
						if transcriptids[e["Parent"]][0] == "within":
							test = classify_overlap(int(start),int(end),e["start"],e["end"],e["strand"])
							if test != "upstream" and test != "downstream":
								exon_relation[e["Parent"]] = test
							else:
								exon_relation[e["Parent"]] = "intronic"
				cds_relation = {}
				for t, tv in transcriptids.items():
					cdsinfo = get_overlapping_features("mus_musculus",chrom, str(tv[1]), str(tv[2]),featurelist=["cds"]).json()
					collatecds = []
					mincds = -1
					maxcds = -1
					for c in cdsinfo:
						if c["Parent"] == t:
							if tv[0] == "within":
								if mincds == -1 or mincds > c["start"]:
									mincds = c["start"]
								if maxcds == -1 or maxcds < c["end"]:
									maxcds = c["end"]
								collatecds.append([c["start"],c["end"]])								
					if mincds != -1 or maxcds != -1 or len(collatecds) != 0:
						if int(start) < mincds and str(strand) == "1":
							cds_relation[t] = "5UTR"
						elif int(end) > maxcds and str(strand) == "-1":
							cds_relation[t] = "5UTR"
						elif int(start) < mincds and str(strand) == "-1":
							cds_relation[t] = "3UTR"
						elif int(end) > maxcds and str(strand) == "1":
							cds_relation[t] = "3UTR"
						else:
							for cf in collatecds:
								test = classify_overlap(int(start),int(end),cf[0],cf[1],int(strand))
								if test == "within" or test.find("overlap") != -1 or test.find("truncation") != -1 or test.find("extension") != -1:
									cds_relation[t] = test
				'''for g, gv in genelist:
					tag = output.write(line + "\t" + f["id"] + "\t" + str(f["start"]) + "\t" + str(f["end"]) + "\tsense\t" + classification + "\t" + genebiotype + "\t")
					transcript_classes = []
					cds_classes = []
					exon_classes = []
					for t in gv:
						transcript_classes.append(transcriptids[t][0])
						if t in exon_relation.keys():
							exon_classes.append(exon_relation[t])
						elif transcriptids[t][0] == "within" or transcriptids[t][0] == "upstream_truncation" or transcriptids[t][0]=="downstream_truncation":
							exon_classes.append("intronic")
						else:
							exon_classes.append("NA")
						if t in cds_relation.keys():
							cds_classes.append(cds_relation[t])
						else:
							cds_classes.append("NA")
					if "within" in transcript_classes or "upstream_truncation" in transcript_classes or "downstream_truncation" in transcript_classes:
						tag += "within\t"
					else:
						tag += "NA\t"
					if "exon'''
				for t,tv in transcriptids.items():
					tag = "strand:sense-coordinate:" + classification + "-biotype:" + genebiotype + "-transcriptrelation:" + tv[0]
					output.write(line + "\t" + f["id"] + "\t" + str(f["start"]) + "\t" + str(f["end"]) + "\tsense\t" + classification + "\t" + genebiotype + "\t" + tv[0] + "\t")
					if t in exon_relation.keys():
						output.write(exon_relation[t] + "\t")
						tag += "-exonrelation:" + exon_relation[t]
					elif tv[0] == "within" or tv == "upstream_truncation" or tv=="downstream_truncation":
						output.write("intronic\t")
					else:
						output.write("NA\t")
					if t in cds_relation.keys():
						output.write(cds_relation[t] + "\n")
						tag += "-cdsrelation:" + cds_relation[t]
					else:
						output.write("NA\n")
					feature_info[line].append(tag)
					print(tag)	
			else:
				feature_info[line] == "antisense-" + f["biotype"]
				output.write(line + "\t" + f["id"] + "\t" + str(f["start"]) + "\t" + str(f["end"]) + "\tanti-sense\tNA\t" + f["biotype"] + "\tNA\tNA\tNA\n")
	else:
		output.write(line + "\tNA\tNA\tNA\tintergenic\tNA\tNA\tNA\tNA\tNA\n")
		feature_info[line].append("intergenic")
output.close()
