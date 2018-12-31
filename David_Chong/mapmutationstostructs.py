'''
Mutations from HGMD and Cosmic are mapped to sorfs based on converted genomic coordinates using Liftover and tblastn.
Mapped mutations that fall on a sorf coding sequence that has a predicted structure are tested for effect based on standard amino acid code.
Affected site on the pdb structure is highlighted in red, accounting for dropped residues during structure prediction
'''
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Alphabet import generic_dna
from pymol import *
import sys, subprocess, os, requests
import Bio.PDB.Polypeptide as pp

def getseq(starts,ends,chrom,species,strand="1"):
     transcriptseq = ""
     for i in range(0,len(starts)):
             requesturl = server + "/sequence/region/" + species + "/" + chrom + ":" + str(starts[i]) + ".." + str(ends[i]) + ":" + strand + "?"
             r = requests.get(requesturl, headers={ "Content-Type" : "text/plain"})
             transcriptseq+=r.text
     return(transcriptseq)

server = "https://rest.ensembl.org"
data_folder = "/media/hdd/Prabakaran_Lab/data"
cosmic_coding = data_folder + "/CosmicCodingVariants_Mapped_to_sorfs_identified_in_B_Tcells_info.csv"
hgmd = data_folder + "/HGMDVariants_Mapped_to_sorfs_identified_in_B_Tcells_info.csv"
cosmic_non_coding = data_folder + "/CosmicNonCodingVariants_Mapped_to_sorfs_identified_in_B_Tcells_info.csv"
mapped_sorfs = data_folder + "/sorfs_identified_in_B_Tcells_info_hg38_OverlapedsOftblastnAndliftOver.csv"

pdbdir = "/media/hdd/Prabakaran_Lab/data/sORFswithCmutations"

if len(sys.argv) > 1:
	pdbfiles = {}
	for p in os.listdir(sys.argv[1]):
		if p.find(".pdb") != -1:
			temp = p.split("_")
			pdbfiles[temp[0]] = {}
			pdbfiles[temp[0]]["pdb_file"] = p
			f = open(sys.argv[1] + "/" + p)
			data = f.read().split("\n")[:-1]
			f.close()
			pdbseq = ""
			pdbstart = 0
			pdbend = 0
			curpos = 0
			for line in data:
				if line[0:4] == "ATOM":
					pdb_fields = list(filter(None,line.split(" ")))
					if curpos == 0:
						curpos = pdb_fields[5]
						pdbstart = pdb_fields[5]
						pdbseq += pp.three_to_one(pdb_fields[3])
					elif pdb_fields[5] != curpos:
						curpos = pdb_fields[5]
						pdbseq += pp.three_to_one(pdb_fields[3])					
			pdbend = curpos
			pdbfiles[temp[0]]["pdb_seq"] = pdbseq
			pdbfiles[temp[0]]["pdb_start"] = pdbstart
			pdbfiles[temp[0]]["pdb_end"] = pdbend
			alignments = subprocess.Popen(["grep", temp[0], mapped_sorfs], stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split("\n")[:-1]
			for a in alignments:
				alignment_fields = a.split(",")
				chrom = alignment_fields[6]
				pdbfiles[temp[0]]["chrom"] = chrom
				nucstart = alignment_fields[16]
				nucend = alignment_fields[17]
				aastart = alignment_fields[14]
				aaend = alignment_fields[15]
				aligned_aaseq = alignment_fields[18]
				if "aa_start" not in pdbfiles[temp[0]].keys():
					pdbfiles[temp[0]]["aa_start"] = int(aastart)
					pdbfiles[temp[0]]["aa_end"] = int(aaend)
					pdbfiles[temp[0]]["genome_start"] = int(nucstart)
					pdbfiles[temp[0]]["genome_end"] = int(nucend)
					pdbfiles[temp[0]]["aligned_aaseq"] = aligned_aaseq
					pdbfiles[temp[0]]["chrom"] = chrom
				else:
					if int(aastart) < pdbfiles[temp[0]]["aa_start"] and int(aaend) > pdbfiles[temp[0]]["aa_end"]:
						pdbfiles[temp[0]]["aa_start"] = int(aastart)
						pdbfiles[temp[0]]["aa_end"] = int(aaend)
						pdbfiles[temp[0]]["genome_start"] = int(nucstart)
						pdbfiles[temp[0]]["genome_end"] = int(nucend)
						pdbfiles[temp[0]]["aligned_aaseq"] = aligned_aaseq
						pdbfiles[temp[0]]["chrom"] = chrom
				
			pdbfiles[temp[0]]["aligned_nucseq"] = getseq([pdbfiles[temp[0]]["genome_start"]], [pdbfiles[temp[0]]["genome_end"]], pdbfiles[temp[0]]["chrom"], "human")

			mutations = subprocess.Popen(["grep", temp[0], cosmic_coding], stdout=subprocess.PIPE)
			mutationids = subprocess.Popen(["cut", "-d", ",", "-f", "12"], stdin=mutations.stdout, stdout=subprocess.PIPE)
			sortedmutationids = subprocess.Popen(["sort"], stdin=mutationids.stdout, stdout=subprocess.PIPE)
			uniquemutationids = subprocess.Popen(["uniq"], stdin=sortedmutationids.stdout, stdout=subprocess.PIPE)
			pdbfiles[temp[0]]["cosmic_coding_mutations"] = uniquemutationids.communicate()[0].decode('utf-8').split("\n")[:-1]

			mutations = subprocess.Popen(["grep", temp[0], cosmic_non_coding], stdout=subprocess.PIPE)
			mutationids = subprocess.Popen(["cut", "-d", ",", "-f", "4"], stdin=mutations.stdout, stdout=subprocess.PIPE)
			sortedmutationids = subprocess.Popen(["sort"], stdin=mutationids.stdout, stdout=subprocess.PIPE)
			uniquemutationids = subprocess.Popen(["uniq"], stdin=sortedmutationids.stdout, stdout=subprocess.PIPE)
			pdbfiles[temp[0]]["cosmic_non_coding_mutations"] = uniquemutationids.communicate()[0].decode('utf-8').split("\n")[:-1]

			mutations = subprocess.Popen(["grep", temp[0], hgmd], stdout=subprocess.PIPE)
			mutationids = subprocess.Popen(["cut", "-d", ",", "-f", "13"], stdin=mutations.stdout, stdout=subprocess.PIPE)
			sortedmutationids = subprocess.Popen(["sort"], stdin=mutationids.stdout, stdout=subprocess.PIPE)
			uniquemutationids = subprocess.Popen(["uniq"], stdin=sortedmutationids.stdout, stdout=subprocess.PIPE)
			pdbfiles[temp[0]]["hgmd_mutations"] = uniquemutationids.communicate()[0].decode('utf-8').split("\n")[:-1]

			mutations = []
			for mut in pdbfiles[temp[0]]["cosmic_coding_mutations"]:
				info = subprocess.Popen(["grep", mut, cosmic_coding], stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split("\n")[0].split(",")
				mut_start = info[1]
				mut_end = info[2]
				mut_chrom = info[0]
				original_base = info[12]
				mutated_base = info[13]
				effect = info[17]
				existing_mutation = False
				for m in mutations:
					if m[0] == mut_start and m[1] == mut_end and m[2] == mut_chrom and m[3] == original_base and m[4] == mutated_base:
						m[5].append(mut + "-" + effect)
						existing_mutation = True
						break
				if not existing_mutation:
					mutations.append([mut_start,mut_end,mut_chrom,original_base,mutated_base,[mut + "-" + effect]])
			for mut in pdbfiles[temp[0]]["cosmic_non_coding_mutations"]:
				info = subprocess.Popen(["grep", mut, cosmic_non_coding], stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split("\n")[0].split(",")
				mut_start = info[1]
				mut_end = info[2]
				mut_chrom = info[0]
				original_base = info[4]
				mutated_base = info[5]
				existing_mutation = False
				for m in mutations:
					if m[0] == mut_start and m[1] == mut_end and m[2] == mut_chrom and m[3] == original_base and m[4] == mutated_base:
						m[5].append(mut)
						existing_mutation = True
						break
				if not existing_mutation:
					mutations.append([mut_start,mut_end,mut_chrom,original_base,mutated_base,[mut]])

			for mut in pdbfiles[temp[0]]["hgmd_mutations"]:
				info = subprocess.Popen(["grep", mut, hgmd], stdout=subprocess.PIPE).communicate()[0].decode('utf-8').split("\n")[0].split(",")
				mut_start = info[1]
				mut_end = info[2]
				mut_chrom = info[0]
				original_base = info[13]
				mutated_base = info[14]
				effect = info[20]
				existing_mutation = False
				for m in mutations:
					if m[0] == mut_start and m[1] == mut_end and m[2] == mut_chrom and m[3] == original_base and m[4] == mutated_base:
						m[5].append(mut + "-" + effect)
						existing_mutation = True
						break
				if not existing_mutation:
					mutations.append([mut_start,mut_end,mut_chrom,original_base,mutated_base,[mut + "-" + effect]])
				
			pdbfiles[temp[0]]["all_mutations"] = mutations
	f = open(sys.argv[1] + ".html","w")
	f.write("<html><body><h2>Mutation Information for " + sys.argv[1] + "</h2>")
	for k,v in pdbfiles.items():
		print(k)
		print(v["aligned_nucseq"])
		print(v["pdb_seq"])
		print(v["pdb_start"])
		f2 = open(sys.argv[1] + "/" + k + "/" + k + "/fold/psipred/" + k + ".ss")
		sinfo = f2.read().split("\n")[:-1]
		f2.close()

		cmd.load(sys.argv[1] + "/" + v["pdb_file"],"cartoon")
		sorfstruct = cmd.get_object_list()[0]
                cmd.hide("all")
		cmd.show("cartoon",sorfstruct)
		cmd.color("cyan",sorfstruct)

		sstruct = []
		mutationlocs = []
		f.write("<table border='1'><tr><th colspan='5'>" + k + "<br /><img src='" + sys.argv[1] + "/" + k + "_mutations.png'></th></tr>")
		f.write("<tr><th>Mutation ID(s)</th><th>Change</th><th>Effect on sORF</th><th>sORF secondary structure</th></tr>")
		for m in v["all_mutations"]:
			if len(m[3]) == 1 and len(m[4]) == 1:
				affected_coord = (int(m[1]) - int(v["genome_start"]) + (int(v["aa_start"]) - 1) * 3)
				seq_coord = affected_coord - (int(v["aa_start"]) - 1) * 3
				if affected_coord < len(v["aligned_nucseq"]) and affected_coord >= 0:
					if int(affected_coord/3) > int(v["pdb_start"]) and int(affected_coord/3) < int(v["pdb_end"]):
						cmd.select("highlight","resi " + str(int(affected_coord/3)-int(v["pdb_start"])) + " in " + sorfstruct)
                    				cmd.color("red","highlight")
					if affected_coord % 3 == 0:
						seq = v["aligned_nucseq"][seq_coord:seq_coord+3]
						mutatedseq = m[4] + seq[1:]
					elif affected_coord % 3 == 1:
						seq = v["aligned_nucseq"][seq_coord-1:seq_coord+2]
						mutatedseq = seq[0] + m[4] + seq[2]
					else:
						seq = v["aligned_nucseq"][seq_coord-2:seq_coord+1]
						mutatedseq = seq[:2] + m[4]
					originalaa = Seq(seq, generic_dna).translate()
					mutatedaa = Seq(mutatedseq, generic_dna).translate()
					if int(affected_coord/3) < len(sinfo):
						if int(affected_coord/3) > int(v["pdb_start"]) and int(affected_coord/3) < int(v["pdb_end"]):
							f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>" + str(originalaa) + str(int(affected_coord/3)) + ">" + str(mutatedaa) + " (PDB " + v["pdb_seq"][int(affected_coord/3)-int(v["pdb_start"]) + 1] + str(int(affected_coord/3) - int(v["pdb_start"]) + 2) + ")</td><td>" + list(filter(None, sinfo[int(affected_coord/3)].split(" ")))[2] + "</td></tr>")
						else:
							f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>" + str(originalaa) + str(int(affected_coord/3)) + ">" + str(mutatedaa) + "</td><td>" + list(filter(None, sinfo[int(affected_coord/3)].split(" ")))[2] + "</td></tr>")
					else:
						if int(affected_coord/3) > int(v["pdb_start"]) and int(affected_coord/3) < int(v["pdb_end"]):
							f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>" + str(originalaa) + str(int(affected_coord/3)) + ">" + str(mutatedaa) + " (PDB " + v["pdb_seq"][int(affected_coord/3)-int(v["pdb_start"]) + 1] + str(int(affected_coord/3) - int(v["pdb_start"]) + 2) + ")</td><td>Out of sORF seq</td></tr>")
						else:
							f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>" + str(originalaa) + str(int(affected_coord/3)) + ">" + str(mutatedaa) + "</td><td>Out of sORF seq</td></tr>")
				else:
					f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td></td><td>Out of nucleotide seq, potential gapped alignment</td></tr>")
			else:
				change = int(m[1]) - int(m[0])
				affected_start = (int(m[0]) - int(v["genome_start"]) + (int(v["aa_start"]) - 1) * 3)/3
				affected_end = (int(m[1]) - int(v["genome_start"]) + (int(v["aa_start"]) - 1) * 3)/3
				if int(affected_start) > int(v["pdb_start"]) and int(affected_start) < int(v["pdb_end"]):
					cmd.select("highlight","resi " + str(int(affected_start) - int(v["pdb_start"])) + " in " + sorfstruct)
                    			cmd.color("red","highlight")
				if change % 3 == 0 and len(m[4]) - len(m[3]) < 0:
					f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>In-frame deletion at " + str(int(affected_start)) + "</td><td></td></tr>")
				elif change % 3 == 0 and len(m[4]) - len(m[3]) > 0:
					f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>In-frame insertion at " + str(int(affected_start)) + "</td><td></td></tr>")
				elif change % 3 != 0 and len(m[4]) - len(m[3]) < 0:
					f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>Frameshift deletion at " + str(int(affected_start)) + "</td><td></td></tr>")
				elif change % 3 != 0 and len(m[4]) - len(m[3]) > 0:
					f.write("<tr><td>" + "<br />".join(m[5]) + "</td><td>" + m[3] + ">" + m[4] + "</td><td>Frameshift insertion at " + str(int(affected_start)) + "</td><td></td></tr>")
		f.write("</table>")
		cmd.orient(sorfstruct)
		cmd.png(sys.argv[1] + "/" + k + "_mutations.png",dpi=300)
		cmd.delete(sorfstruct)
	f.close()
