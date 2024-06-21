import argparse
import pysam
import os
import subprocess
import re
import pandas

from Bio.Seq import Seq


yeast_genome="/media/theboocock/Data/Dropbox/Postdoc/projects/crispr_coupling/data/tn5_96/saccer3_plk88.fasta" 
ALIGNMENT_BUFFER=20
class Region(object):

    def __init__(self, chrom, start, end):
        self._chrom = chrom
        self._start = start
        self._end = end 
    def start(self):
        return self._start
    def end(self):
        return self._end
    def chrom(self):
        return self._chrom


def filter_bam_and_run_spades(bam,out_prefix, oligo_table,log_f):
    out_bam =  out_prefix + "_plk88.bam"
    samfile = pysam.AlignmentFile(bam,"rb")
    outfile = pysam.AlignmentFile(out_bam,"wb",template=samfile)
    # TODO: Read pairs 
    for read in samfile.fetch("plk88"):
        if(read.mapping_quality > 20):
            outfile.write(read)
    ### Run spades ### 
    outfile.close()
    out_path = (outfile.filename.decode())
    out_fastq = out_prefix + "_plk88.fastq"
    cmd = f"samtools bam2fq {out_path} > {out_fastq}"
    subprocess.check_call(cmd, shell=True) 
    output_spades = "spades"
    cmd_spades = f"spades.py -s {out_fastq} -o {output_spades}"
    # Final GFA #
    subprocess.check_call(cmd_spades, shell=True)
    out_gfa = f"{output_spades}/assembly_graph_after_simplification.gfa"
    matched_strings = {}
    lengths = []
    i = 0
    max_len = 0
    with open(out_gfa) as in_gfa:
        for line in in_gfa:
            l_s = line.split("\t")
            type_str = l_s[0]
            if type_str == "S":
                seq= l_s[2]
                matched_strings[str(i)]=seq
                lengths.append(len(seq))
                if(len(seq) > max_len):
                    max_idx = i
                i = i + 1
    if(len(matched_strings) == 0):
        # No assembly # 
        log_f.write("Failed at the assembly step\n")
        return False
    a = (matched_strings[str(max_idx)])
    replacement = "$1_$2_$3"
    str_match_structural_region=r"GCAGTGAAAGATAGGTGACC(.{20})(.*)(.{90})TTCCCGACGAGAGTAAATGGCGAGGATACGTTCTCTATGG(.{30})"
    str_rep = re.search(str_match_structural_region,a)
    # B02 BAM
    a_rev=Seq(a).reverse_complement()
    a_rev=(str(a_rev))
    str_rep_rev = re.search(str_match_structural_region,a_rev)
    if(str_rep is None):
        if(str_rep_rev is not None):
            str_rep = str_rep_rev
    out_match = out_prefix + "_plasmid.txt"
    if (str_rep is not None):
        #print(str_rep)
        with open(out_match, 'w') as out_f:
            out_f.write("guide\tstruct\trepair\tbarcode\n")
            if(str_rep is not None):
                a1 = str_rep.groups()
                guide = a1[0]
                struct = a1[1]
                repair = a1[2]
                barcode=a1[3]
                out_f.write(guide + "\t" + struct + "\t" + repair + "\t" + barcode+ "\n")
    else:
        log_f.write("Failed to find structural region\n")
        return False 
        # Don't run the Rscript if you don't have 
    ### RSCRIPT to check info ###
    # Oligo table 
    #print("Run rscript")
    dir_path = os.path.dirname(os.path.realpath(__file__))
    out_rscript = os.path.join(dir_path, "check_match.R")
    rscript_cmd = f"/usr/bin/Rscript {out_rscript} {out_match} {oligo_table} {out_prefix}" 
    subprocess.check_call(rscript_cmd, shell=True)
    out_rscript = out_prefix + "_combo_match.txt"
    regions_to_remove = list()
    with open(out_rscript) as in_f:
        line = in_f.readline().strip()
        if(line == "no_match"):
            log_f.write("Incorrect oligo combination")
            return False 
        else:
            # Has a match # 
            l_s = line.split()
            chrom = l_s[0]
            #print("Line 83")
            #print(chrom)
            out_fasta = out_prefix + "_repair.fasta"
            with open(out_fasta, "w") as out_f:
               out_f.write(">rep\n")
               out_f.write(repair + "\n")

            makeblastdb_cmd=f"makeblastdb -in {yeast_genome} -parse_seqids -dbtype nucl"
            #
            #subprocess.check_call(makeblastdb_cmd, shell=True)
            out_blast=out_prefix  + "_blast.txt"
            blastn_cmd=f"blastn -query {out_fasta} -db {yeast_genome} -outfmt 6 -task blastn > {out_blast}"
            subprocess.check_call(blastn_cmd, shell=True)
            got_match=False
            with open(out_blast) as blast_in:
                for line in blast_in:
                    line_s = line.strip().split("\t")
                    chrom_blast = line_s[1]
                    if chrom_blast == chrom:
                        got_match = True
                    length_match = int(line_s[3])
                    if (length_match < 80):
                        continue
                    start = int(line_s[8])
                    end = int(line_s[9])
                    if(start > end):
                        tmp = start
                        end = start
                        start = end
                    region = chrom_blast +":" + str(start) + "-" + str(end)
                    start_tmp = start - ALIGNMENT_BUFFER
                    if(start_tmp < 0):
                        start_tmp = 0

                    reg_tmp = Region(chrom_blast, start-ALIGNMENT_BUFFER, end + ALIGNMENT_BUFFER)
                    regions_to_remove.append(reg_tmp)
    return regions_to_remove

MAPQ_THRESHOLD=20
def filter_pairs_one_mapping_to_plasmid(bam, out_prefix,regions, out_f):
    """
        
    """
    out_bam = out_prefix+"_bam_filtered_scuffed.bam"
    #bam_sort_qname = os.path.join(out_prefix,"qname_sort.bam") 
    #samtools_sort_by_name = f"samtools sort -n {bam} > {bam_sort_qname}"
    #subprocess.check_call(samtools_sort_by_name, shell=True)
    samfile= pysam.AlignmentFile(bam,"rb")
    outfile = pysam.AlignmentFile(out_bam,"wb",template=samfile)
    bamiter = samfile.fetch(until_eof=True)
    if(len(regions) > 0):
        read_list = {}
        i=1
        qname_filt = list()
        for read in bamiter:
            qname = read.qname
            read1 = read.is_read1
            try:
                XAs = read.get_tag("XA")
                chrom_xa = [x.split(",")[0] ==  "plk88" for x in XAs.split(";")]
                if(any(chrom_xa)):
                    if(qname not in read_list):
                        read_list[qname]={}
                    read_list[qname]["XA_plk88"] = True 
            except:
                pass
            if read1:
                try:
                    read_list[qname]["r1"].append(read)
                except:
                    if(qname not in read_list):
                        read_list[qname]={}
                    read_list[qname]["r1"] = [read]
            else: 
                try:
                    read_list[qname]["r2"].append(read)
                except:
                    if(qname not in read_list):
                        read_list[qname]={}
                    read_list[qname]["r2"] = [read]
                i=i+1

        #print(read_list["VH00763:48:AAFKNYMM5:1:1402:55667:31196"])
        
        for qname in read_list:
            read_l_tmp = read_list[qname]
            try:
                read_ones = read_l_tmp["r1"]
            except:
                next
            reference_namesr1=([read.reference_name == "plk88" for read in read_ones])
            try:
                read_twos = read_l_tmp["r2"]
            except:
                next 
            reference_namesr2=([read.reference_name == "plk88" for read in read_twos])
            try:
                XA_plk88 = read_l_tmp["XA_plk88"]
                continue
            except:
                pass
            if(not any(reference_namesr1) and not any(reference_namesr2)):
                #print("HERE")
                chrom1 = [read1.reference_name for read1 in read_ones]
                chrom2 = [read2.reference_name for read2 in read_twos]
                start1 = [read1.reference_start for read1 in read_ones]
                start2 = [read2.reference_start for read2 in read_twos]
                end1 = [read1.reference_end for read1 in read_ones]
                end2 = [read2.reference_end for read2 in read_twos]
                for region in regions:
                    r1_match = False
                    r2_match = False
                    chrom =region.chrom()
                    left = region.start() 
                    right = region.end()
                    i= 0
                    for chromr1 in chrom1:
                        if chrom == chromr1:
                            try:
                                if start1[i] >= left and end1[i] <= right: 
                                    r1_match = True
                            except:
                                r1_match = True
                        i = i + 1
                    i = 0
                    for chromr2 in chrom2: 
                        if chrom == chromr2:
                            try:
                                if start2[i] >=left and end2i[i] <= right:
                                    r2_match = True
                            except:
                                r2_match = True
                        i = i + 1 
                    if r1_match and r2_match:
                        #print("HERE")
                        break
                if not (r1_match and r2_match):
                    for read in read_ones + read_twos:
                        if read.mapping_quality > MAPQ_THRESHOLD:
                            outfile.write(read)
        outfile.close()
        final_bam = out_prefix + "_final.bam"
        pysam.sort("-o",final_bam, out_bam)
        pysam.index(final_bam)
    else:
        final_bam = out_prefix + "_final.bam"
        pysam.sort("-o", final_bam, bam)
        pysam.index(final_bam)
    out_f.write("Completed pipeline with no errors \n") 

def main():
    parser = argparse.ArgumentParser(description="Filters reads for spades")
    parser.add_argument("-o","--out-prefix", help="Output prefix",default="test/")
    parser.add_argument("--oligo-table",default="/media/theboocock/Data/Dropbox/Postdoc/projects/crispr_coupling/data/tn5_96/otableLong.RDS")
    parser.add_argument(dest="bam", help="Input bam")
    args = parser.parse_args()
    with open(args.out_prefix +"_stats.txt","w") as out_f:
        regions = filter_bam_and_run_spades(args.bam, args.out_prefix, args.oligo_table,out_f)
        # Returns false if it di
        if regions != False:
        # Regions might be zero when randomer is the location # 
            filter_pairs_one_mapping_to_plasmid(args.bam, args.out_prefix,regions,out_f)


if __name__=="__main__":
    main()
