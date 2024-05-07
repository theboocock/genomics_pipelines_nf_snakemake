import argparse
import pysam
import os
import subprocess
import re
import pandas

def filter_bam_and_run_spades(bam,out_prefix):
    out_bam = os.path.join(out_prefix,"plk88.bam")
    samfile = pysam.AlignmentFile(bam,"rb")
    print(out_bam)
    outfile = pysam.AlignmentFile(out_bam,"wb",template=samfile)
    # TODO: Read pairs 
    for read in samfile.fetch("plk88"):
        if(read.mapping_quality > 20):
            outfile.write(read)
    ### Run spades ### 
    outfile.close()
    out_path = (outfile.filename.decode())
    out_fastq = os.path.join(out_prefix,"plk88.fastq")
    cmd = f"samtools bam2fq {out_path} > {out_fastq}"
    subprocess.check_call(cmd, shell=True) 
    output_spades = "test/spades"
    cmd_spades = f"plasmidspades.py -s {out_fastq} -o {output_spades}"
    # Final GFA #
    #subprocess.check_call(cmd_spades, shell=True)
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
    a = (matched_strings[str(max_idx)])
    replacement = "$1_$2_$3"
    str_match_structural_region=r"GCAGTGAAAGATAGGTGACC(.{20})(.*)(.{30})GGTACCCAATTCGCC"
    str_rep = re.search(str_match_structural_region,a)
    if(str_rep is not None):
        a1 = re.search(str_match_structural_region,a).groups()
        guide = a1[0]
        struct = a1[1]
        repair = a1[2]
        print("guide: ", guide, "struct: ", struct, "repair:", repair)


def filter_pairs_one_mapping_to_plasmid(bam, out_prefix):
    

    out_bam = os.path.join(out_prefix,"bam_filtered.bam")
    samfile= pysam.AlignmentFile(bam,"rb")
    outfile = pysam.AlignmentFile(out_bam,"wb",template=samfile)
    for read in samfile:
        print(read)


def main():
    parser = argparse.ArgumentParser(description="Filters reads for spades")
    parser.add_argument("-o","--out-prefix", help="Output prefix",default="test/")
    parser.add_argument(dest="bam", help="Input bam")

    args = parser.parse_args()
    filter_bam_and_run_spades(args.bam, args.out_prefix)
    filter_pairs_one_mapping_to_plasmid(args.bam, args.out_prefix)

if __name__=="__main__":
    main()
