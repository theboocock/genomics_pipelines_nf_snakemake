suppressMessages(library(ampliconParser))
suppressMessages(library(stringr))
args = commandArgs(trailingOnly=T)
print(args)
plasmid_in_location = args[1]
otable = args[2]
print(plasmid_in_location)
out_prefix= args[3]
otable = suppressMessages(readRDS(otable))
uguides=suppressMessages(as.character(attr(otable, 'uguides')))
#expected unique repair templates 
utemps=as.character(attr(otable, 'utemps'))
#expected structural sequence 
srt=as.character(attr(otable,'struct'))


in_f = read.delim(plasmid_in_location,sep="\t",header=T)


print(in_f)
print("HERE2")
maxDist.gRNA=2 
#lv edit distance for repair template
maxDist.repTemp=5 
#maximum distance for structural region 
maxDist.struct=5 
#total number of threads
nthreads=8

#ampliconParser::seqtrie_match(
if(nrow(in_f) >0 ){
    grna.matches = ampliconParser::seqtrie_match(in_f$guide[1],uguides,maxDist=maxDist.gRNA,lv.nthreads=8,name="gRNA.") 
    repTemp.matches= ampliconParser::seqtrie_match(in_f$repair[1],utemps,maxDist=maxDist.repTemp,name="repTemp.") 
    structReg.matches = ampliconParser::seqtrie_match(in_f$struct[1],srt,maxDist=maxDist.struct,name="struct.") 
    print(grna.matches)
    print(repTemp.matches)
    print(structReg.matches)
}
expected_combos = attr(otable,"ecombos.str")
observed_combo= paste0(grna.matches$gRNA.expInd,":",repTemp.matches$repTemp.expInd)


matches = sum(observed_combo %in% expected_combos)
matches_idx = which(expected_combos %in% observed_combo)
print(matches_idx)
if(matches == 1){
    # Matches an expected combo 
    idx_chrom = otable$guideIndex[matches_idx]
    chrom = str_split(idx_chrom, ":")[[1]][1]
    out_df = data.frame(chrom=chrom,observed_combos = observed_combo)
}else{
    out_df = data.frame(observed_combos = "no_match")
}
out_file_test = paste(out_prefix,"combo_match.txt",sep="/")
print(out_file_test)
write.table(out_df, file=out_file_test, col.names=F,row.names=F,quote=F)


