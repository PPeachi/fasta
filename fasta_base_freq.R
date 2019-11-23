#source("http://bioconductor.org/biocLite.R")
#biocLite("Biostrings") #安装Biostrings
#library(Biostrings) ;
#s1 = readDNAStringSet("flu_seq.fas")
#readBStringSet(filepath, format="fasta", nrec=-1L, skip=0L, seek.first.rec=FALSE, use.names=TRUE)
#install.packages("seqinr")
#library(seqinr); 
#fastafile1 <- read.fasta(file = "flu_seq.fas")
#library(ape);
#read.dna();
#read.FASTA

###1
library(seqinr)
fastafile1 <- read.fasta(file = "flu_seq.fas")
fastafile2 <- read.fasta(file = "flu_seq_v2.fas")
fastafile1
fastafile2
id1<-attr(fastafile1,"name")
id2<-attr(fastafile2,"name")
length(id1)
length(id2)
id1[1]
unlist(strsplit(id1[1],split = "[|]"))
#x<-"a.b.c";y<-gsub("[.]","-",x)
x<-c()
y<-c()
for (i in 1:length(id1)){
  x[i]<-unlist(strsplit(id1[i],split = "[|]"))[1]
  y[i]<-unlist(strsplit(id2[i],split = "[|]"))[1]
}
fastafile1[[1]]
length(fastafile1[[1]])
seq1<-c()
seq2<-c()
for (i in 1:length(fastafile1)){
  seq1[i]<-paste(unlist(fastafile1[[i]]),collapse = "")
  seq2[i]<-paste(unlist(fastafile2[[i]]),collapse = "")
}
df1<-data.frame(id=x,seq=seq1)
df2<-data.frame(id=y,seq=seq2)
df1$id<-as.character(df1$id)
df2$id<-as.character(df2$id)
df1$seq<-as.character(df1$seq)
df2$seq<-as.character(df2$seq)
for (i in 1:11){
  count_A<-0
  count_T<-0
  count_G<-0
  count_C<-0
  s<-unlist(strsplit(df1$seq[i],""))
  for (j in 1:length(s)){
    if(s[j]=='a'){count_A=count_A+1}
    if(s[j]=='t'){count_T=count_T+1}
    if(s[j]=='g'){count_G=count_G+1}
    if(s[j]=='c'){count_C=count_C+1}
  }
  m<-length(s)
  df1$a[i]<-count_A/m
  df1$t[i]<-count_T/m
  df1$g[i]<-count_G/m
  df1$c[i]<-count_C/m
}
for (i in 1:11){
  count_A<-0
  count_T<-0
  count_G<-0
  count_C<-0
  s<-unlist(strsplit(df2$seq[i],""))
  for (j in 1:length(s)){
    if(s[j]=='a'){count_A=count_A+1}
    if(s[j]=='t'){count_T=count_T+1}
    if(s[j]=='g'){count_G=count_G+1}
    if(s[j]=='c'){count_C=count_C+1}
  }
  m<-length(s)
  df2$a[i]<-count_A/m
  df2$t[i]<-count_T/m
  df2$g[i]<-count_G/m
  df2$c[i]<-count_C/m
}
tb1<-df1
tb2<-df2
tb1<-tb1[,-2]
tb2<-tb2[,-2]
rownames(tb1)<-tb1[,1]
rownames(tb2)<-tb2[,1]
tb1<-tb1[,-1]
tb2<-tb2[,-1]
pheatmap::pheatmap(tb1)
pheatmap::pheatmap(tb2)

###2
library(seqinr) 
fastafile <- read.fasta(file = "flu_seq.fas")
base_freq<-function(fasta){
  id<-attr(fasta,"name")
  ids<-c()
  for (i in 1:length(id)){
    ids[i]<-unlist(strsplit(id[i],split = "[|]"))[1]
  }
  seqs<-c()
  for (i in 1:length(fasta)){
    seqs[i]<-paste(unlist(fasta[[i]]),collapse = "")
  }
  df<-data.frame(id=ids,seq=seqs)
  df$id<-as.character(df$id)
  df$seq<-as.character(df$seq)
  for (i in 1:length(df$id)){
    count_A<-0
    count_T<-0
    count_G<-0
    count_C<-0
    s<-unlist(strsplit(df$seq[i],""))
    for (j in 1:length(s)){
      if(s[j]=='a'){count_A=count_A+1}
      if(s[j]=='t'){count_T=count_T+1}
      if(s[j]=='g'){count_G=count_G+1}
      if(s[j]=='c'){count_C=count_C+1}
    }
    m<-length(s)
    df$a[i]<-count_A/m
    df$t[i]<-count_T/m
    df$g[i]<-count_G/m
    df$c[i]<-count_C/m
  }
  rownames(df)<-df[,1]
  tb<-df[,-(1:2)]
  return (tb)
}
res<-base_freq(fastafile)
res
pheatmap::pheatmap(res)

###function1
library(seqinr)
base_freq <- function(fasta){
  #fastafile <- read.fasta(fasta)
  id<-attr(fasta,"name")
  ids<-c()
    for (i in 1:length(id)){
      ids[i]<-unlist(strsplit(id[i],split = "[|]"))[1]
    }
  seqs<-c()
  len<-c()
  for (i in 1:length(fasta)){
    seqs[i]<-paste(unlist(fasta[[i]]),collapse = "")
    len[i]<-nchar(seqs[i])
  }
  df<-data.frame(id=ids,seq=seqs,length=len)
  df$id<-as.character(df$id)
  df$seq<-as.character(df$seq)
  df$length<-as.numeric(df$length)
  for (i in 1:length(df$id)){
    count_A<-0
    count_T<-0
    count_G<-0
    count_C<-0
    s<-unlist(strsplit(df$seq[i],""))
    for (j in 1:length(s)){
      if(s[j]=='a'){count_A=count_A+1}
      if(s[j]=='t'){count_T=count_T+1}
      if(s[j]=='g'){count_G=count_G+1}
      if(s[j]=='c'){count_C=count_C+1}
    }
    m<-length(s)
    df$a[i]<-count_A/m
    df$t[i]<-count_T/m
    df$g[i]<-count_G/m
    df$c[i]<-count_C/m
  }
  rownames(df)<-df[,1]
  tb<-df[,-(1:2)]
  res<-df[,-(1:3)]
  writeLines(paste0("#Base composition:"))
  print(tb)
  return (res)
}
read_fasta<-function(file){
  fastafile <- read.fasta(file)
  seq_num <- length(fastafile)
  cat(seq_num,"sequences in total","\n","\n")
  id<-attr(fastafile,"name")
  cat("#Labels:","\n")
  writeLines(id)
  cat("\n")
  (res<-base_freq(fastafile))
  writeLines(paste0("\n","#id:"))
  writeLines(paste(1:nrow(res),rownames(res),sep = "."))
  return (res)
}
f<-read_fasta("flu_seq.fas")
pheatmap::pheatmap(f)



###function2
line_1<-readLines("flu_seq.fas")
line_2<-readLines("flu_seq_v2.fas")
i_1<-grep(">",line_1)
i_2<-grep(">",line_2)
id_1<-sub(">","",line_1[i_1])
id_2<-sub(">","",line_2[i_2])
start_1<-(i_1+1)
start_2<-(i_2+1)
end_1<-c(i_1[-1]-1,length(line_1))
end_2<-c(i_2[-1]-1,length(line_2))
seq_1<-sapply(seq_along(start_1), function(i) paste0(line_1[start_1[i]:end_1[i]],collapse = ""))
seq_2<-sapply(seq_along(start_2), function(i) paste0(line_2[start_2[i]:end_2[i]],collapse = ""))
df1<-data.frame(id=id_1,seq=seq_1)
df2<-data.frame(id=id_2,seq=seq_2)
cat(line_1,sep="\n")
cat(line_2,sep="\n")



###final: function3
read_fas<-function(file){
  line<-readLines(file)
  i<-grep(">",line)
  id<-sub(">","",line[i])
  start<-(i+1)
  end<-c(i[-1]-1,length(line))
  seq<-sapply(seq_along(start), function(i) paste0(line[start[i]:end[i]],collapse = ""))
  df<-data.frame(id=id,seq=seq)
  df$id<-as.character(df$id)
  df$seq<-as.character(df$seq)
  return (df)
}
x<-read_fas("flu_seq.fas")
base_freq<-function(fasta){
  seq_num <- length(fasta$id)
  cat(seq_num,"sequences in total","\n","\n")
  cat("#Labels:","\n")
  id<-fasta$id
  writeLines(id)
  ids<-c()
  len<-c()
  for (i in 1:seq_num){
    ids[i]<-unlist(strsplit(id[i],split = "[|]"))[1]
    len[i]<-nchar(fasta$seq[i])
  }
  df<-data.frame(id=ids,length=len)
  for (i in 1:seq_num){
    count_A<-0
    count_T<-0
    count_G<-0
    count_C<-0
    s<-toupper(unlist(strsplit(fasta$seq[i],"")))
    for (j in 1:length(s)){
      if(s[j]=='A'){count_A=count_A+1}
      if(s[j]=='T'){count_T=count_T+1}
      if(s[j]=='G'){count_G=count_G+1}
      if(s[j]=='C'){count_C=count_C+1}
    }
    m<-length(s)
    df$A[i]<-count_A/m
    df$T[i]<-count_T/m
    df$G[i]<-count_G/m
    df$C[i]<-count_C/m
  }
  rownames(df)<-df[,1]
  tb<-df[,-1]
  res<-df[,-(1:2)]
  writeLines(paste0("\n","#Base composition:"))
  print(tb)
  return (res)
}
y<-base_freq(x)
y
pheatmap::pheatmap(y)



#if we need to merge 2 or more sequences
x1<-read_fas("AB115403.fas")
x2<-read_fas("AJ534526.fas")
x3<-read_fas("AJ534527.fas")
total<-rbind(x1,x2,x3)
y<-base_freq(total)
y
pheatmap::pheatmap(y)



###attention
res<-treeio::read.fasta("flu_seq.fas")
ape::base.freq(res)
write.fasta <- function(x, file) {
  ape::write.FASTA(x, file)
}