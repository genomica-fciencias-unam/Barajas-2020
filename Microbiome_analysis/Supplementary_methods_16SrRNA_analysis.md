
#Bioinformatic supplementary methods. Barajas et al., 2020.
#16S rRNA amplicon sequence processing and analysis.

#Sequence preprocessing.
Raw sequence quality analysis, assembly, filtering, and indexing of sequences.


```python
#Raw sequence quality was inspected using FastxToolkit.
fastx_quality_stats -i muestra.foo.fastq.gz -o muestra.foo.fastq.qstat
fastq_quality_boxplot_graph.sh -i muestra.foo.qstat -t muestra -o muestra.foo.qstat.png
```


```python
#In order to avoid the low quality in the 3' end  of the forward and reverse sequencing reads, sequences were trimmed to a 250 bp length
#Sequence trimming was performed with FastxToolkit to a length of 250bp.
#Paired end assembly was performed with Pandaseq with the following script.
#Assembly parameters are set as follows: 95% threshold of correct base assignation, over a 250bp minimal length and 470bp maximal length, and a minimum overlap of 15bp. 

#!/bin/bash

SEQS=/home/ldalcaraz/run_ago_15
SALIDAS=/home/ldalcaraz/salidas
BIN=/home/ldalcaraz/bin
BIN2=/home/ldalcaraz/bin
COUNT=0

for FAA in `cat lista`   #lista is a a text file with a list of raw reads files.
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

#Trimming Forward reads
echo "zcat $SEQS/$FAA"_R1.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr
_R1.fastq" &" >>$*.$COUNT.scr
#Trimming Reverse reads
echo "zcat $SEQS/$FAA"_R2.fastq.gz" | $BIN2/fastx_trimmer -l 250 > $SEQS/$FAA"tr
_R2.fastq"" >>$*.$COUNT.scr
#Assembly of Forward and Reverse reads.
echo "$BIN/pandaseq -B -f $SEQS/$FAA"tr_R1.fastq" -r $SEQS/$FAA"tr_R2.fastq" -t
0.95 -l 250 -L 470 -o 15 -w $SALIDAS/$FAA"_4nov15.fasta" -G $SALIDAS/$FAA"-4nov
15.log.bz2"" >>$*.$COUNT.scr    #Assembly of Forward and Reverse reads. 

perl remove_small.pl 250 assembled_seqs.fasta >assembled_seqs.f250.fasta   #Remove sequences shorter than 250bp.
```


```python
#Assembled sequences of each sample were indexed with the following script: header.fasta.numbers.pl

#!/usr/bin/perl
# este script cambia cada uno de los identificadores de un archivo fasta a >prefix_# donde prefix es el
# primer argumento del script y # es una numeración de cada una de las secuencias que componen el archivo el segundo argumento es el nombre del archivo fasta a renombrar.
# este script es útil si se quiere usar QIIME en archivos fasta desmultiplexados o individuales. 
# uso: ./header.fasta.numbers.pl prefix nombre_del_archivo.fasta
# Luis David Alcaraz 2013-04-11

my $prefix = $ARGV[0]; chomp $prefix; 
my $f =  1;

my $fasta_file = $ARGV [1]; chomp $fasta_file;

my $fh;
open($fh, $fasta_file) or die "can't open $fasta_file: $!\n";
open(OUT, ">$fasta_file.numbered.fas") || die "can't open $fasta_file.numbered.fas\n";

my %sequence_data;
while (read_fasta_sequence($fh, \%sequence_data)) {
   print OUT ">$sequence_data{header}\n$sequence_data{seq}\n";
}

close $fh; 
close OUT;

sub read_fasta_sequence {
   my ($fh, $seq_info) = @_;

   $seq_info->{seq} = undef; # clear out previous sequence

   # put the header into place
   $seq_info->{header} = $seq_info->{next_header} if $seq_info->{next_header};

   my $file_not_empty = 0; 
   while (<$fh>) {
      $file_not_empty = 1;
      next if /^\s*$/;  # skip blank lines
      chomp;    

      if (/^>/) { # fasta header line
         my $h = $_;    
         $h =~ s/>/$prefix\_$f\ /; 
    $f++; 
         if ($seq_info->{header}) {
            $seq_info->{next_header} = $h;
            return $seq_info;
   
         }              
         else { # first time through only
            $seq_info->{header} = $h;
         }              
      }         
      else {    
         s/\s+//;  # remove any white space
         s/\n\n/\n/; 
         $seq_info->{seq} .= $_;
      }         
   }    

   if ($file_not_empty) {
      return $seq_info;
   }    
   else {
      # clean everything up
      $seq_info->{header} = $seq_info->{seq} = $seq_info->{next_header} = undef;

      return;   
   }    
}
```


```python
#Indexing 16S rRNA gene sequences in every sample with the latter script.
./header.fasta.numbers.pl SAMPLE_NAME Sample.fasta
```

#16S rRNA Data analysis. Taxonomy assignments and OTU tables


```python
#Indexed sequences from every sample were concatenated into a single file and clustered into OTUs at 97% sequence identity with cd-hit-est.
cat sample_1.numbered.fasta sample_2.numbered.fasta sample_n.fasta >>soil_tom_loc.fasta
cd-hit-est -i soil_tom_loc.fasta -o soil_tom_loc.fasta.cl -c 0.97 -M 4000

#Cluster table was edited in order to get a table in which the first column is the cluster number and the second is the representative sequence identifier.
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_\n"; exit}' soil_tom_loc.fasta.cl.clstr >soil_tom_loc.otu
sed -i '1d' soil_tom_loc.otu

#Parsing the otu table with the indexed fasta file containing all 16S rRNA gene sequences from all sampels to obtain the representative sequences of the clusters using QIIME scripts.
pick_rep_set.py -i soil_tom_loc.otu -f soil_tom_loc.fasta -o rep_set_stl.fasta
```


```python
#Taxonomy assignments of the representative sequences were performed with QIIME against the Greengenes database.
parallel_assign_taxonomy_blast.py -i rep_set_stl.fasta -o taxonomy -r /qiime/gg_otus-13_8-release/rep_set/97_otus.fasta -t /qiime/gg_otus-13_8-release/taxonomy/97_otu_taxonomy.txt

#Taxonomy table edition. Substitute spaces with underscores.
sed -i 's/ /_/g;s/;_/; /g;s/__/_/g;s/\]//g;s/\[//g'  rep_set_stl_tax_assignments.txt
```


```python
#A OTU table was created with the edited cluster table and taxonomy table using QIIME scripts.
make_otu_table.py soil_tom_loc.otu -t rep_set_stl_tax_assignments.txt -o soil_tom_loc.biom

#Removing singletons from the OTU table in BIOM format using the QIIME script.
filter_otus_from_otu_table.py -i soil_tom_loc.biom -o soil_tom_loc_wos.biom -n 2
    
#OTU table was converted from BIOM format to tabular format.
biom convert --to-tsv -i soil_tom_loc_wos.biom -o soil_tom_loc_wos.txt --table-type "Taxon table" --header-key=taxonomy

#Tabular OTU table was filtered to discard OTUs assigned as mitochondria, chloroplast, and sequences without hits in the greengenes database.
grep -iv "no_blast_hit" soil_tom_loc_wos.txt | grep -iv "chloroplast" | grep -iv "mitochondria" >>soil_tom_loc_wos.clean.txt

#Creating the taxonomy table and the OTU count table from the filtered table in last step.
#Taxonomy table.
perl -pe 's/\; /\;/g' soil_tom_loc_wos.clean.txt | awk '{print $1,"\t",$NF}' | perl -pe 's/\;/\t/g' >soil_tom_loc_wos.clean.tax
#OTU count table.
rev soil_tom_loc_wos.clean.txt | cut -f2- | rev >soil_tom_loc_wos.clean.otu.txt

#Cleaning up headers in the OTU counts table.
grep -v "biom file" soil_tom_loc_wos.clean.otu.txt | perl -pe 's/\#OTU\ ID/OTU_ID/g' >soil_tom_loc_wos.clean.otu1.txt

```

#Microbiome profiling and diversity analyses.
Loading datasets in R for taxonomic diversity and abundance analysis with phyloseq


```python
#Taxonomic profiling of soil, tomatoes, and ruderal plants microbiomes through QIIME scripts.


#Create an abundance table at the phylum level with QIIME script using the singleton and contaminants depleted file. The resulting phylum level table was edited in libreOffice Calc before loading data into R kernel. Proteobacteria abundances were divided into Class level.
biom convert -i soil_tom_loc_wos.clean.txt -o soil_tom_loc_wos.clean.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
summarize_taxa.py -i soil_tom_loc_wos.clean.biom -o barplot

#Load required libraries in R.
library(ggplot2)
library(RColorBrewer)

#Samplewise relative abundance barplots
#Create a data frame with a phylum/class abundance table
phylum3<-read.table("phylum_barplot_proteoclass.txt", sep="\t", header=T, row.names=1)
abund_table3<-t(phylum3)

#Ordering the abundances in decreasing order.
x3 <- abund_table3[,order(colSums(abund_table3), decreasing=T)]
for (i in 1:dim(x3)[2]){
tmp3 <- data.frame(row.names=NULL,Sample=rownames(x3),Taxa=rep(colnames(x3)[i],dim(x3)[1]),Value=x3[,i])
if(i==1){df3<-tmp3} else {df3<-rbind(df3,tmp3)}
}
#Set the order of appereance of samples in X and plot a stacked barchart of Phyla and Class (for Proteobacteria) .
df3$orden <- factor(df3$Sample, levels = c("ZAC1ECT", "ZAC1RFT", "DGO1RZ", "NAY3RFT", "NAY3ECT", "NAY2ECT", "ZAC1SF", "ZAC1US", "ZAC1S0", "NAY2RZ", "NAY2EC", "NAY3RZ", "NAY3US", "NAY3SF", "NAY3S0", "NAY2US", "NAY2SF", "NAY2S0", "NAY2RFT", "JAL2US", "JAL2SF", "JAL2RFT", "JAL5US", "JAL4US", "JAL4SF", "JAL4RFT", "JAL1US", "JAL1SF", "JAL1RFT", "DGO1SF", "DGO1US", "SIN2SF", "SIN2S0", "SIN2US", "GTO2US", "GTO2SF", "GTO2S0", "GTO1US", "GTO1SF", "GTO1S0", "JAL5S0", "JAL5SF", "AGS1US", "SLP1SF", "AGS1S0", "GTO3S0", "GTO3US", "GTO3SF", "SIN1S0", "SIN1SF", "SIN1US", "JAL2S0", "DGO1S0", "JAL3US", "JAL3SF", "JAL3S0", "JAL4S0", "JAL1S0", "SLP1S0", "GTO1RZ", "GTO1EC", "SIN1RZ", "GTO3RZ", "JAL4RZ", "JAL1RZ", "JAL3EC", "JAL3RZ", "JAL5RZ", "JAL5EC", "GTO3EC", "AGS1EC", "ZAC1RZ", "JAL1EC", "SLP1RZ", "SLP1EC", "SLP1US", "SIN2EC", "SIN1EC", "SIN2RZ", "JAL2RZ", "JAL2EC", "GTO2RZ", "GTO2EC", "JAL3ECT", "JAL1ECT", "DGO1ECT", "JAL5ECT", "JAL5RFT", "GTO2ECT", "AGS1ECT", "GTO3RFT", "GTO3ECT", "GTO1ECT", "GTO1RFT", "JAL2ECT", "JAL4ECT", "SIN2ECT", "SIN2RFT", "AGS1RFT", "GTO2RFT", "SLP1ECT", "SLP1RFT", "DGO1RFT", "JAL3RFT", "SIN1ECT", "SIN1RFT"))
p.todos5 <- ggplot (df3, aes(x=df3$orden, Value, fill=Taxa)) + geom_bar(stat="identity")
p.todos5 + scale_fill_manual(values=rev(colorRampPalette(brewer.pal(16, "Paired"))(16))) + theme_bw() +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

#Grouped relative abundance barplots
#Create a data frame with a phylum abundance table
phylum<-read.table("phylum_barplot_groups.txt", sep="\t", header=T, row.names=1)
abund_table<-t(phylum)

#Ordering the abundances in decreasing order.
x <- abund_table[,order(colSums(abund_table), decreasing=T)]
for (i in 1:dim(x)[2]){
tmp <- data.frame(row.names=NULL,Sample=rownames(x),Taxa=rep(colnames(x)[i],dim(x)[1]),Value=x[,i])
if(i==1){df<-tmp} else {df<-rbind(df,tmp)}
}
#Set the order of appereance of samples in X and plot a stacked barchart of Phyla and Class (for Proteobacteria) .
df$orden <- factor(df$Sample, levels = c("SI", "ZAC1RFT", "DGO1RZ", "NAY3RFT", "NAY3ECT", "NAY2ECT", "ZAC1SF", "ZAC1US", "ZAC1S0", "NAY2RZ", "NAY2EC", "NAY3RZ", "NAY3US", "NAY3SF", "NAY3S0", "NAY2US", "NAY2SF", "NAY2S0", "NAY2RFT", "JAL2US", "JAL2SF", "JAL2RFT", "JAL5US", "JAL4US", "JAL4SF", "JAL4RFT", "JAL1US", "JAL1SF", "JAL1RFT", "DGO1SF", "DGO1US", "SIN2SF", "SIN2S0", "SIN2US", "GTO2US", "GTO2SF", "GTO2S0", "GTO1US", "GTO1SF", "GTO1S0", "JAL5S0", "JAL5SF", "AGS1US", "SLP1SF", "AGS1S0", "GTO3S0", "GTO3US", "GTO3SF", "SIN1S0", "SIN1SF", "SIN1US", "JAL2S0", "DGO1S0", "JAL3US", "JAL3SF", "JAL3S0", "JAL4S0", "JAL1S0", "SLP1S0", "GTO1RZ", "GTO1EC", "SIN1RZ", "GTO3RZ", "JAL4RZ", "JAL1RZ", "JAL3EC", "JAL3RZ", "JAL5RZ", "JAL5EC", "GTO3EC", "AGS1EC", "ZAC1RZ", "JAL1EC", "SLP1RZ", "SLP1EC", "SLP1US", "SIN2EC", "SIN1EC", "SIN2RZ", "JAL2RZ", "JAL2EC", "GTO2RZ", "GTO2EC", "JAL3ECT", "JAL1ECT", "DGO1ECT", "JAL5ECT", "JAL5RFT", "GTO2ECT", "AGS1ECT", "GTO3RFT", "GTO3ECT", "GTO1ECT", "GTO1RFT", "JAL2ECT", "JAL4ECT", "SIN2ECT", "SIN2RFT", "AGS1RFT", "GTO2RFT", "SLP1ECT", "SLP1RFT", "DGO1RFT", "JAL3RFT", "SIN1ECT", "SIN1RFT"))
p.todos <- ggplot (df, aes(x=df$orden, Value, fill=Taxa)) + geom_bar(stat="identity")
p.todos + scale_fill_manual(values=rev(colorRampPalette(brewer.pal(16, "Paired"))(16))) + theme_bw() +theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))


```

Phylogenetic tree of OTUs


```python
#Calculation of the 16S rRNA phylogenetic tree using OTU reference sequences. 
#Sequence identifiers of the singletons filtered dataset were obtained and sequences were extracted from the multifasta file containing all samples sequences.

#Extract sequence identifiers from OTU tax table.
cut -f1 soil_tom_loc_wos.clean.tax >ids_para_unifrac.txt

#Extract OTU representative sequences from representative sequences fasta file.
perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1wif @ARGV' ids_para_unifrac.txt rep_set_stl.fasta >>seqs_unifrac.fasta

#Alignment of sequences was performed with SSU-align. First a script was created with ssu-prep to run the job in parallel on an HPC. The resulting script was run
ssu-prep seqs_unifrac.fasta align 100 /home/hugo/download_programs/ssu-align/ssu-align-0.1/tutorial/janelia-cluster-presuf.txt
#Convert the Stockholm format into an aligned FASTA format.
ssu-mask --stk2afa align.bacteria.stk /home/hugo/soil_tom_loc/alignment/align/align.bacteria.afa
#Multiple alignment file clean up. Susbtituting dots for dashes.
sed -i 's/\./\-/g' align.bacteria.afa
#Calculation of the phylogenetic tree using Fasttree
fasttree -nt -fastest align.bacteria.afa > stl_wos_otus.bacteria.tre

```

Load OTU count table, tax table, metadata and phylogenetic tree for downsteam analysis.


```python
#R Libraries required for the analysis. Phyloseq, ggplot2, limma, ape, vegan, colorRamps

#Load required libraries.
source("https://bioconductor.org/biocLite.R")
library(phyloseq)
library(ggplot2)
library(ape)
biocLite("limma")

#Load OTU count table and create phyloseq otu object.
otu<-as.matrix(read.table("soil_tom_loc_wos.clean.otu1.txt", header=T, row.names=1))
OTU = otu_table(otu, taxa_are_rows=T)

#Load taxonomy table and create phyloseq taxonomy object.
taximat = as.matrix(read.table("soil_tom_loc_wos.clean.tax", header=T, row.names=1))
taxi=tax_table(taximat)

#Create phyloseq object from OTU table and taxonomy table.
suelo.phy<-phyloseq(OTU, taxi)

#Load metadata table and create sample data object.
suelos_data=read.table("metadatos_tomates_ordenado.txt", header=TRUE, row.names=1, sep="\t")
sampledata=sample_data(data.frame(id=suelos_data$id.1, sitio=suelos_data$site, estado=suelos_data$estado, nombre=suelos_data$name, muestra=suelos_data$sample_type, sc=suelos_data$soil_class, tiempo=suelos_data$tiempo, ph=suelos_data$ph, nt=suelos_data$Nt, pt=suelos_data$Pt, ct=suelos_data$Ct, c_n=suelos_data$C_N, c_p=suelos_data$C_P, n_p=suelos_data$N_P, bm=suelos_data$biomass, alt=suelos_data$altitude, ia=suelos_data$IA, latitud=suelos_data$latitude, longitud=suelos_data$longitude, pol=suelos_data$polimerase, conc_amp=suelos_data$conc_amp, abs_amp=suelos_data$abs_amp, conc_dna=suelos_data$conc_dna, abs_dna=suelos_data$abs_dna, row.names=sample_names(suelo.phy)))

#Loading the tree in R kernel and calculating UniFrac distances 
read.tree("stl_wos_otus.bacteria.tre")->stl_wos_tree

#Add metadata to the phyloseq object.
stl.phy<-phyloseq(OTU, taxi, sampledata, stl_wos_tree)

#The final soil (SF) sampĺe of site AGS1 was discarded due to sequencing errors. 
stl.ss.phy<-subset_samples(suelo.phy, id!="AGS1SF")
```

#Alpha diversity 


```python
#Diversity metrics were calculated with Phyloseq R package.

#Obtain alpha diversity metrics.
rich_all=estimate_richness(stl.ss.phy, measures=c("Observed", "Chao1", "Simpson", "Shannon"))

#Shannon diversity was plotted as a heatmap using a table in which samlples are ordered according to the beta-diversity dendrogram (See beta-viversity calculations).
shan <- read.table("shannon_ordered.txt", header=TRUE, row.names=1)
x <- shan$SHAN
shan$id.1<-rownames(shan)
x <- "Sample"
y <- shan$id.1
y$name <- factor(shan$id.1, levels= shan$id.1)
ggplot(shan, aes(x, y = y$name, fill=SHAN)) + geom_tile() + scale_y_discrete(limits=rev(levels(y))) + scale_fill_gradient(low="white", high="red")

#Plot diversity boxplots for each sample type groups using a custom table with groups metadata.
read.table("div_amplicones_bplot.tsv", header=TRUE, sep="\t", row.names=1)->div.amplicones
ggplot(data=div.amplicones, aes(x=Grupo, y=Value)) + geom_boxplot(aes(fill=Variable)) +geom_point(aes(y=Value, group=Label), position=position_dodge(width=0.75)) + theme_bw()

#Alpha diversity comparison between soil, tomatoes and ruderal plants. The ggsignif package was used to estimate pairwise significance differences between grous using a T-test. 
plot_richness(stl.ss.phy, measures=c("Observed", "Shannon", "Chao1"), x="samp" ) + geom_boxplot() + geom_point(size=3) + geom_signif(comparisons=list(c("EC", "ECT"), c("RFT", "RZ"), c("RFT", "S0"), c("RZ","S0"), c("SF","S0")), map_signif_level = TRUE) + theme_light() + theme(axis.text.x=element_text(angle=90, hjust=1)) -> p.samp.rich.bplot.sign.stl.ss.phy

```


```python
#Comparison of Phylum level abundance between soils, ruderal plants and tomatoes

#Agglomeration of OTUs based on phylum and relative abundance calculation of counts.
Phylum_stl <- tax_glom(stl.phy, taxrank="Phylum")
Phylum_stl_rel <- transform_sample_counts(Phylum_stl, function(x) x / sum(x))

#Join tax table, OTU table, and metadata table
cphylum_stl_rel <- data.frame(cbind(tax_table(Phylum_stl_rel), otu_table(Phylum_stl_rel)), row.names=NULL)
cphylum_stl_rel <- cphylum_stl_rel[, colSums(is.na(cphylum_stl_rel)) != nrow(cphylum_stl_rel)]
tcphylum_stl_rel <- data.frame(t(cphylum_stl_rel))
tcphylum_stl_rel <- tcphylum_stl_rel[-1,]
colnames(tcphylum_stl_rel) <- NULL

#Export and load table
write.table(tcphylum_stl_rel, "tcphylum_stl_rel.tsv", quote = TRUE, sep = "\t")
tcphylum_stl_rel <- read.table("tcphylum_stl_rel.tsv", header = TRUE, sep = "\t")

#Add metadata.
microbiome <- data.frame(suelos_data2$microbiome)
mtcphylum_stl_rel <- cbind(microbiome, tcphylum_stl_rel, row.names = NULL)
mtcphylum_stl_rel


#Plot relative abundance boxplots
melt_phylum_stl_rel <- melt(mtcphylum_stl_rel)
colnames(melt_phylum_stl_rel) <- c("Microbiome_type", "Sample", "Phylum", "Relative_abundance")

#Change variable names.
melt_phylum_stl_rel$Microbiome_type <- sub("Plant_local", "Ruderal plants", melt_phylum_stl_rel$Microbiome_type)
melt_phylum_stl_rel$Microbiome_type <- sub("Plant", "S. lycopersicum", melt_phylum_stl_rel$Microbiome_type)
melt_phylum_stl_rel$Microbiome_type <- sub("Soil", "Soils", melt_phylum_stl_rel$Microbiome_type)

gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

plot_phylum_stl_rel <- ggplot(melt_phylum_stl_rel, aes(Phylum, Relative_abundance, fill=Microbiome_type)) + 
                            geom_boxplot() + theme_light(base_size = 22) +
                            theme(axis.text.y = element_text(size = 20), 
                                  axis.text.x = element_text(angle = 90, size = 20)) +

ggsave("plot_phylum_stl_rel.pdf", width=80, height=40, units="cm")

#ANOVA for each group
apr <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Proteobacteria"))
ave <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Verrucomicrobia"))
afi <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Firmicutes"))
aac <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Actinobacteria"))
aar <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Armatimonadetes"))
ach <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Chloroflexi"))
aci <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Acidobacteria"))
apl <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Planctomycetes"))
aba <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Bacteroidetes"))
age <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Gemmatimonadetes"))
atm <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_TM7"))
acy <- aov(Relative_abundance~Microbiome_type, filter(fmelt_phylum_stl_rel, Phylum == "p_Cyanobacteria"))

#Tukey test for each comparison.
options(digits=22)
TukeyHSD(apr)
TukeyHSD(ave)
TukeyHSD(afi)
TukeyHSD(aac)
TukeyHSD(aar)
TukeyHSD(ach)
TukeyHSD(aci)
TukeyHSD(apl)
TukeyHSD(aba)
TukeyHSD(age)
TukeyHSD(atm)
TukeyHSD(acy)
options(digits=5)


#Comparison of Class level abundance for Proteobacteria between soils, ruderal plants and tomatoes
Class_stl <- tax_glom(stl.phy, taxrank="Class")
Class_stl_rel <- transform_sample_counts(Class_stl, function(x) x / sum(x))

#Join tax table, OTU table, and metadata table for Proteobacteria classes.
cclass_stl_rel <- data.frame(cbind(tax_table(Class_stl_rel), otu_table(Class_stl_rel)), row.names=NULL)
cclass_stl_rel <- cclass_stl_rel[, colSums(is.na(cclass_stl_rel)) != nrow(cclass_stl_rel)]

tcclass_stl_rel <- data.frame(t(cclass_stl_rel))
tcclass_stl_rel <- tcclass_stl_rel[-1,]
tcclass_stl_rel <- tcclass_stl_rel[-1,]
colnames(tcclass_stl_rel) <- NULL

#Export and load table 
write.table(tcclass_stl_rel, "tcclass_stl_rel.tsv", quote = TRUE, sep = "\t")
tcclass_stl_rel <- read.table("tcclass_stl_rel.tsv", header = TRUE, sep = "\t")
#Add metadata
microbiome <- data.frame(suelos_data2$microbiome)
mtcclass_stl_rel <- cbind(microbiome, tcclass_stl_rel, row.names = NULL)
mtcclass_stl_rel


#Plot relative abundance of Alphaproteobacteria.
melt_class_stl_rel <- melt(mtcclass_stl_rel)
colnames(melt_class_stl_rel) <- c("Microbiome_type", "Sample", "Class", "Relative_abundance")

class_alpha <- filter(melt_class_stl_rel, Class == "c_Alphaproteobacteria" | 
                               Class == "c_Betaproteobacteria" | Class == "c_Deltaproteobacteria" |
                               Class == "c_Gammaproteobacteria")
#change variable names
class_alpha$Microbiome_type <- sub("Plant_local", "Ruderal plants", class_alpha$Microbiome_type)
class_alpha$Microbiome_type <- sub("Plant", "S. lycopersicum", class_alpha$Microbiome_type)
class_alpha$Microbiome_type <- sub("Soil", "Soils", class_alpha$Microbiome_type)

#Order appeareance of plots.
class_alpha$Class_o = factor(class_alpha$Class, levels=c('c_Alphaproteobacteria', 
                                        'c_Betaproteobacteria', 'c_Deltaproteobacteria', 
                                        'c_Gammaproteobacteria'))

#Plot graphics
plot_class_proteo <- ggplot(class_alpha, aes(Class, Relative_abundance, fill=Microbiome_type)) + 
                            geom_boxplot() + facet_wrap(~ Class_o, scales="free", nrow = 1) + 
                            scale_fill_manual(values = c("#009E73", "#E69F00", "#CC79A7"))+
                            theme_bw(base_size = 30) + theme(strip.background = element_blank(), 
                                                             strip.text = element_blank()) +
                            theme(axis.text.y = element_text(size = 22), 
                                  axis.text.x = element_text(size = 22)) 

ggsave("plot_class_proteo.pdf", width=80, height=13.3, units="cm")


#ANOVA for each group (Proteobacteria Classes)
aal <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Alphaproteobacteria"))
abe <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Betaproteobacteria"))
ade <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Deltaproteobacteria"))
aga <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Gammaproteobacteria"))

#Tukey test for each ANOVA
options(digits=22)
TukeyHSD(aal)
TukeyHSD(ade)
TukeyHSD(abe)
TukeyHSD(aga)
options(digits=5)
```

#Soil, tomato and ruderal plants Core microbiome analysis. 


```python
#Venn diagram comparison using UpSetR package. Venn diagrams were used to describe the core microbiome of ruderal plants and tomatoes at the genus level.

#Load library
library(UpSetR)

#Agglomerate OTUs based on genus
stl.ss.phy.genus<-tax_glom(stl.ss.phy, taxrank="Genus")

#Make a subset of each sample type taking into account the complete root system(tomato, soils, or ruderal plants). Filter taxa present in all samples
tomate.phy<-filter_taxa(subset_samples(stl.ss.phy.genus, mb=="Plant"), function (x) {sum(x > 0) > 0}, prune=TRUE)
locales.phy<-filter_taxa(subset_samples(stl.ss.phy.genus, mb=="Plant_local"), function (x) {sum(x > 0) > 0}, prune=TRUE)
suelos.phy<-filter_taxa(subset_samples(stl.ss.phy.genus, mb=="Soil"), function (x) {sum(x > 0) > 0}, prune=TRUE)

#Convert OTU abundance table into presence-absence table in each sample type table.
otomatef<-otu_table(tomate.phy)
otomatef[otomatef>0] = 1

olocalesf<-otu_table(locales.phy)
olocalesf[olocalesf>0] = 1

osuelosf<-otu_table(suelos.phy)
osuelosf[olocalesf>0] = 1

#Perform Venn diagram analysis using UpSetR package
otomatef_ups <- upset(as.data.frame(otomatef), nsets = 256,  order.by="freq", nintersects=NA, number.angles = 0, point.size = 2.5, line.size = 0.5, mainbar.y.label = "Genus intersections", sets.x.label = "Genus per sample")
pdf("tomates_upset.pdf", height=10, width=25)
otomatef_ups
dev.off()

olocalesf_ups <- upset(as.data.frame(olocalesf), nsets = 256,  order.by="freq", nintersects=NA, number.angles = 0, point.size = 2.5, line.size = 0.5, mainbar.y.label = "Genus intersections", sets.x.label = "Genus per sample")
pdf("locales_upset.pdf", height=10, width=25)
olocalesf_ups
dev.off()

osuelosf_ups <- upset(as.data.frame(osuelosf), nsets = 256,  order.by="freq", nintersects=NA, number.angles = 0, point.size = 2.5, line.size = 0.5, mainbar.y.label = "Genus intersections", sets.x.label = "Genus per sample")
pdf("suelos_upset.pdf", height=10, width=25)
osuelosf_ups
dev.off()


#Create get_intersect_members function.
get_intersect_members <- function (x, ...){
 require(dplyr)
 require(tibble)
 x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
 n <- names(x)
 x %>% rownames_to_column() -> x
 l <- c(...)
 a <- intersect(names(x), l)
 ar <- vector('list',length(n)+1)
 ar[[1]] <- x
 i=2
 for (item in n) {
   if (item %in% a){
     if (class(x[[item]])=='integer'){
       ar[[i]] <- paste(item, '>= 1')
       i <- i + 1
     }
   } else {
     if (class(x[[item]])=='integer'){
       ar[[i]] <- paste(item, '== 0')
       i <- i + 1
     }
   }
 }
 do.call(filter_, ar) %>% column_to_rownames() -> x
 return(x)
}

#Obtain tomato core microbiome members.
int.tomates<-as.table(otu_table(tomate.phy))
int.tomates <- replace(int.tomates, int.tomates>0, 1)
write.table(int.tomates, "int.tomates.tmp")
tomates.so <- read.table("int.tomates.tmp", header=TRUE, row.names = 1)
get_intersect_members(int.tomates, c(colnames(int.tomates)))->int.tomates.core



#Obtain lists of OTUs from all genera in each subset (tomato, soils or ruderal plants)
#The OTU lists were compared using the web interface  https://bioinformatics.psb.ugent.be/webtools/Venn/ to produce the venn diagrams
write.table(rownames(otu_table(tomate.phy)), "tomate_genus.txt", quote=F)
write.table(rownames(otu_table(locales.phy)), "locales_genus.txt", quote=F)
write.table(rownames(otu_table(suelos.phy)), "suelos_genus.txt", quote=F)

#Obtaining taxonomy of OTUs. sample_type is used as a generalization any type of sample (initial soils, tomatoes, ruderal plants)
for i in `cut -f2` sampletype_genus.txt; do grep -w $i soil_tom_loc_wos.clean.tax >>sampletype_genus_venn_taxonomy.txt; done

```

#Beta diversity analysis


```python
#calculation of weighted Unifrac distances.
UniFrac(stl.ss.phy, weighted=TRUE, normalized=FALSE, parallel=FALSE, fast=TRUE)->stl.weightUF

#Clustering of samples based on weighted UniFrac distances using the hclust method with default parameters.
hclust(stl.weightUF)->stl.weightUF.clus
plot(stl.weightUF.clus)

#Perform anosim on weighted Unifrac distance matrix based on the order of samples in clusters. s=soil sample; t=tomato sample; a=ruderal plant sample.
stl.groups<-as.factor(c("a", "s", "s", "s", "s", "s", "s", "s", "s", "s", "s", "s", "a", "a", "a", "s", "s", "s", "s", "s", "s", "s", "s", "s", "s", "s", "s", "t", "s", "s", "t", "s", "s", "s", "s", "s", "s", "s", "s", "a", "s", "s", "s", "s", "s", "a", "a", "s", "s", "s", "s", "t", "s", "t", "t", "t", "t", "a", "t", "t", "t", "s", "s", "s", "s", "s", "s", "a", "a", "a", "s", "s", "s", "s", "a", "a", "s", "s", "s", "s", "s", "a", "a", "s", "t", "t", "s", "s", "t", "t", "t", "t", "t", "t", "t", "t", "s", "s", "s", "t", "a", "t", "t", "a", "a", "a"))
anosim(stl.weightUF, stl.groups)

#Obtain samplewise Unifrac distances dendrogram. 
samples.wUF.tre<-as.phylo(stl.weightUF.clus)
write.tree(phy=samples.wUF.tre, file="16S_wUF_dend.newick")

#Calculate cophenetic distances of samples based on the 16S rRNA dendrogram.
16s.nwk <- read.tree(“16S_wUF_dend.newick”)
write.table(cophenetic(16s.nwk), file=‘dendrogram_distances.tsv’, quote=FALSE, sep=‘\t’)
```

#Differential OTU abundance


```python
#Differential abundance analysis was performed using DESeq2 package in R.
#Comparison between intial soil and endosphere.
alpha = 0.01
#Subseting the initial soil and endosphere samples and merging into a different phyloseq object.
merge_phyloseq(subset_samples(suelo.phy.ss, muestra=="S0"), subset_samples(suelo.phy.ss, muestra=="EC"))->so.ec.phy
#Converting the phyloseq obaject to deseq format and running DESeq.
so.ec.ds = phyloseq_to_deseq2(so.ec.phy, ~ muestra)
so.ec.ds<-DESeq(so.ec.ds, test="Wald", fitType="local")
so.ec.ds.res<-results(so.ec.ds, cooksCutoff=FALSE)
sigtab.so.ec<-so.ec.ds.res[which(so.ec.ds.res$padj < alpha), ]
sigtab.so.ec<-cbind(as(sigtab.so.ec, "data.frame"), as(tax_table(so.ec.phy)[rownames(sigtab.so.ec), ], "matrix"))
sigtab.so.ec.x=tapply(sigtab.so.ec$log2FoldChange, sigtab.so.ec$Genus, function(x) max(x))
sigtab.so.ec.x=sort(sigtab.so.ec.x, TRUE)
sigtab.so.ec$Genus = factor(as.character(sigtab.so.ec$Genus), levels=names(sigtab.so.ec.x))
ggplot(sigtab.so.ec, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5, alpha=0.85)+theme_light() +theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("EC vs S0; a=0.001; 33 OTUs")

#Comparison between initial soil and rhizosphere.
#Subseting the initial soil and rhizosphere samples and merging into a different phyloseq object.
merge_phyloseq(subset_samples(suelo.phy.ss, muestra=="S0"), subset_samples(suelo.phy.ss, muestra=="RF"))->so.ec.phy
so.rf.ds = phyloseq_to_deseq2(so.rf.phy, ~ muestra)
so.rf.ds<-DESeq(so.rf.ds, test="Wald", fitType="local")
so.rf.ds.res<-results(so.rf.ds, cooksCutoff=FALSE)
sigtab.so.rf<-so.rf.ds.res[which(so.rf.ds.res$padj < alpha), ]
sigtab.so.rf<-cbind(as(sigtab.so.rf, "data.frame"), as(tax_table(so.rf.phy)[rownames(sigtab.so.rf), ], "matrix"))
sigtab.so.rf.x=tapply(sigtab.so.rf$log2FoldChange, sigtab.so.rf$Genus, function(x) max(x))
sigtab.so.rf.x=sort(sigtab.so.rf.x, TRUE)
sigtab.so.rf$Genus = factor(as.character(sigtab.so.rf$Genus), levels=names(sigtab.so.rf.x))
ggplot(sigtab.so.rf, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5, alpha=0.58)+theme_light() +theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("RF vs S0; a=0.01, 53 OTUs")
    
#Comparison between final soil and endosphere.
#Subseting the initial soil and endosphere samples and merging into a different phyloseq object.
merge_phyloseq(subset_samples(suelo.phy.ss, muestra=="SF"), subset_samples(suelo.phy.ss, muestra=="EC"))->so.ec.phy
sf.ec.ds = phyloseq_to_deseq2(sf.ec.phy, ~ muestra)
sf.ec.ds<-DESeq(sf.ec.ds, test="Wald", fitType="local")
sf.ec.ds.res<-results(sf.ec.ds, cooksCutoff=FALSE)
sigtab.sf.ec<-sf.ec.ds.res[which(sf.ec.ds.res$padj < alpha), ]
sigtab.sf.ec
sigtab.sf.ec<-cbind(as(sigtab.sf.ec, "data.frame"), as(tax_table(sf.ec.phy)[rownames(sigtab.sf.ec), ], "matrix"))
sigtab.sf.ec.x=tapply(sigtab.sf.ec$log2FoldChange, sigtab.sf.ec$Genus, function(x) max(x))
sigtab.sf.ec.x=sort(sigtab.sf.ec.x, TRUE)
sigtab.sf.ec$Genus = factor(as.character(sigtab.sf.ec$Genus), levels=names(sigtab.sf.ec.x))
ggplot(sigtab.sf.ec, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=5, alpha=0.85)+theme_light() +theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("EC vs SF; a=0.01; 33 OTUs")
```

#Initial and final soils geochemichal diversity and plant biomass production (Common garden experiment)


```python
#Soils geochemichal data ordination plots were performed with MetaMDS function using the NMDS ordination method. 

#Load the metadata in R
sisf.df<-read.table("metadata_suelo_ord_bgq.txt", header=TRUE, sep="\t", row.names=1)w
#Apply metaMDS function on physicochemical data
sisf.df.nmds<-metaMDS(sisf.df)
#Get the scores of the NMDS analysis and load as a data frame.
sisf.df.nmds.score<-as.data.frame(scores(sisf.df.nmds))
#Add a column with the site names.
sisf.df.nmds.score$site<-rownames(sisf.df.nmds.score)
#Add a time column to the scores data frame.
sisf.df.nmds.score$time<-c("SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SF","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI","SI")
#Plotting the ordination plot with ggplot2, coloring by time.
p.time<-ggplot(sisf.df.nmds.score,aes(x=NMDS1,y=NMDS2,colour=time)) + theme_bw() +geom_point(aes(size=4))
p.time + geom_text_repel(sisf.df.nmds.score, mapping = aes(x=NMDS1,y=NMDS2, label=id), size=3)+ scale_color_brewer(palette="Dark2")

#Barplot for biomass production.

biomasa.df<-read.table("biomasa_stdev.txt", header=T, sep="\t", row.names=1)
p <- ggplot(data = biomasa.df, aes(x = nombre, y = mean, fill=nombre))
p + geom_bar(stat = "identity", position = dodge) + geom_errorbar(limits, position = dodge, width = 0.25) + theme_bw()

```
