
#Paired end assembly and quality filtering


```python
#Obtaining quality filtered, paired and unpaired sequences using Trimmomatic. This is performed for every paired end WSM libraries.

java -jar trimmomatic.jar PE -threads 20 -phred33 -trimlog trim.log SAMPLE_1.fastq.gz SAMPLE_2.fastq.gz SAMPLE_paired_1.fastq.gz SAMPLE_unpaired_1.fastq.gz SAMPLE_paired_2.fastq.gz SAMPLE_unpaired_2.fastq.gz ILLUMINACLIP:/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

#Host sequences filtering.


```python
#Tomato, rice, and arabidopsis genomes were used to filter out host sequences in the metagenomes.
#Tomato was used for tomato rhizosphere samples. Rice and arabidopsis were used for local plants rhizopsheres. 

#Create host genome index databse
bowtie2-build A_thaliana.fna A_thaliana
bowtie2-build O_sativa_genomic.fna O_sativa_genomic
bowtie2-build tomato_genome.fna tomato_genome

#Map tomtato rhizosphere metagenomes against host (tomato) genome with Bowtie2.

/usr/bin/bowtie2 -x tomato_genome -1 SAMPLE_1.fastq.gz -2 SAMPLE_2.fastq.gz -S SAMPLE_hf.sam --quiet -p 4 --very-sensitive
/usr/bin/samtools view -bS SAMPLE_hf.sam > SAMPLE_hf.bam
#filter unmapped reads
/usr/bin/samtools view -b -f 12 -F 256 SAMPLE_hf.bam > SAMPLE_bothEndsUnmapped_hf.bam 
#Sort and split paired sequences into R1 and R2
/usr/bin/samtools sort -n SAMPLE_bothEndsUnmapped_hf.bam SAMPLE_bothEndsUnmapped_sorted_hf 
/usr/bin/bedtools bamtofastq -i SAMPLE_bothEndsUnmapped_sorted_hf.bam -fq SAMPLE_hf_r1.fastq -fq2 SAMPLE_hf_r2.fastq

#Map ruderal plants and soils metagenomes against Oriza sativa.
/usr/bin/bowtie2 -x O_sativa_genomic -1 /SAMPLE_1.fastq.gz -2 SAMPLE_2.fastq.gz -S /SAMPLE_OShf.sam --quiet -p 4 --very-sensitive
/usr/bin/samtools view -bS SAMPLE_OShf.sam > SAMPLE_OShf.bam
#filter unmapped reads
/usr/bin/samtools view -b -f 12 -F 256 SAMPLE_OShf.bam > SAMPLE_bothEndsUnmapped_OShf.bam 
#Sort and split paired sequences into R1 and R2
/usr/bin/samtools sort -n SAMPLE_OShf.bam SAMPLE_bothEndsUnmapped_sorted_OShf 
/usr/bin/bedtools bamtofastq -i SAMPLE_bothEndsUnmapped_sorted_OShf.bam -fq SAMPLE_OShf_r1.fastq -fq2 SAMPLE_OShf_r2.fastq

#Map O. sativa filtered ruderal plants metagenomes againstand Arabidospsis thaliana genomes.
/usr/bin/bowtie2 -x A_thaliana -1 SAMPLE_OShf_r1.fastq -2 SAMPLE_OShf_r2.fastq -S SAMPLE_ATOShf.sam --quiet -p 4 --very-sensitive
/usr/bin/samtools view -bS SAMPLE_ATOShf.sam > SAMPLE_ATOShf.bam
#filter unmapped reads
/usr/bin/samtools view -b -f 12 -F 256 SAMPLE_ATOShf.bam > SAMPLE_bothEndsUnmapped_ATOShf.bam 
#Sort and split paired sequences into R1 and R2
/usr/bin/samtools sort -n SAMPLE_bothEndsUnmapped_ATOShf.bam SAMPLE_bothEndsUnmapped_sorted_ATOShf
/usr/bin/bedtools bamtofastq -i SAMPLE_bothEndsUnmapped_sorted_ATOShf.bam -fq SAMPLE_ATOShf_r1.fastq -fq2 SAMPLE_ATOShf_r2.fastq
```

#Hybrid Assembly using Spades and Velvet.


```python
#A mixed assembly strategy was used. MetaSpades was used to assemble initial contigs. Reads were mapped against resulting contigs, then unmapped reads were assembled using Velvet.

#Spades assembly.
spades.py --meta -1 SAMPLE_r1.fastq.gz  -2 SAMPLE_r2.fastq.gz  -o /home/spades_assembly/

#Mapping paired reads to Spades contigs with BBMap and obtaining unmapped reads for assembly with Velvet
bbwrap.sh ref=SAMPLE_spades_contigs.fasta in=SAMPLE_r1.fastq.gz in2=SAMPLE_r2.fastq.gz out=SAMPLE.sam.gz kfilter=22 subfilter=15 maxindel=80
pileup.sh in=SAMPLE.sam.gz out=SAMPLE_cov.txt

#Separate Spades unmapped reads.
samtools view -u -f12 -F256 SAMPLE.sam.gz | samtools bam2fq - | gzip > SAMPLE_unmapped.pe.fq.gz

#Separation of unmapped reads into forward and reverse.
zcat SAMPLE_unmapped.pe.fq.gz | grep '^@.*/1$' -A 3 --no-group-separator  > SAMPLE_um_r1.fastq
zcat SAMPLE_unmapped.pe.fq.gz | grep '^@.*/2$' -A 3 --no-group-separator  > SAMPLE_um_r2.fastq



#Velvet assembly of reads that did not map to Spades contigs.
velveth SAMPLE_assembly31 31 -shortPaired -fastq SAMPLE_um_r1.fastq SAMPLE_um_r2.fastq
velvetg SAMPLE_assembly31 -exp_cov 2 -ins_length 350

#Mapping paired reads to Velvet contigs.
bbwrap.sh ref=SAMPLE_contigs_velv.fa in=SAMPLE_um_r1.fastq in2=/SAMPLE_um_r2.fastq out=SAMPLE_velvmap.sam.gz kfilter=22 subfilter=15 maxindel=80
pileup.sh in=SAMPLE_velvmap.sam.gz out=SAMPLE_velvcov.txt

#Obtaining final unmapped paired reads.
samtools view -u -f4 SAMPLE_velvmap.sam.gz | samtools bam2fq - | gzip > SAMPLE_unmapped_velv.pe.fq.gz
#Convert final unmapped reads from FASTQ to FASTA using Fastxtoolkit.
fastq_to_fasta -i SAMPLE_unmapped_velv.pe.fq  -o SAMPLE_unmapped_velv.pe.fasta 


#Spades and Velvet contigs were merged. muestras.txt is a file with a list of sample names
for i in `cat muestras`; do echo "cat "$i"_spades_contigs.fasta "$i"_velv_contigs.fasta >"$i"_SVcontigs.fasta"| sh; done

#Merged contigs were length filtered to 100bp.
for i in `cat muestras`; do echo "perl /home/hugo/scripts/remove_small.pl 100 "$i"_SVcontigs.fasta >"$i"_SVc_lf.fasta" | sh; done

#Final contigs were indexed according to sample name. muestras is a text file containing sample names.
for i in `cat muestras`; do echo " perl /home/hugo/scripts/header.fasta.numbers.pl "$i" "$i"_SVc_lf.fasta"| sh; done
```

#ORFs and protein prediction.


```python
#Scripts: prediccion.sh 
#!/bin/bash

SEQS=/path/to/sequences
PROD=/path/to/prodigal

COUNT=0
for FAA in `cat muestras`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr

echo ""$PROD"  -i "$SEQS/$FAA"_SVc_lf.fasta.numbered.fas -a "$SEQS/$FAA"_SVclf.faa -d "$SEQS/$FAA"_SVclf.genes -p meta -o "$SEQS/$FAA".prodigal" >> $*.$COUNT.scr
chmod +x *.scr
done

```

#Metagenome annotation against M5NR databse using Diamond.


```python
#!/bin/bash
#bash nombre_shipt.sh <job_name>

SEQS=/path/to/sequences
SALIDA=path/to/output
BIN2=/path/to/diamond_binaries

COUNT=0
for FAA in `ls *.faa`
do
let COUNT=COUNT+1
echo "#!/bin/bash" >$*.$COUNT.scr
echo "#$ -cwd" >>$*.$COUNT.scr
echo "#$ -j y" >>$*.$COUNT.scr
echo "#$ -S /bin/bash" >>$*.$COUNT.scr
echo "#$ -l h_vmem=12G" >>$*.$COUNT.scr

#diamond alignments against M5NR. 
echo  "$BIN2 blastp -d /databases/nr/nr.dmnd -q $SEQS/$FAA -f 6 -e 1e-10 -k 10 -p 1 -o $SALIDA/$FAA.dout" >>$*.$COUNT.scr
done
```

#Calculation of protein abundance tables


```python
#Gene abundances.
#Make databases for every ORF prediction file of every metagenome.
/usr/bin/bowtie2-build  SAMPLE_SVclf.fna SIN2S0_SVclf

#Align reads to reference ORFs using mapped and unmapped high quality reads.
/usr/bin/bowtie2 -f -x SAMPLE_SVclf -U SAMPLE_rn.fasta -S SAMPLE.sam --quiet -p 20 --very-sensitive

#Obtain only quality alignments to the reference seqeunces and their frequency 
grep -v '^@' SAMPLE.sam | awk '{if($5 == "42") print $3}' | sort | uniq -c > SAMPLE.hits

#Make an OTU table. Every predicted protein is considered an OTU
for FAA in `ls *.faa.numbered.fas | sed -e 's/.faa.numbered.fas//g'`; do grep ">" "$FAA".faa.numbered.fas | sed 's/#/\t/g;s/>//g' | cut -f1 | sed 's/ /\t/g' | awk 'BEGIN{i=0} /.*/{printf "%d\t% s\n",i,$0; i++}' | cut -f1,2 > "$FAA".otu; done

#Get the best hits from diamond output tables based on bitscore values.
for FAA in `ls *.faa.dout | sed -e 's/.faa.dout//g'`; do cat "$FAA".faa.dout |  perl -pe ' $name_col=0; $score_col=11; while(<>) { s/\r?\n//; @F=split /\t/, $_; ($n, $s) = @F[$name_col, $score_col]; if (! exists($max{$n})) { push @names, $n }; if (! exists($max{$n}) || $s > $max{$n}) { $max{$n} = $s; $best{$n} = () }; if ($s == $max{$n}) { $best{$n} .= "$_\n" }; } for $n (@names) { print $best{$n} } ' > best;  perl -e ' $column=0; $unique=0; while(<>) { s/\r?\n//; @F=split /\t/, $_; if (! ($save{$F[$column]}++)) { print "$_\n"; $unique++ } } ' best > "$FAA"_best_uniq; rm best; done

#Simplify output of best hits files.
for i in `cat muestras.txt`; do awk '{print $1"\t"$3"\t"$11"\t"$12"\t"$2}' $i_best_uniq > $i_best.simple.tsv; done

#Using hitter.py script to obtain abundances of each OTU (hits to M5NR) parsing diamond output tables and OTU tables. 
#We used an inhouse pipeline to obtain gene counts available at https://github.com/genomica-fciencias-unam/tabla_genes to obtain 
for N in `cat muestras`; do echo "python3 hitter.py "$N"_best.simple.tsv "$N".hits "$N".tmp.clstrs "$N"" | sh ;done &   #muestras is a text file with a list of sample names.
sed 's/\(.*\)_/\1./' *.hits

#Make a list of all output files from hitter.py
ls *.hout > lista
#Joining hits tables with the script hitter_table.py. available at https://github.com/genomica-fciencias-unam/tabla_genes.
python3 hitter_table.py lista so_rf_rz

#Obtain unique MD5 identifiers for annotation retrieval against Subsystems and Refseq databases using m5nr-tools.pl script.
cut -f1 so_rf_rz.tsv | sed '1d' >ids_md5_all.txt
/databases/m5nr16may17/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source Subsystems --md5 /home/hugo/m5_anot/ids_md5_all.txt  > /home/hugo/m5_anot/Subsystems_all
/databases/m5nr16may17/m5nr-tools.pl --api http://api.metagenomics.anl.gov --option annotation --source Refseq --md5 /home/hugo/m5_anot/ids_md5_all.txt  > /home/hugo/m5_anot/Refseq_all

#Parse the ouput of hitter_table.py with the annotations file (Refseq or SEED Subsystems)
perl -e ' $col1=0; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' so_tom_loc.tsv Subsystems_all >so_tom_loc_subsyst.tsv
perl -e ' $col1=0; $col2=1; ($f1,$f2)=@ARGV; open(F2,$f2); while (<F2>) { s/\r?\n//; @F=split /\t/, $_; $line2{$F[$col2]} .= "$_\n" }; $count2 = $.; open(F1,$f1); while (<F1>) { s/\r?\n//; @F=split /\t/, $_; $x = $line2{$F[$col1]}; if ($x) { $num_changes = ($x =~ s/^/$_\t/gm); print $x; $merged += $num_changes } } warn "\nJoining $f1 column $col1 with $f2 column $col2\n$f1: $. lines\n$f2: $count2 lines\nMerged file: $merged lines\n"; ' so_tom_loc.tsv refseq_all >stl_refseq.tsv

#Remove duplicate md5 identifiers for Refseq tables.
sort -k1,1 stl_refseq.tsv | cut -f1 | uniq >refseq_all_nr_md5s
for i in `cat refseq_all_nr_md5`; do grep "$i" stl_refseq.tsv  | head -n1 >>stl_refseq_clean.tsv; done
#Remove Oranism column from Refseq annotation.
 awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$19}' stl_refseq_clean.tsv > stl_refseq_norganism.tsv

#Move identifiers columns in Subsystems annotated table.
awk -F "\t" '{print $20"-"$19"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26}' so_tom_loc_subsyst_anotado.tsv >so_tom_loc_anotado_completo.tsv

#Remove duplicate md5 identifiers  
cut -f2 so_tom_loc_subsyst_anotado.tsv | sort | uniq >md5_nr
for i in `cat md5_nr` ; do grep $i so_tom_loc_subsistemas_anotado_completo.tsv | head -n1 >> so_tom_loc_subsistemas_anotado_completo_nr.tsv ; done

#Obtain Subsystems ontolgy (tax table) and counts table
awk -F "\t" '{print $1"\t"$21"\t"$22"\t"$23"\t"$24}' so_tom_loc_subsistemas_anotado_completo_nr.tsv | sed "s/\#//g;s/;//g;s/&/_/g;s/'//g" >so_tom_loc_ontology_nr.tsv
cut -f1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19  so_tom_loc_subsistemas_anotado_completo_nr.tsv >so_tom_loc_counts_nr.tsv


```


```python
#Creating a table with hypothetical predicted sequences that did not had matches against M5NR database. 
#Obtain all sequences with hits against M5NR database.
for i in `cat muestras`; do echo "awk '{print $1}' $i_best.simple.tsv > $i_anotados.txt"; done 

#Obtain a list of every sequence.
for i in `cat muestras`; do echo "grep '>' "$i"_SVcontigs.faa | sed 's/>//g' | cut -f1 -d ' ' > "$i"_todos.txt"; done

#Obtain identifiers of non matching sequences.
for i in `cat muestras`; do echo "cat "$i"_anotados.txt "$i"_todos.txt | sort | uniq -c | grep '1 ' | awk '{print \$2}' > "$i"_na.txt"; done

#Extract non matching sequences.
for i in `cat muestras`; do echo "seqtk subseq "$i"_SVcontigs.faa "$i"_na.txt > "$i"_na.faa"; done

#merge non-matching sequences from all samples
cat *_na.faa > todos_na.faa

#Cluster into gene families at 70% protein identity.
cd-hit -i todos_na.faa -o todos70 -c 0.70 -n 4 -aL 0.7 -d 0 -M 3000 -T 2 > todos70.cdhit.out

#Convert cdhit output to otu table
perl -pne 's/\t//g;s/^.*,//g;s/\.\.\..*$//g;s/\n/\t/g;s/\>Cluster\ /\n/g;s/\>//g; eof && do{chomp; print "$_ \n"; exit}' todos70.clstr > todos70.otu
sed -i.bak 's/_/./g;s/RFT./RFT_/g;s/RZPS./RZPS_/g;s/S0./S0_/g' todos70.otu

#Create a list with mapping alignment files.
ls *.hits > para_na_list.txt

#Using hitter_na.py script to obtain abundances of each OTU (hits to M5NR) parsing hits tables and OTU tables. 
#We used an inhouse pipeline to obtain gene counts available at https://github.com/genomica-fciencias-unam/tabla_genes to obtain 
python3 hitter_na.py hits.txt todos70.otu stl_na_todos

#Add two columns of increasing numbers of hipothetical proteins as annotation.
awk '{print $0"\t""hp_"NR-1"\t""hp_"NR-1}' stl_na_todos.tsv >stl_na_todos_numerado.tsv



#Unnanotated sequences tables were merged with Refseq annotation for complete diversity estimation in the metagenomes.
#Remove header from unnnanotated sequences counts table.
sed -1d stl_na_todos_numerado.tsv

#Remove Oranism column from Refseq annotation.
 awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$19}' stl_refseq.tsv > stl_refseq_norganism.tsv

#merge refseq + hypothetical proteins tables. 
cat stl_refseq_norganism.tsv stl_na_todos_numerado.tsv > stl_refseq_na_numerado.tsv

#Refseq + hypothetical proteins counts table.
awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"}' stl_refseq_na_numerado.tsv > stl_refseq_hp_counts.tsv
#Refseq + hypothetical proteins annotation table.
awk -F "\t" '{print $1"\t"$19"\t"$20}' stl_refseq_na_numerado.tsv >stl_refseq_ontol_numerada.tsv

#Edit special characters of ontology table.
sed -i "s/\//-/g" stl_refseq_ontol_numerada.tsv
sed -i 's/"//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\.//g' stl_refseq_ontol_numerada.tsv
sed -i 's/[//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\[//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\]//g' stl_refseq_ontol_numerada.tsv
sed -i 's/;//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\+//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\#//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\(//g' stl_refseq_ontol_numerada.tsv
sed -i 's/(//g' stl_refseq_ontol_numerada.tsv
sed -i 's/)//g' stl_refseq_ontol_numerada.tsv
sed -i 's/://g' stl_refseq_ontol_numerada.tsv
sed -i 's/&//g' stl_refseq_ontol_numerada.tsv
sed -i 's/*//g' stl_refseq_ontol_numerada.tsv
sed -i 's/=//g' stl_refseq_ontol_numerada.tsv
sed -i 's/?//g' stl_refseq_ontol_numerada.tsv
sed -i 's/^//g' stl_refseq_ontol_numerada.tsv
sed -i 's/\^//g' stl_refseq_ontol_numerada.tsv
sed -i 's/,//g' stl_refseq_ontol_numerada.tsv
sed -i "s/'//g" stl_refseq_ontol_numerada.tsv

```


```python
#Loading Refseq data in phyloseq for diversity analysis.
library(phyloseq)
library(vegan)
library(ggplot2)

md5_count<-read.table("stl_refseq_hp_counts.tsv", header=T, row.names=1, sep="\t")
otumat = as(md5_count, "matrix")
OTU = otu_table(otumat, taxa_are_rows=TRUE)
suelos_data <- read.table("metadatos.txt",header=T,row.names=1, sep="\t")
metagenomas <- phyloseq(OTU)
sampledata = sample_data(data.frame(id=suelos_data$id, sitio=suelos_data$Site, edo=suelos_data$estado, muestra=suelos_data$sample_type, sc=suelos_data$soil_class, ph=suelos_data$ph, nt=suelos_data$Nt, pt=suelos_data$Pt, ct=suelos_data$Ct, c_n=suelos_data$C_N, c_p=suelos_data$C_P, n_p=suelos_data$N_P, bm=suelos_data$biomass, alt=suelos_data$altitude, ia=suelos_data$IA, latitud=suelos_data$latitude, longitud=suelos_data$longitude, row.names=sample_names(metagenomas)))
taximat = as.matrix(read.table("stl_refseq_ontol_numerada.tsv", header=T, row.names=1, sep = "\t"))
taxi=tax_table(taximat)
metagenomas = phyloseq(OTU,sampledata,taxi)

#Diversity boxplots
#Using phyloseq we plot different diversity metrics.
boxplot.ggsignif.plot <- plot_richness(metagenomas,  measures=c("Observed", "Shannon", "Chao1"), x="muestra" )  + geom_boxplot() + geom_point(size=4) + geom_signif(comparisons=list(c("RFT", "RZ"), c("RFT", "S0"), c("RZ", "S0")), map_signif_level = TRUE, test=t.test) + theme_light() + theme(axis.text.x=element_text(angle=90, hjust=1)) 
#Obtaining Shannon diversity and plotting with ggplot2 from a shannnon diversity table.
write.table(estimate.richness(metagenomas, "Shannon"), "div_metagenomas_bplot.tsv", quote=F)
read.table("div_metagenomas_bplot.tsv", header=TRUE, sep="\t", row.names=1)->div.metagenoma
ggplot(data=div.metagenoma, aes(x=Variable, y=Value)) + geom_boxplot(aes(fill=Variable)) +geom_point(aes(y=Value), position=position_dodge(width=0.75)) + theme_bw() +ylim(5,13)

#Diversity plots
p <- plot_richness(suelos, x="type", measures=c("Shannon")) + geom_boxplot()+geom_point(size=5, alpha=1)+theme_light() + theme(axis.text.x=element_text(angle=90, hjust=1))
#p$data
a <- p$data

metaa <- read.table("Shannon.csv", sep="\t", header=T, row.names = 1)
head(metaa)
hist(metaa$shannon)
summary(metaa$shannon)
aov.shannon = aov(shannon~type, data=metaa)
summary(aov.shannon)
TukeyHSD(aov.shannon)

#CAP analysis was made using the complete set of predicted hypothetical proteins without matches agains M5NR database and proteins with hits in the database.
metagenomas_cap<-ordinate(metagenomas, "CAP", "bray", ~muestra + ph + nt + pt + ct)
cap_plot <- plot_ordination(metagenomas, metagenomas.cap, color="type", axes =c(1,2)) + geom_point(aes(colour = type), alpha=0.4, size=4) + geom_point(colour = "grey90", size=1.5) + theme_bw()
arrowmat <- vegan::scores(metagenomas.cap, display="bp")
arrowdf <- data.frame(labels=rownames(arrowmat), arrowmat)
# Define the arrow aesthetic mapping
arrow_map <- aes(xend = CAP1,
    yend = CAP2,
    x = 0,
    y = 0,
    shape = NULL,
    color = NULL,
    label = labels)
label_map <- aes(x = 1.3 * CAP1,
    y = 1.3 * CAP2,
    shape = NULL,
    color = NULL,
    label = labels)
arrowhead = arrow(length = unit(0.02, "npc"))
cap_plot  +
  geom_segment(
    mapping = arrow_map,
    size = .5,
    data = arrowdf,
    color = "gray",
    arrow = arrowhead
  ) +
  geom_text(
    mapping = label_map,
    size = 5,
    data = arrowdf,
    show.legend = TRUE
  ) +
geom_point() +geom_text(aes(label=id,hjust=0,vjust=0), size = 2.5)

#Perform ANOVA and ANOSIM analysis 
an <- anova(metagenomas.cap, permutations=9999)
tipo <- get_variable(suelos, "type")
suelos_anosim <- anosim(distance(suelos, "bray"), tipo)
suelos_anosim
```

#Differentially abundant proteins (Refseq + predicted hypothetical proteins).


```python
#Differentially abundant proteins (only refseq annotated) were calculated using Deseq2. 
#Load data into R and create phyloseq objects
otu <- as.matrix(read.table("stl_refseq_hp_counts.tsv", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("stl_refseq_ontol_numerada.tsv", header=T, row.names=1)) 
taxi=tax_table(taximat)
metagenomas = phyloseq(OTU, taxi)
sample_names(suelos)
data= read.table("metadatos_metagenomas.txt", header=T, row.names=1, sep="\t")
sampledata = sample_data(data.frame(id=data$id,site=data$Site,type=data$sample_type, ph=data$ph, pt=data$Pt, ct=data$Ct, nt=data$Nt, row.names=sample_names(suelos)))
metagenomas = phyloseq (OTU, sampledata, taxi)

#Subset for pairwise comparisons and convert phyloseq objects into deseq objects. In each object we exclude one of the sample types to obtain pairwise comparisons (soils vs tomatoes; soils vs ruderals; tomatoes vs ruderals)
so.rft.ds = phyloseq_to_deseq2((subset_samples(metagenomas, muestra!="RZ")), ~muestra)
so.rz.ds = phyloseq_to_deseq2((subset_samples(metagenomas, muestra!="RFT")), ~muestra)
rz.rft.ds = phyloseq_to_deseq2((subset_samples(metagenomas, muestra!="S0")), ~muestra)

#Apply deseq function for each comparison 
so.rft.ds<-DESeq(so.rft.ds, test="Wald", fitType="local")
so.rz.ds<-DESeq(so.rz.ds, test="Wald", fitType="local")
rz.rft.ds<-DESeq(rz.rft.ds, test="Wald", fitType="local")

#Obtain results from deseq objects
so.rft.ds.res<-results(so.rft.ds, cooksCutoff=FALSE)
so.rz.ds.res<-results(so.rz.ds, cooksCutoff=FALSE)
rz.rft.ds.res<-results(rz.rft.ds, cooksCutoff=FALSE)

#Obtain significantly enriched proteins (alpha<0.0001)
sigtab.so.rft<-so.rft.ds.res[which(so.rft.ds.res$padj < 0.001), ]
sigtab.so.rz<-so.rz.ds.res[which(so.rz.ds.res$padj < 0.001), ]
sigtab.rz.rft<-rz.rft.ds.res[which(rz.rft.ds.res$padj < 0.001), ]

#join annotation to results.
sigtab.so.rft<-cbind(as(sigtab.so.rft, "data.frame"), as(tax_table(metagenomas)[rownames(sigtab.so.rft), ], "matrix"))
sigtab.so.rz<-cbind(as(sigtab.so.rz, "data.frame"), as(tax_table(metagenomas)[rownames(sigtab.so.rz), ], "matrix"))
sigtab.rz.rft<-cbind(as(sigtab.rz.rft, "data.frame"), as(tax_table(metagenomas)[rownames(sigtab.rz.rft), ], "matrix"))

#Sort results in the table
sigtab.so.rft.x=tapply(sigtab.so.rft$log2FoldChange, sigtab.so.rft$Protein, function(x) max(x))
sigtab.so.rz.x=tapply(sigtab.so.rz$log2FoldChange, sigtab.so.rz$Protein, function(x) max(x))
sigtab.rz.rft.x=tapply(sigtab.rz.rft$log2FoldChange, sigtab.rz.rft$Protein, function(x) max(x))
sigtab.so.rft.x=sort(sigtab.so.rft.x, TRUE)
sigtab.so.rz.x=sort(sigtab.so.rz.x, TRUE)
sigtab.rz.rft.x=sort(sigtab.rz.rft.x, TRUE)

#Add annotation column
sigtab.so.rft$Protein = factor(as.character(sigtab.so.rft$Protein), levels=names(sigtab.so.rft.x))
sigtab.so.rz$Protein = factor(as.character(sigtab.so.rz$Protein), levels=names(sigtab.so.rz.x))
sigtab.rz.rft$Protein = factor(as.character(sigtab.rz.rft$Protein), levels=names(sigtab.rz.rft.x))

#Plot L2FC 
ggplot(sigtab.so.rft, aes(x=Protein, y=log2FoldChange)) + geom_point(size=5, alpha=0.85)+theme_light() +theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("S0 vs RFT; a=0.07")
ggplot(sigtab.rz.rft, aes(x=Protein, y=log2FoldChange)) + geom_point(size=5, alpha=0.85)+theme_light() +theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("RZ vs RFT; a=0.07"
ggplot(sigtab.so.rz, aes(x=Protein, y=log2FoldChange)) + geom_point(size=5, alpha=0.85)+theme_light() +theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + ggtitle("RZ vs RFT; a=0.07"

write.table(sigtab.rz.rft, "sigtab_rz_rft_metagenoma.txt", sep="\t")
write.table(sigtab.so.rft, "sigtab_so_rft_metagenoma.txt", sep="\t")
write.table(sigtab.so.rft, "sigtab_so_rz_metagenoma.txt", sep="\t")

```

#Soil, tomato and ruderal plants Core metagenomes


```python
#Unique and shared proteins between all samples.
#Load required libraries
library(tidyverse)
library(ggplot2)
library(vegan)
library(phyloseq)
library (ape)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(UpSetR)
library(grid)

#Load data.
otu <- as.matrix(read.table("stl_refseq_hp_counts.tsv", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("stl_refseq_ontol_numerada.tsv", header=T, row.names=1)) 
taxi=tax_table(taximat)
suelos = phyloseq(OTU, taxi)
sample_names(suelos)
data= read.table("metadatos_metagenomas.txt", header=T, row.names=1, sep="\t")
head(data)
sampledata = sample_data(data.frame(id=data$id,site=data$Site,type=data$sample_type, ph=data$ph, pt=data$Pt, ct=data$Ct, nt=data$Nt, row.names=sample_names(suelos)))
suelos = phyloseq (OTU, sampledata, taxi)


#merge sampels based on sample type (soils, tomato rhizospheres, ruderal plants rhizospheres). If a protein is in ony one sample amongst a sample type it is considered present for the sample type (soils, tomtato or ruderals). 
merge2 <- merge_samples(suelos,"type")
x2 <- merge2
merge2 <- as.table(t(otu_table(merge2)))
#Convert counts table to presence absence
merge2 <- replace(merge2, merge2>0, 1)
write.table(merge2,"merge2.tmp")

merge2 <- read.table("merge2.tmp", header=TRUE, row.names = 1)
#Peform venn diagram analysis through upset function.
upset(merge2, sets=c('RFT','RZ','S0'),
     sets.bar.color ="#000000", order.by="freq", empty.intersections="on")


#Obtain proteins in intersection 
intersect0 <- get_intersect_members(merge2, c('RFT','RZ','S0'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "rft-rz-so.csv", sep="\t")

#Obtain unique proteins in tomato rhizospheres  
intersect0 <- get_intersect_members(merge2, c('RFT'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "rft.csv", sep="\t")

#Obtain unique proteins in ruderals rhizospheres  
intersect0 <- get_intersect_members(merge2, c('RZ'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "rz.csv", sep="\t")

#Obtain unique proteins in soils.  
intersect0 <- get_intersect_members(merge2, c('S0'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "so.csv", sep="\t")

#Obtain proteins in soil ruderals intersection 
intersect0 <- get_intersect_members(merge2, c('RZ','S0'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "rz-so.csv", sep="\t")

#Obtain proteins in tomatoes-ruderals intersection 
intersect0 <- get_intersect_members(merge2, c('RFT','RZ'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "rft-rz.csv", sep="\t")

#Obtain proteins in tomatoes-soils intersection 
intersect0 <- get_intersect_members(merge2, c('RFT','S0'))
intersect0 <- (row.names(intersect0)) #identidades del core
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
# tax_table(physeqsubsetOTU)
write.table((tax_table(physeqsubsetOTU)), "rft-so.csv", sep="\t")



#Core metagenomes for each sample type.
#Tomato core metagenome
tomatoes <- subset_samples(suelos, type=="RFT")
tomatoes_clean <- (prune_taxa(taxa_sums(tomatoes) >1, tomatoes))

merge2 <- tomatoes_clean
x2 <- merge2
merge2 <- as.table((otu_table(merge2))) 
merge2 <- replace(merge2, merge2>0, 1)
write.table(merge2,"merge2.tmp")
merge2 <- read.table("merge2.tmp", header=TRUE, row.names = 1)
#Plot Upset diagram and display the number of proteins in each intersection.
upset(merge2, sets=c("AGS1RFT","DGO1RFT","JAL5RFT","SIN2RFT","SLP1RFT"),
      sets.bar.color ="#000000", order.by="freq", empty.intersections="on")

#Obtain the tomato core metagenome proteins annotations(shared proteins between all tomato rhizosphere samples)
intersect0 <- get_intersect_members(merge2, c("AGS1RFT","DGO1RFT","JAL5RFT","SIN2RFT","SLP1RFT"))
intersect0 <- (row.names(intersect0)) 
vec0 <- setNames(nm=c(intersect0))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec0)
physeqsubsetOTU
write.table((tax_table(physeqsubsetOTU)), "Tomato_Core.csv", sep="\t")



#Soils core metagenome
soils <- subset_samples(suelos, type=="S0")
soils_clean <- (prune_taxa(taxa_sums(soils) >1, soils))
         
merge2 <- soils_clean
x2 <- merge2
merge2 <- as.table((otu_table(merge2))) 
merge2 <- replace(merge2, merge2>0, 1)
write.table(merge2,"merge2.tmp")
merge2 <- read.table("merge2.tmp", header=TRUE, row.names = 1)
#Plot Upset diagram and display the number of proteins in each intersection.
upset(merge2, sets=c("AGS1S0","DGO1S0","JAL5S0","SIN2S0","SLP1S0","NAY2S0"),
      sets.bar.color ="#000000", order.by="freq", empty.intersections="on")

#Obtain the soil core metagenome proteins annotations(shared proteins between all  soil samples)
intersect1 <- get_intersect_members(merge2, c("AGS1S0","DGO1S0","JAL5S0","SIN2S0","SLP1S0","NAY2S0"))
intersect1 <- (row.names(intersect1)) 
vec1 <- setNames(nm=c(intersect1))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec1)
write.table((tax_table(physeqsubsetOTU)), "Soil_Core.csv", sep="\t")



#Ruderal plants core metagenome
ruderals <- subset_samples(suelos, type=="RZ")
ruderals_clean <- (prune_taxa(taxa_sums(ruderals) >1, soils))
         
merge2 <- ruderals_clean
x2 <- merge2
merge2 <- as.table((otu_table(merge2)))
merge2 <- replace(merge2, merge2>0, 1)
write.table(merge2,"merge2.tmp")
merge2 <- read.table("merge2.tmp", header=TRUE, row.names = 1)

#Obtain the soil core metagenome proteins annotations(shared proteins between all  soil samples)
upset(merge2, sets=c("AGS1RZPS","DGO1RZPS","JAL5RZPS","SIN2RZPS","SLP1RZPS","NAY2RZPS"),
      sets.bar.color ="#000000", order.by="freq", empty.intersections="on")

#Obtain the ruderals core metagenome proteins annotations(shared proteins between all ruderal samples)
intersect2 <- get_intersect_members(merge2, c("AGS1RZPS","DGO1RZPS","JAL5RZPS","SIN2RZPS","SLP1RZPS","NAY2RZPS"))
intersect2 <- (row.names(intersect2)) 
vec2 <- setNames(nm=c(intersect2))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec2)
physeqsubsetOTU
write.table((tax_table(physeqsubsetOTU)), "RZ_Core.csv", sep="\t")


#Comparison between every sample type Core proteins. Core shared and unique proteins.
otu <- as.matrix(read.table("only_matches.txt", header=T, row.names=1)) 
OTU = otu_table(otu, taxa_are_rows=T)
taximat = as.matrix(read.table("only_matches_ontol.txt", header=T, row.names=1)) 
taxi=tax_table(taximat)
suelos = phyloseq(OTU, taxi)
sample_names(suelos)
data= read.table("metadatos_metagenomas.txt", header=T, row.names=1, sep="\t")#los metadatos deben de estar en el mismo orden en el que estan en la tabla de OTUs
sampledata = sample_data(data.frame(id=data$id,site=data$Site,type=data$sample_type, row.names=sample_names(suelos)))
suelos = phyloseq (OTU, sampledata, taxi)


suelos
merge2 <- merge_samples(suelos,"type")
x2 <- merge2
merge2 <- as.table(t(otu_table(merge2)))
merge2 <- replace(merge2, merge2>0, 1)
write.table(merge2,"merge2.tmp")

##Shared core between tomatoes, ruderals & soils
## in this analysis the shared proteins between the sets are shown, the 308 RFT are exclusive to 
## tomato, the core of the tomato is 2762, the sum of tomato proteins in this plot is larger
## because if a single sample shares genes with soil or rft it would add to the count of tomato

vec <- c(vec0,vec1, vec2)
physeqsubsetOTU <- subset_taxa(suelos, rownames(tax_table(suelos)) %in% vec)
physeqsubsetOTU


merge2 <- merge_samples(physeqsubsetOTU,"type")
x2 <- merge2
merge2 <- as.table(t(otu_table(merge2)))
merge2 <- replace(merge2, merge2>0, 1)
write.table(merge2,"merge2.tmp")

merge2 <- read.table("merge2.tmp", header=TRUE, row.names = 1)

upset(merge2, sets=c('RFT','RZ','S0'),
      sets.bar.color ="#000000", order.by="freq", empty.intersections="on")

#Extract tomato metagenomes protein core members
intersect2 <- get_intersect_members(merge2, c("RFT"))
intersect2 <- (row.names(intersect2)) 
vec2 <- setNames(nm=c(intersect2))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec2)
physeqsubsetOTU

write.table((tax_table(physeqsubsetOTU)), "RFT_ONLY_Core.csv", sep="\t")

#Extract tomato-ruderals  metagenomes shared protein core members
intersect2 <- get_intersect_members(merge2, c("RFT","RZ"))
intersect2 <- (row.names(intersect2)) 
vec2 <- setNames(nm=c(intersect2))
physeqsubsetOTU <- subset_taxa(x2, rownames(tax_table(x2)) %in% vec2)
physeqsubsetOTU
#Obtain tomatoes and ruderal plants Shared core proteins
write.table((tax_table(physeqsubsetOTU)), "RFT_RZ-shared_Core.csv", sep="\t")
```

#Subsystems summary


```python
#Metagenome summary SEED subsystems annotation.

otu <- as.matrix(read.table("so_tom_loc_counts_nr.tsv", header=T, row.names=1)) #tabla de OTUs sin singletons, formato tabular.  eliminados con: http://qiime.org/scripts/filter_otus_from_otu_table.html
OTU = otu_table(otu, taxa_are_rows=T)
head(OTU)
taximat = as.matrix(read.table("so_tom_loc_ontology_nr.tsv", header=T, row.names=1)) #revisar los encabezados
taxi=tax_table(taximat)
suelos = phyloseq(OTU, taxi)
sample_names(suelos)
data= read.table("metadatos_metagenomas.txt", header=T, row.names=1, sep="\t")#los metadatos deben de estar en el mismo orden en el que estan en la tabla de OTUs:  sample_data(camaron)
head(data)
sampledata = sample_data(data.frame(site=data$Site, type=data$sample_type, row.names=sample_names(suelos)))
suelos = phyloseq (OTU, sampledata, taxi)

merge1 = merge_samples(suelos, "type")
merge1
x2 <- tax_glom(merge1, taxrank="level1")
x2

t <-(tax_table(x2)) 
write.table(t, "R_cal1.txt") 
t <- as.data.frame(read.table("R_taxvennupset.txt", header=T, row.names = 1))
# system("paste vennupset.txt taxvennupset.txt >coreupset.txt")
t1 <- as.data.frame(colSums(otu_table(x2))) #Oredr by abundance
head(t1)

t3 <- merge(t,t1, by=0, all=TRUE )

names(t3)[6] <-"abundance" 

t4 <- t3[order(-t3$abundance),]
write.table(t4, "Rcal2.txt")
 t5 <-   as.vector(t4$Row.names)

p <- plot_bar(x2, "level1")
yx <- p$data
    
    yx <- as.data.frame(yx)
head(yx)
yx <- yx[order(-yx$Abundance),]

head(yx)
levels1=c("RFT","S0","RZ")
levels1
levels2 = rev(t5)

x3 = transform_sample_counts(x2, function(x) x/sum(x))

yx$orden <- factor(yx$Sample, levels = levels1)
yx$abun <- factor(yx$OTU, levels= levels2)
head(yx$abun)
                                   
    
ylabvec = as(tax_table(x2)[,"level1"], "character")
names(ylabvec) <- taxa_names(x2)
ylabvec[is.na(ylabvec)] <- ""
    
plot_heatmap(x3, "NULL", sample.order=levels1, taxa.order=rev(t5)) + scale_y_discrete(labels=ylabvec)  + theme_minimal () +
   scale_fill_gradient(low="#ffffff", high="#000000", na.value = "white", trans = "log10")  +
   theme(axis.text.x = element_text( angle = 90, hjust = 1),
         axis.text.y = element_text(size = 10))
                             
                             
svg("heatmap_l1_merged.svg")
plot_heatmap(x3, "NULL", sample.order=levels1, taxa.order=rev(t5)) + scale_y_discrete(labels=ylabvec)  + theme_minimal () +
   scale_fill_gradient(low="#ffffff", high="#000000", na.value = "white", trans = "log10")  +
   theme(axis.text.x = element_text( angle = 90, hjust = 1),
         axis.text.y = element_text(size = 10))
dev.off()



ta <- (tax_table(x3)[,1])
ot <- t(otu_table(x3))
ttt <- (psmelt(x3))[,2]
#join tax and OTU table with annotation.
cphylum_stl_rel <- data.frame(cbind(ta,ot))
cphylum_stl_rel <- cphylum_stl_rel[, colSums(is.na(cphylum_stl_rel)) != nrow(cphylum_stl_rel)]   
tcphylum_stl_rel <- data.frame(t(cphylum_stl_rel))
colnames(tcphylum_stl_rel) <- NULL

#Esport table to load without tax_glom
write.table(tcphylum_stl_rel, "tcphylum_stl_rel.tsv", quote = TRUE, sep = "\t")

#Load table
tcphylum_stl_rel <- read.table("tcphylum_stl_rel.tsv", header = TRUE, sep = "\t")
                             
#Add metadata
microbiome <- ttt
# microbiome
mtcphylum_stl_rel <- cbind(microbiome, tcphylum_stl_rel, row.names = NULL)
mtcphylum_stl_rel

#Testing ANOVA and TUkey for each categorycategoría

melt_phylum_stl_rel <- melt(mtcphylum_stl_rel)
# head(melt_phylum_stl_rel)
tested_cat <- filter(melt_phylum_stl_rel, variable == "Potassium_metabolism")
tested_cat
# Perform ANOVA
abac <- aov(value~microbiome, tested_cat)
abac

summary(abac)

TukeyHSD(abac)


categ <- (melt_phylum_stl_rel)[,3]
categ <- categ[!duplicated(categ)]
n <- length(categ)
n
categ <- as.vector(categ)
categ
categ[27] 

anovota <- lapply(categ, function(x){
    tested_cat <- filter(melt_phylum_stl_rel, variable == as.name(x))
     aov(value~microbiome, tested_cat)
})
                  

lapply(anovota, summary)
lapply(anovota, TukeyHSD)

```
