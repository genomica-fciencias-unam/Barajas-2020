

```R
#Cargar bibliotecas
library(phyloseq)
library(ggplot2)
library(ape)
library(reshape2)
library(tidyr)
library(dplyr)
library(agricolae)
```


```R
#Cargar tablas de OTUs
otu<-as.matrix(read.table("soil_tom_loc_wos.clean.otu1.txt", header=T, row.names=1))
OTU = otu_table(otu, taxa_are_rows=T)

#Carga de taxonomía
taximat = as.matrix(read.table("soil_tom_loc_wos.clean.tax", header=T, row.names=1))
taxi=tax_table(taximat)

#Objeto phyloseq
stl.phy<-phyloseq(OTU, taxi)
suelos_data2=read.table("metadatos2_stl.txt", header=TRUE, row.names=1, sep="\t")

#Metadatos
sampledata2=sample_data(data.frame(id=suelos_data2$id.1, sitio=suelos_data2$site, estado=suelos_data2$estado, mb=suelos_data2$microbiome, muestra=suelos_data2$sample_type, samp=suelos_data2$muestra, sc=suelos_data2$soil_class, tiempo=suelos_data2$tiempo, ph=suelos_data2$ph, nt=suelos_data2$Nt, pt=suelos_data2$Pt, ct=suelos_data2$Ct, c_n=suelos_data2$C_N, c_p=suelos_data2$C_P, n_p=suelos_data2$N_P, bm=suelos_data2$biomass, alt=suelos_data2$altitude, ia=suelos_data2$IA, latitud=suelos_data2$latitude, longitud=suelos_data2$longitude, pol=suelos_data2$polimerase, conc_amp=suelos_data2$conc_amp, abs_amp=suelos_data2$abs_amp, conc_dna=suelos_data2$conc_dna, abs_dna=suelos_data2$abs_dna, row.names=sample_names(stl.phy)))

#Cargar arbol
read.tree("stl_wos_otus.bacteria.tre")-> stl_wos_tree

#Objeto phyloseq con metadatos
stl.phy<-phyloseq(OTU, taxi, sampledata2, stl_wos_tree)
```


```R
#Aglomerar phylum
Phylum_stl <- tax_glom(stl.phy, taxrank="Phylum")
```


```R
#Obtener la abundancia relativa de los Phylum
Phylum_stl_rel <- transform_sample_counts(Phylum_stl, function(x) x / sum(x))
```


```R
#Unir la tabla de taxonomía con los otus y con los metadatos
cphylum_stl_rel <- data.frame(cbind(tax_table(Phylum_stl_rel), otu_table(Phylum_stl_rel)), row.names=NULL)
cphylum_stl_rel <- cphylum_stl_rel[, colSums(is.na(cphylum_stl_rel)) != nrow(cphylum_stl_rel)]

tcphylum_stl_rel <- data.frame(t(cphylum_stl_rel))
tcphylum_stl_rel <- tcphylum_stl_rel[-1,]
colnames(tcphylum_stl_rel) <- NULL

#Exportar tabla para cargar directamente sin tax_glom
write.table(tcphylum_stl_rel, "tcphylum_stl_rel.tsv", quote = TRUE, sep = "\t")

#Cargar tabla
tcphylum_stl_rel <- read.table("tcphylum_stl_rel.tsv", header = TRUE, sep = "\t")


#Agregar los metadatos
microbiome <- data.frame(suelos_data2$microbiome)
mtcphylum_stl_rel <- cbind(microbiome, tcphylum_stl_rel, row.names = NULL) %>% drop_na()

```


```R
#Hacer gráfica
melt_phylum_stl_rel <- melt(mtcphylum_stl_rel)
colnames(melt_phylum_stl_rel) <- c("Microbiome_type", "Sample", "Phylum", "Relative_abundance")

#Cambiar los nombres de las variables que se quiere agregar a las gráficas
melt_phylum_stl_rel$Microbiome_type <- sub("Plant_local", "Ruderal plants", melt_phylum_stl_rel$Microbiome_type)
melt_phylum_stl_rel$Microbiome_type <- sub("Plant", "S. lycopersicum", melt_phylum_stl_rel$Microbiome_type)
melt_phylum_stl_rel$Microbiome_type <- sub("Soil", "Soils", melt_phylum_stl_rel$Microbiome_type)

gsub(pattern, replacement, x, ignore.case = FALSE, perl = FALSE,
     fixed = FALSE, useBytes = FALSE)

plot_phylum_stl_rel <- ggplot(melt_phylum_stl_rel, aes(Phylum, Relative_abundance, fill=Microbiome_type)) + 
                            geom_boxplot() + 
                            scale_fill_manual(values = c("#4883b2ff", "#fd0000ff", "#007f00ff")) +
                            theme_light(base_size = 22) +
                            theme(axis.text.y = element_text(size = 20), 
                                  axis.text.x = element_text(angle = 90, size = 20)) 
                                    

ggsave("plot_phylum_stl_rel.pdf", width=80, height=40, units="cm")
```


```R
#Gráfica con Phyla de mayor abundancia
melt_phylum_stl_rel <- melt(mtcphylum_stl_rel)
colnames(melt_phylum_stl_rel) <- c("Microbiome_type", "Sample", "Phylum", "Relative_abundance")
fmelt_phylum_stl_rel <- filter(melt_phylum_stl_rel, Phylum == "p_Proteobacteria" |
                               Phylum == "p_Actinobacteria" | Phylum == "p_Acidobacteria" |
                               Phylum == "p_Planctomycetes" | Phylum == "p_Bacteroidetes" |
                               Phylum == "p_Chloroflexi" | Phylum == "p_Cyanobacteria" |
                               Phylum == "p_Verrucomicrobia" | 
                               Phylum == "p_TM7" | Phylum == "p_Armatimonadetes" |
                               Phylum == "p_Firmicutes" | Phylum == "p_Gemmatimonadetes")


#Cambiar los nombres de las variables que se quiere agregar a las gráficas
fmelt_phylum_stl_rel$Microbiome_type <- sub("Plant_local", "Ruderal plants", fmelt_phylum_stl_rel$Microbiome_type)
fmelt_phylum_stl_rel$Microbiome_type <- sub("Plant", "S. lycopersicum", fmelt_phylum_stl_rel$Microbiome_type)
fmelt_phylum_stl_rel$Microbiome_type <- sub("Soil", "Soils", fmelt_phylum_stl_rel$Microbiome_type)

#Ordenar como se quieren graficar
fmelt_phylum_stl_rel$Phylum_o = factor(fmelt_phylum_stl_rel$Phylum, levels=c('p_Proteobacteria', 
                                        'p_Verrucomicrobia', 'p_Firmicutes', 'p_Actinobacteria', 
                                        'p_Armatimonadetes', 'p_Chloroflexi', 'p_Acidobacteria', 
                                        'p_Planctomycetes', 'p_Bacteroidetes','p_Gemmatimonadetes',  'p_TM7',
                                        'p_Cyanobacteria'))
                                        
#Crear la gráfica
plot_fphylum_stl_rel <- ggplot(fmelt_phylum_stl_rel, aes(Phylum, Relative_abundance, fill=Microbiome_type)) + 
                            geom_boxplot() + facet_wrap(~ Phylum_o, scales="free") + 
                            scale_fill_manual(values = c("#009E73", "#E69F00", "#CC79A7"))+
                            theme_bw(base_size = 30) + theme(strip.background = element_blank(), 
                                                             strip.text = element_blank()) +
                            theme(axis.text.y = element_text(size = 22), 
                                  axis.text.x = element_text(size = 22)) 
                                    

ggsave("plot_fphylum_stl_rel.pdf", width=80, height=40, units="cm")
```


```R
#ANOVA para cada grupo
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
```


```R
#Luego se hace sumary para cada anova y se obtiene lo siguiente
> summary(apr)
                 Df Sum Sq Mean Sq F value  Pr(>F)    
Microbiome_type   2  0.481   0.241      48 1.8e-15 ***
Residuals       103  0.516   0.005                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(ave)
                 Df  Sum Sq  Mean Sq F value Pr(>F)  
Microbiome_type   2 0.00091 0.000456     3.9  0.023 *
Residuals       103 0.01207 0.000117                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(afi)
                 Df  Sum Sq  Mean Sq F value Pr(>F)   
Microbiome_type   2 0.00079 0.000394    5.14 0.0074 **
Residuals       103 0.00789 0.000077                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(aac)
                 Df Sum Sq Mean Sq F value Pr(>F)    
Microbiome_type   2  0.611  0.3054    64.1 <2e-16 ***
Residuals       103  0.491  0.0048                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

> summary(aar)
                 Df  Sum Sq  Mean Sq F value Pr(>F)  
Microbiome_type   2 0.00028 1.39e-04    4.06   0.02 *
Residuals       103 0.00353 3.43e-05                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(ach)
                 Df Sum Sq Mean Sq F value Pr(>F)   
Microbiome_type   2 0.0295 0.01476    7.23 0.0012 **
Residuals       103 0.2104 0.00204                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(aci)
                 Df Sum Sq Mean Sq F value Pr(>F)   
Microbiome_type   2 0.0196 0.00980    5.11 0.0076 **
Residuals       103 0.1975 0.00192                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(apl)
                 Df Sum Sq Mean Sq F value  Pr(>F)    
Microbiome_type   2 0.0734  0.0367    23.1 5.3e-09 ***
Residuals       103 0.1638  0.0016                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(aba)
                 Df Sum Sq Mean Sq F value  Pr(>F)    
Microbiome_type   2 0.0917  0.0459    48.6 1.3e-15 ***
Residuals       103 0.0971  0.0009                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(age)
                 Df Sum Sq  Mean Sq F value Pr(>F)
Microbiome_type   2 0.0016 0.000794    1.72   0.18
Residuals       103 0.0475 0.000461               
> summary(atm)
                 Df Sum Sq  Mean Sq F value Pr(>F)  
Microbiome_type   2 0.0026 0.001292    4.15  0.019 *
Residuals       103 0.0321 0.000312                 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(acy)
                 Df Sum Sq  Mean Sq F value Pr(>F)
Microbiome_type   2 0.0018 0.000921    1.69   0.19
Residuals       103 0.0560 0.000543 
```


```R
#Prueba de Tukey para cada comparación
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
```


```R
#Resultados prueba de Tukey
> TukeyHSD(apr)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Proteobacteria"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants  0.09078108698433995371602
Soils-Ruderal plants           -0.06816022387657327818999
Soils-S. lycopersicum          -0.15894131086091323190601
                                                      lwr
S. lycopersicum-Ruderal plants  0.04680584294240299619849
Soils-Ruderal plants           -0.10879746094195177907871
Soils-S. lycopersicum          -0.19750945585733131748540
                                                      upr
S. lycopersicum-Ruderal plants  0.13475633102627690429465
Soils-Ruderal plants           -0.02752298681119477036239
Soils-S. lycopersicum          -0.12037316586449514632662
                                                     p adj
S. lycopersicum-Ruderal plants 1.023261034138744918209e-05
Soils-Ruderal plants           3.637819189995505908541e-04
Soils-S. lycopersicum          1.854072450000000121804e-14

> TukeyHSD(ave)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Verrucomicrobia"))

$Microbiome_type
                                                       diff
S. lycopersicum-Ruderal plants  0.0006576295084960472048241
Soils-Ruderal plants           -0.0055298446421053158716941
Soils-S. lycopersicum          -0.0061874741506013630765182
                                                       lwr
S. lycopersicum-Ruderal plants -0.006068866872567787139336
Soils-Ruderal plants           -0.011745756208076963203468
Soils-S. lycopersicum          -0.012086895368033329783986
                                                        upr
S. lycopersicum-Ruderal plants  0.0073841258895598815489847
Soils-Ruderal plants            0.0006860669238663314600801
Soils-S. lycopersicum          -0.0002880529331693963690508
                                                   p adj
S. lycopersicum-Ruderal plants 0.97064952515674840505255
Soils-Ruderal plants           0.09164308562795142609758
Soils-S. lycopersicum          0.03751187681850987054588

> TukeyHSD(afi)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Firmicutes"))

$Microbiome_type
                                                       diff
S. lycopersicum-Ruderal plants -0.0008210727180562550284182
Soils-Ruderal plants            0.0050098274855546491221570
Soils-S. lycopersicum           0.0058309002036109041505751
                                                        lwr
S. lycopersicum-Ruderal plants -6.261211185184806758686e-03
Soils-Ruderal plants           -1.736909683071503895935e-05
Soils-S. lycopersicum           1.059669167332521165503e-03
                                                      upr
S. lycopersicum-Ruderal plants 0.004619065749072296701849
Soils-Ruderal plants           0.010037024067940013283273
Soils-S. lycopersicum          0.010602131239889286268285
                                                   p adj
S. lycopersicum-Ruderal plants 0.93150378513183551554278
Soils-Ruderal plants           0.05100453666737447822044
Soils-S. lycopersicum          0.01235173250893206731149

> TukeyHSD(aac)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Actinobacteria"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants -0.19737951372513898307126
Soils-Ruderal plants           -0.06781879434934484640607
Soils-S. lycopersicum           0.12956071937579413666519
                                                      lwr
S. lycopersicum-Ruderal plants -0.24028651640121712063092
Soils-Ruderal plants           -0.10746887651488812309708
Soils-S. lycopersicum           0.09192946714227520632612
                                                      upr
S. lycopersicum-Ruderal plants -0.15447251104906084551160
Soils-Ruderal plants           -0.02816871218380157665395
Soils-S. lycopersicum           0.16719197160931306700427
                                                     p adj
S. lycopersicum-Ruderal plants 6.883382799999999734377e-15
Soils-Ruderal plants           2.729803361829485908174e-04
Soils-S. lycopersicum          2.288391698399999927649e-12

> TukeyHSD(aar)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Armatimonadetes"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants  0.004271724460469040682287
Soils-Ruderal plants            0.001663923371640654326475
Soils-S. lycopersicum          -0.002607801088828386355811
                                                        lwr
S. lycopersicum-Ruderal plants  0.0006341514553898780492958
Soils-Ruderal plants           -0.0016975341886415339261174
Soils-S. lycopersicum          -0.0057981061377590471447863
                                                       upr
S. lycopersicum-Ruderal plants 0.0079092974655482033152776
Soils-Ruderal plants           0.0050253809319228421453873
Soils-S. lycopersicum          0.0005825039601022744331638
                                                   p adj
S. lycopersicum-Ruderal plants 0.01701035225403535466882
Soils-Ruderal plants           0.46945926355297828980184
Soils-S. lycopersicum          0.13168979824315629745257

> TukeyHSD(ach)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Chloroflexi"))

$Microbiome_type
                                                    diff
S. lycopersicum-Ruderal plants 0.01027678193434432285525
Soils-Ruderal plants           0.03827738209022168991424
Soils-S. lycopersicum          0.02800060015587736705900
                                                       lwr
S. lycopersicum-Ruderal plants -0.017808186204712920380766
Soils-Ruderal plants            0.012324245658458490154663
Soils-S. lycopersicum           0.003368897774101757558451
                                                     upr
S. lycopersicum-Ruderal plants 0.03836175007340156262181
Soils-Ruderal plants           0.06423051852198488620438
Soils-S. lycopersicum          0.05263230253765297655955
                                                    p adj
S. lycopersicum-Ruderal plants 0.660258560998454924195755
Soils-Ruderal plants           0.001928839701890949992827
Soils-S. lycopersicum          0.021719861422787811733315

> TukeyHSD(aci)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Acidobacteria"))

$Microbiome_type
                                                    diff
S. lycopersicum-Ruderal plants 0.01031398398596121424031
Soils-Ruderal plants           0.03185712443515881231448
Soils-S. lycopersicum          0.02154314044919759807417
                                                       lwr
S. lycopersicum-Ruderal plants -0.016898268118100124440817
Soils-Ruderal plants            0.006710459216405024324148
Soils-S. lycopersicum          -0.002323153136976507937561
                                                     upr
S. lycopersicum-Ruderal plants 0.03752623609002254945199
Soils-Ruderal plants           0.05700378965391260377427
Soils-S. lycopersicum          0.04540943403537170408590
                                                    p adj
S. lycopersicum-Ruderal plants 0.640694082350955040894291
Soils-Ruderal plants           0.009067444849775840864936
Soils-S. lycopersicum          0.085591175392197738069910

> TukeyHSD(apl)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Planctomycetes"))

$Microbiome_type
                                                    diff
S. lycopersicum-Ruderal plants 0.00251495265626312802798
Soils-Ruderal plants           0.05428647530048966163285
Soils-S. lycopersicum          0.05177152264422653360487
                                                      lwr
S. lycopersicum-Ruderal plants -0.02226838449960967178387
Soils-Ruderal plants            0.03138435459107503239240
Soils-S. lycopersicum           0.03003548997977388904101
                                                     upr
S. lycopersicum-Ruderal plants 0.02729828981213592783983
Soils-Ruderal plants           0.07718859600990429781220
Soils-S. lycopersicum          0.07350755530867918163818
                                                     p adj
S. lycopersicum-Ruderal plants 9.684165710458656484150e-01
Soils-Ruderal plants           4.543009758384287085461e-07
Soils-S. lycopersicum          4.025835741439109938256e-07

> TukeyHSD(aba)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Bacteroidetes"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants  0.05424779061476089209171
Soils-Ruderal plants           -0.01370182907534609612066
Soils-S. lycopersicum          -0.06794961969010698821236
                                                      lwr
S. lycopersicum-Ruderal plants  0.03516342947445400685114
Soils-Ruderal plants           -0.03133756310698351810196
Soils-S. lycopersicum          -0.08468740986408375415184
                                                      upr
S. lycopersicum-Ruderal plants  0.07333215175506777039338
Soils-Ruderal plants            0.00393390495629132586064
Soils-S. lycopersicum          -0.05121182951613022227288
                                                     p adj
S. lycopersicum-Ruderal plants 2.549380928584799795605e-09
Soils-Ruderal plants           1.594317049816702525078e-01
Soils-S. lycopersicum          1.942890290000000088305e-14

> TukeyHSD(age)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Gemmatimonadetes"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants -0.004589798803163333551502
Soils-Ruderal plants            0.004503298759339843804206
Soils-S. lycopersicum           0.009093097562503177355708
                                                       lwr
S. lycopersicum-Ruderal plants -0.017933413016534421879911
Soils-Ruderal plants           -0.007827448366943323687650
Soils-S. lycopersicum          -0.002609815268788376677045
                                                      upr
S. lycopersicum-Ruderal plants 0.008753815410207754776906
Soils-Ruderal plants           0.016834045885623011296062
Soils-S. lycopersicum          0.020796010393794729653738
                                                  p adj
S. lycopersicum-Ruderal plants 0.6927697073340020050836
Soils-Ruderal plants           0.6613140305723815837169
Soils-S. lycopersicum          0.1593884426287284217238

> TukeyHSD(atm)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_TM7"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants  0.012995191820297756490898
Soils-Ruderal plants            0.004991639647612068977911
Soils-S. lycopersicum          -0.008003552172685687512987
                                                       lwr
S. lycopersicum-Ruderal plants  0.002024595092746687627350
Soils-Ruderal plants           -0.005146217450827397787982
Soils-S. lycopersicum          -0.017625228515900877301892
                                                      upr
S. lycopersicum-Ruderal plants 0.023965788547848823619724
Soils-Ruderal plants           0.015129496746051535743804
Soils-S. lycopersicum          0.001618124170529500541194
                                                   p adj
S. lycopersicum-Ruderal plants 0.01589890661870774213327
Soils-Ruderal plants           0.47319069196295726609236
Soils-S. lycopersicum          0.12276672971707792980656

> TukeyHSD(acy)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Cyanobacteria"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants  0.011190279617914310889004
Soils-Ruderal plants            0.006599922080886841749270
Soils-S. lycopersicum          -0.004590357537027469139734
                                                       lwr
S. lycopersicum-Ruderal plants -0.003294230998081194503047
Soils-Ruderal plants           -0.006785119992259433518988
Soils-S. lycopersicum          -0.017293884664861588346119
                                                     upr
S. lycopersicum-Ruderal plants 0.02567479023390981801578
Soils-Ruderal plants           0.01998496415403311701753
Soils-S. lycopersicum          0.00811316959080665006665
                                                  p adj
S. lycopersicum-Ruderal plants 0.1626683616813151544989
Soils-Ruderal plants           0.4721862655968949162855
Soils-S. lycopersicum          0.6670633530799907617848
> TukeyHSD(apr)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Proteobacteria"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants  0.09078108698433995371602
Soils-Ruderal plants           -0.06816022387657327818999
Soils-S. lycopersicum          -0.15894131086091323190601
                                                      lwr
S. lycopersicum-Ruderal plants  0.04680584294240299619849
Soils-Ruderal plants           -0.10879746094195177907871
Soils-S. lycopersicum          -0.19750945585733131748540
                                                      upr
S. lycopersicum-Ruderal plants  0.13475633102627690429465
Soils-Ruderal plants           -0.02752298681119477036239
Soils-S. lycopersicum          -0.12037316586449514632662
                                                     p adj
S. lycopersicum-Ruderal plants 1.023261034138744918209e-05
Soils-Ruderal plants           3.637819189995505908541e-04
Soils-S. lycopersicum          1.854072450000000121804e-14

> TukeyHSD(ave)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Verrucomicrobia"))

$Microbiome_type
                                                       diff
S. lycopersicum-Ruderal plants  0.0006576295084960472048241
Soils-Ruderal plants           -0.0055298446421053158716941
Soils-S. lycopersicum          -0.0061874741506013630765182
                                                       lwr
S. lycopersicum-Ruderal plants -0.006068866872567787139336
Soils-Ruderal plants           -0.011745756208076963203468
Soils-S. lycopersicum          -0.012086895368033329783986
                                                        upr
S. lycopersicum-Ruderal plants  0.0073841258895598815489847
Soils-Ruderal plants            0.0006860669238663314600801
Soils-S. lycopersicum          -0.0002880529331693963690508
                                                   p adj
S. lycopersicum-Ruderal plants 0.97064952515674840505255
Soils-Ruderal plants           0.09164308562795142609758
Soils-S. lycopersicum          0.03751187681850987054588

> TukeyHSD(afi)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Firmicutes"))

$Microbiome_type
                                                       diff
S. lycopersicum-Ruderal plants -0.0008210727180562550284182
Soils-Ruderal plants            0.0050098274855546491221570
Soils-S. lycopersicum           0.0058309002036109041505751
                                                        lwr
S. lycopersicum-Ruderal plants -6.261211185184806758686e-03
Soils-Ruderal plants           -1.736909683071503895935e-05
Soils-S. lycopersicum           1.059669167332521165503e-03
                                                      upr
S. lycopersicum-Ruderal plants 0.004619065749072296701849
Soils-Ruderal plants           0.010037024067940013283273
Soils-S. lycopersicum          0.010602131239889286268285
                                                   p adj
S. lycopersicum-Ruderal plants 0.93150378513183551554278
Soils-Ruderal plants           0.05100453666737447822044
Soils-S. lycopersicum          0.01235173250893206731149

> TukeyHSD(aac)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Actinobacteria"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants -0.19737951372513898307126
Soils-Ruderal plants           -0.06781879434934484640607
Soils-S. lycopersicum           0.12956071937579413666519
                                                      lwr
S. lycopersicum-Ruderal plants -0.24028651640121712063092
Soils-Ruderal plants           -0.10746887651488812309708
Soils-S. lycopersicum           0.09192946714227520632612
                                                      upr
S. lycopersicum-Ruderal plants -0.15447251104906084551160
Soils-Ruderal plants           -0.02816871218380157665395
Soils-S. lycopersicum           0.16719197160931306700427
                                                     p adj
S. lycopersicum-Ruderal plants 6.883382799999999734377e-15
Soils-Ruderal plants           2.729803361829485908174e-04
Soils-S. lycopersicum          2.288391698399999927649e-12

> TukeyHSD(aar)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Armatimonadetes"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants  0.004271724460469040682287
Soils-Ruderal plants            0.001663923371640654326475
Soils-S. lycopersicum          -0.002607801088828386355811
                                                        lwr
S. lycopersicum-Ruderal plants  0.0006341514553898780492958
Soils-Ruderal plants           -0.0016975341886415339261174
Soils-S. lycopersicum          -0.0057981061377590471447863
                                                       upr
S. lycopersicum-Ruderal plants 0.0079092974655482033152776
Soils-Ruderal plants           0.0050253809319228421453873
Soils-S. lycopersicum          0.0005825039601022744331638
                                                   p adj
S. lycopersicum-Ruderal plants 0.01701035225403535466882
Soils-Ruderal plants           0.46945926355297828980184
Soils-S. lycopersicum          0.13168979824315629745257

> TukeyHSD(ach)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Chloroflexi"))

$Microbiome_type
                                                    diff
S. lycopersicum-Ruderal plants 0.01027678193434432285525
Soils-Ruderal plants           0.03827738209022168991424
Soils-S. lycopersicum          0.02800060015587736705900
                                                       lwr
S. lycopersicum-Ruderal plants -0.017808186204712920380766
Soils-Ruderal plants            0.012324245658458490154663
Soils-S. lycopersicum           0.003368897774101757558451
                                                     upr
S. lycopersicum-Ruderal plants 0.03836175007340156262181
Soils-Ruderal plants           0.06423051852198488620438
Soils-S. lycopersicum          0.05263230253765297655955
                                                    p adj
S. lycopersicum-Ruderal plants 0.660258560998454924195755
Soils-Ruderal plants           0.001928839701890949992827
Soils-S. lycopersicum          0.021719861422787811733315

> TukeyHSD(aci)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Acidobacteria"))

$Microbiome_type
                                                    diff
S. lycopersicum-Ruderal plants 0.01031398398596121424031
Soils-Ruderal plants           0.03185712443515881231448
Soils-S. lycopersicum          0.02154314044919759807417
                                                       lwr
S. lycopersicum-Ruderal plants -0.016898268118100124440817
Soils-Ruderal plants            0.006710459216405024324148
Soils-S. lycopersicum          -0.002323153136976507937561
                                                     upr
S. lycopersicum-Ruderal plants 0.03752623609002254945199
Soils-Ruderal plants           0.05700378965391260377427
Soils-S. lycopersicum          0.04540943403537170408590
                                                    p adj
S. lycopersicum-Ruderal plants 0.640694082350955040894291
Soils-Ruderal plants           0.009067444849775840864936
Soils-S. lycopersicum          0.085591175392197738069910

> TukeyHSD(apl)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Planctomycetes"))

$Microbiome_type
                                                    diff
S. lycopersicum-Ruderal plants 0.00251495265626312802798
Soils-Ruderal plants           0.05428647530048966163285
Soils-S. lycopersicum          0.05177152264422653360487
                                                      lwr
S. lycopersicum-Ruderal plants -0.02226838449960967178387
Soils-Ruderal plants            0.03138435459107503239240
Soils-S. lycopersicum           0.03003548997977388904101
                                                     upr
S. lycopersicum-Ruderal plants 0.02729828981213592783983
Soils-Ruderal plants           0.07718859600990429781220
Soils-S. lycopersicum          0.07350755530867918163818
                                                     p adj
S. lycopersicum-Ruderal plants 9.684165710458656484150e-01
Soils-Ruderal plants           4.543009758384287085461e-07
Soils-S. lycopersicum          4.025835741439109938256e-07

> TukeyHSD(aba)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Bacteroidetes"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants  0.05424779061476089209171
Soils-Ruderal plants           -0.01370182907534609612066
Soils-S. lycopersicum          -0.06794961969010698821236
                                                      lwr
S. lycopersicum-Ruderal plants  0.03516342947445400685114
Soils-Ruderal plants           -0.03133756310698351810196
Soils-S. lycopersicum          -0.08468740986408375415184
                                                      upr
S. lycopersicum-Ruderal plants  0.07333215175506777039338
Soils-Ruderal plants            0.00393390495629132586064
Soils-S. lycopersicum          -0.05121182951613022227288
                                                     p adj
S. lycopersicum-Ruderal plants 2.549380928584799795605e-09
Soils-Ruderal plants           1.594317049816702525078e-01
Soils-S. lycopersicum          1.942890290000000088305e-14

> TukeyHSD(age)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Gemmatimonadetes"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants -0.004589798803163333551502
Soils-Ruderal plants            0.004503298759339843804206
Soils-S. lycopersicum           0.009093097562503177355708
                                                       lwr
S. lycopersicum-Ruderal plants -0.017933413016534421879911
Soils-Ruderal plants           -0.007827448366943323687650
Soils-S. lycopersicum          -0.002609815268788376677045
                                                      upr
S. lycopersicum-Ruderal plants 0.008753815410207754776906
Soils-Ruderal plants           0.016834045885623011296062
Soils-S. lycopersicum          0.020796010393794729653738
                                                  p adj
S. lycopersicum-Ruderal plants 0.6927697073340020050836
Soils-Ruderal plants           0.6613140305723815837169
Soils-S. lycopersicum          0.1593884426287284217238

> TukeyHSD(atm)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_TM7"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants  0.012995191820297756490898
Soils-Ruderal plants            0.004991639647612068977911
Soils-S. lycopersicum          -0.008003552172685687512987
                                                       lwr
S. lycopersicum-Ruderal plants  0.002024595092746687627350
Soils-Ruderal plants           -0.005146217450827397787982
Soils-S. lycopersicum          -0.017625228515900877301892
                                                      upr
S. lycopersicum-Ruderal plants 0.023965788547848823619724
Soils-Ruderal plants           0.015129496746051535743804
Soils-S. lycopersicum          0.001618124170529500541194
                                                   p adj
S. lycopersicum-Ruderal plants 0.01589890661870774213327
Soils-Ruderal plants           0.47319069196295726609236
Soils-S. lycopersicum          0.12276672971707792980656

> TukeyHSD(acy)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(fmelt_phylum_stl_rel, Phylum == "p_Cyanobacteria"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants  0.011190279617914310889004
Soils-Ruderal plants            0.006599922080886841749270
Soils-S. lycopersicum          -0.004590357537027469139734
                                                       lwr
S. lycopersicum-Ruderal plants -0.003294230998081194503047
Soils-Ruderal plants           -0.006785119992259433518988
Soils-S. lycopersicum          -0.017293884664861588346119
                                                     upr
S. lycopersicum-Ruderal plants 0.02567479023390981801578
Soils-Ruderal plants           0.01998496415403311701753
Soils-S. lycopersicum          0.00811316959080665006665
                                                  p adj
S. lycopersicum-Ruderal plants 0.1626683616813151544989
Soils-Ruderal plants           0.4721862655968949162855
Soils-S. lycopersicum          0.6670633530799907617848
```


```R
#Con la función HSD.test, se pueden definir los grupos con la opción, group=TRUE. Si se pone en false, 
#se puede mostrar el valor de P, los resultados son los mismos que al usar Tukey.HSD. 
#Estos grupos serán agregadas a la figura de abundancia relativa.
tapr <- HSD.test(apr, "Microbiome_type", group=TRUE)
tave <- HSD.test(ave, "Microbiome_type", group=TRUE)
tafi <- HSD.test(afi, "Microbiome_type", group=TRUE)
taac <- HSD.test(aac, "Microbiome_type", group=TRUE)
taar <- HSD.test(aar, "Microbiome_type", group=TRUE)
tach <- HSD.test(ach, "Microbiome_type", group=TRUE)
taci <- HSD.test(aci, "Microbiome_type", group=TRUE)
tapl <- HSD.test(apl, "Microbiome_type", group=TRUE)
taba <- HSD.test(aba, "Microbiome_type", group=TRUE)
tage <- HSD.test(age, "Microbiome_type", group=TRUE)
tatm <- HSD.test(atm, "Microbiome_type", group=TRUE)
tacy <- HSD.test(acy, "Microbiome_type", group=TRUE)
```


```R
#Grupos de la prueba de Tukey
> tapr
$statistics
    MSerror  Df    Mean     CV
  0.0050076 103 0.34036 20.791

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r     Min     Max     Q25     Q50
Ruderal plants             0.34318 0.080400 27 0.21043 0.56858 0.29076 0.32703
S. lycopersicum            0.43396 0.078734 32 0.26425 0.58457 0.37717 0.44207
Soils                      0.27502 0.058150 47 0.16565 0.43466 0.24547 0.27673
                    Q75
Ruderal plants  0.38640
S. lycopersicum 0.49067
Soils           0.30473

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum            0.43396      a
Ruderal plants             0.34318      b
Soils                      0.27502      c

attr(,"class")
[1] "group"
> tave
$statistics
     MSerror  Df    Mean     CV
  0.00011716 103 0.01898 57.031

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance       std  r       Min      Max      Q25
Ruderal plants            0.021233 0.0086478 27 0.0017985 0.038391 0.016316
S. lycopersicum           0.021891 0.0132283 32 0.0040071 0.064041 0.013198
Soils                     0.015703 0.0101068 47 0.0016358 0.057673 0.010310
                     Q50      Q75
Ruderal plants  0.017884 0.024251
S. lycopersicum 0.020274 0.026949
Soils           0.013343 0.018442

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum           0.021891      a
Ruderal plants            0.021233     ab
Soils                     0.015703      b

attr(,"class")
[1] "group"
> tafi
$statistics
     MSerror  Df      Mean     CV
  7.6636e-05 103 0.0083917 104.32

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance       std  r        Min      Max       Q25
Ruderal plants           0.0064183 0.0079682 27 0.00039966 0.038014 0.0019869
S. lycopersicum          0.0055972 0.0062823 32 0.00023571 0.024912 0.0017680
Soils                    0.0114281 0.0104458 47 0.00049092 0.058859 0.0055968
                      Q50       Q75
Ruderal plants  0.0035141 0.0067939
S. lycopersicum 0.0029750 0.0074273
Soils           0.0075503 0.0147822

$comparison
NULL

$groups
                Relative_abundance groups
Soils                    0.0114281      a
Ruderal plants           0.0064183     ab
S. lycopersicum          0.0055972      b

attr(,"class")
[1] "group"
> taac
$statistics
    MSerror  Df    Mean    CV
  0.0047673 103 0.23549 29.32

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r      Min     Max      Q25
Ruderal plants             0.32515 0.074635 27 0.195345 0.45715 0.265647
S. lycopersicum            0.12777 0.058241 32 0.036285 0.29417 0.099265
Soils                      0.25733 0.072389 47 0.143275 0.45638 0.202107
                    Q50     Q75
Ruderal plants  0.30493 0.38088
S. lycopersicum 0.12294 0.14969
Soils           0.23424 0.31120

$comparison
NULL

$groups
                Relative_abundance groups
Ruderal plants             0.32515      a
Soils                      0.25733      b
S. lycopersicum            0.12777      c

attr(,"class")
[1] "group"
> taar
$statistics
     MSerror  Df      Mean    CV
  3.4264e-05 103 0.0080539 72.68

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance       std  r        Min      Max       Q25
Ruderal plants           0.0060265 0.0037200 27 0.00097940 0.015429 0.0036493
S. lycopersicum          0.0102983 0.0080568 32 0.00169238 0.038311 0.0052836
Soils                    0.0076905 0.0050154 47 0.00044883 0.029176 0.0046757
                      Q50       Q75
Ruderal plants  0.0051326 0.0085738
S. lycopersicum 0.0080990 0.0132758
Soils           0.0065789 0.0102311

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum          0.0102983      a
Soils                    0.0076905     ab
Ruderal plants           0.0060265      b

attr(,"class")
[1] "group"
> tach
$statistics
    MSerror  Df     Mean     CV
  0.0020425 103 0.064337 70.245

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r       Min     Max      Q25
Ruderal plants            0.044263 0.023351 27 0.0077648 0.13229 0.030851
S. lycopersicum           0.054540 0.034242 32 0.0141426 0.16149 0.035897
Soils                     0.082540 0.058949 47 0.0318853 0.31077 0.048937
                     Q50      Q75
Ruderal plants  0.041182 0.050232
S. lycopersicum 0.045650 0.057796
Soils           0.067783 0.081518

$comparison
NULL

$groups
                Relative_abundance groups
Soils                     0.082540      a
S. lycopersicum           0.054540      b
Ruderal plants            0.044263      b

attr(,"class")
[1] "group"
> taci
$statistics
    MSerror  Df    Mean     CV
  0.0019175 103 0.11224 39.013

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r      Min     Max      Q25
Ruderal plants            0.095003 0.042556 27 0.023294 0.21287 0.065714
S. lycopersicum           0.105317 0.047691 32 0.016971 0.19753 0.081414
Soils                     0.126860 0.041680 47 0.062657 0.23219 0.093132
                     Q50     Q75
Ruderal plants  0.089692 0.12744
S. lycopersicum 0.103401 0.13428
Soils           0.124822 0.15638

$comparison
NULL

$groups
                Relative_abundance groups
Soils                     0.126860      a
S. lycopersicum           0.105317     ab
Ruderal plants            0.095003      b

attr(,"class")
[1] "group"
> tapl
$statistics
    MSerror  Df     Mean    CV
  0.0015905 103 0.091011 43.82

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r      Min     Max      Q25
Ruderal plants            0.066181 0.022616 27 0.025932 0.10198 0.049406
S. lycopersicum           0.068696 0.045288 32 0.010652 0.15147 0.028160
Soils                     0.120468 0.043475 47 0.037017 0.23208 0.093157
                     Q50      Q75
Ruderal plants  0.073185 0.081141
S. lycopersicum 0.060919 0.106956
Soils           0.119271 0.152176

$comparison
NULL

$groups
                Relative_abundance groups
Soils                     0.120468      a
S. lycopersicum           0.068696      b
Ruderal plants            0.066181      b

attr(,"class")
[1] "group"
> taba
$statistics
     MSerror  Df     Mean     CV
  0.00094312 103 0.055289 55.545

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r       Min      Max      Q25
Ruderal plants            0.044988 0.018592 27 0.0098378 0.094910 0.031542
S. lycopersicum           0.099235 0.050293 32 0.0210961 0.261388 0.065250
Soils                     0.031286 0.014553 47 0.0081788 0.062737 0.019721
                     Q50      Q75
Ruderal plants  0.045335 0.057275
S. lycopersicum 0.091334 0.119732
Soils           0.028655 0.042870

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum           0.099235      a
Ruderal plants            0.044988      b
Soils                     0.031286      b

attr(,"class")
[1] "group"
> tage
$statistics
     MSerror  Df     Mean     CV
  0.00046106 103 0.028931 74.218

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r       Min      Max      Q25
Ruderal plants            0.028320 0.021351 27 0.0093502 0.084395 0.012683
S. lycopersicum           0.023730 0.015898 32 0.0085745 0.093459 0.013987
Soils                     0.032824 0.024584 47 0.0081739 0.165678 0.019369
                     Q50      Q75
Ruderal plants  0.019560 0.038992
S. lycopersicum 0.020572 0.029463
Soils           0.027879 0.038347

$comparison
NULL

$groups
                Relative_abundance groups
Soils                     0.032824      a
Ruderal plants            0.028320      a
S. lycopersicum           0.023730      a

attr(,"class")
[1] "group"
> tatm
$statistics
     MSerror  Df     Mean     CV
  0.00031165 103 0.015153 116.51

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance       std  r       Min      Max       Q25
Ruderal plants           0.0090163 0.0049373 27 0.0026108 0.023380 0.0050146
S. lycopersicum          0.0220115 0.0164662 32 0.0021834 0.083727 0.0120982
Soils                    0.0140079 0.0223905 47 0.0013676 0.143492 0.0056966
                      Q50      Q75
Ruderal plants  0.0086242 0.010596
S. lycopersicum 0.0161930 0.029875
Soils           0.0079365 0.013033

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum          0.0220115      a
Soils                    0.0140079     ab
Ruderal plants           0.0090163      b

attr(,"class")
[1] "group"
> tacy
$statistics
     MSerror  Df     Mean     CV
  0.00054328 103 0.011094 210.11

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r        Min      Max        Q25
Ruderal plants            0.004789 0.013219 27 0.00025615 0.065562 0.00049005
S. lycopersicum           0.015979 0.023697 32 0.00071189 0.115828 0.00465300
Soils                     0.011389 0.027189 47 0.00000000 0.172595 0.00131386
                       Q50       Q75
Ruderal plants  0.00080638 0.0022538
S. lycopersicum 0.00781763 0.0179700
Soils           0.00316807 0.0092832

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum           0.015979      a
Soils                     0.011389      a
Ruderal plants            0.004789      a

attr(,"class")
[1] "group"

```


```R
#Analizar las diferencias entre Proteobacterias, se utiliza un script similar al anterior
Class_stl <- tax_glom(stl.phy, taxrank="Class")
```


```R
#Obetener abundancia relativa
Class_stl_rel <- transform_sample_counts(Class_stl, function(x) x / sum(x))
```


```R
#Unir la tabla de taxonomía con los otus y con los metadatos
cclass_stl_rel <- data.frame(cbind(tax_table(Class_stl_rel), otu_table(Class_stl_rel)), row.names=NULL)
cclass_stl_rel <- cclass_stl_rel[, colSums(is.na(cclass_stl_rel)) != nrow(cclass_stl_rel)]

tcclass_stl_rel <- data.frame(t(cclass_stl_rel))
tcclass_stl_rel <- tcclass_stl_rel[-1,]
tcclass_stl_rel <- tcclass_stl_rel[-1,]
colnames(tcclass_stl_rel) <- NULL

#Exportar tabla para cargar directamente sin tax_glom
write.table(tcclass_stl_rel, "tcclass_stl_rel.tsv", quote = TRUE, sep = "\t")

#Cargar tabla
tcclass_stl_rel <- read.table("tcclass_stl_rel.tsv", header = TRUE, sep = "\t")


#Agregar los metadatos
microbiome <- data.frame(suelos_data2$microbiome)

mtcclass_stl_rel <- cbind(microbiome, tcclass_stl_rel, row.names = NULL)
mtcclass_stl_rel <- cbind(microbiome, tcclass_stl_rel, row.names = NULL) %>% drop_na()
```


```R
#Graficar abundancia relativa de Alphaproteobacteria
melt_class_stl_rel <- melt(mtcclass_stl_rel)
colnames(melt_class_stl_rel) <- c("Microbiome_type", "Sample", "Class", "Relative_abundance")

class_alpha <- filter(melt_class_stl_rel, Class == "c_Alphaproteobacteria" | 
                               Class == "c_Betaproteobacteria" | Class == "c_Deltaproteobacteria" |
                               Class == "c_Gammaproteobacteria")

#Cambiar los nombres de las variables que se quiere agregar a las gráficas
class_alpha$Microbiome_type <- sub("Plant_local", "Ruderal plants", class_alpha$Microbiome_type)
class_alpha$Microbiome_type <- sub("Plant", "S. lycopersicum", class_alpha$Microbiome_type)
class_alpha$Microbiome_type <- sub("Soil", "Soils", class_alpha$Microbiome_type)

#Ordenar como se quieren graficar
class_alpha$Class_o = factor(class_alpha$Class, levels=c('c_Alphaproteobacteria', 
                                        'c_Betaproteobacteria', 'c_Deltaproteobacteria', 
                                        'c_Gammaproteobacteria'))

#Crear la gráfica
plot_class_proteo <- ggplot(class_alpha, aes(Class, Relative_abundance, fill=Microbiome_type)) + 
                            geom_boxplot() + facet_wrap(~ Class_o, scales="free", nrow = 1) + 
                            scale_fill_manual(values = c("#009E73", "#E69F00", "#CC79A7"))+
                            theme_bw(base_size = 30) + theme(strip.background = element_blank(), 
                                                             strip.text = element_blank()) +
                            theme(axis.text.y = element_text(size = 22), 
                                  axis.text.x = element_text(size = 22)) 

ggsave("plot_class_proteo.pdf", width=80, height=13.3, units="cm")
```


```R
#ANOVA para cada grupo (Clase de Proteobacteria)
aal <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Alphaproteobacteria"))
abe <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Betaproteobacteria"))
ade <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Deltaproteobacteria"))
aga <- aov(Relative_abundance~Microbiome_type, filter(class_alpha, Class == "c_Gammaproteobacteria"))
```


```R
#Ver los resultados de cada ANOVA
> summary(aal)
                 Df Sum Sq Mean Sq F value Pr(>F)    
Microbiome_type   2  0.445  0.2225    58.5 <2e-16 ***
Residuals       103  0.392  0.0038                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(abe)
                 Df Sum Sq Mean Sq F value Pr(>F)    
Microbiome_type   2 0.0382 0.01908    67.8 <2e-16 ***
Residuals       103 0.0290 0.00028                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(ade)
                 Df  Sum Sq  Mean Sq F value  Pr(>F)    
Microbiome_type   2 0.00537 0.002685    19.9 4.9e-08 ***
Residuals       103 0.01390 0.000135                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> summary(aga)
                 Df Sum Sq Mean Sq F value  Pr(>F)    
Microbiome_type   2  0.026 0.01298    12.8 1.1e-05 ***
Residuals       103  0.104 0.00101                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


```R
#Hacer una prueba de Tukey para cada ANOVA
options(digits=22)
TukeyHSD(aal)
TukeyHSD(ade)
TukeyHSD(abe)
TukeyHSD(aga)
options(digits=5)
```


```R
> TukeyHSD(aal)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(class_alpha, Class == "c_Alphaproteobacteria"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants  0.15734986556514246602845
Soils-Ruderal plants            0.02918334298263297599263
Soils-S. lycopersicum          -0.12816652258250949003582
                                                       lwr
S. lycopersicum-Ruderal plants  0.119019957253232472327298
Soils-Ruderal plants           -0.006237076044033125488131
Soils-S. lycopersicum          -0.161783469882753327961922
                                                     upr
S. lycopersicum-Ruderal plants  0.1956797738770524597296
Soils-Ruderal plants            0.0646037620092990705345
Soils-S. lycopersicum          -0.0945495752822656659875
                                                     p adj
S. lycopersicum-Ruderal plants 1.953992519999999886902e-14
Soils-Ruderal plants           1.276056624576952236438e-01
Soils-S. lycopersicum          4.518607709999999737111e-14

> TukeyHSD(ade)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(class_alpha, Class == "c_Deltaproteobacteria"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants -0.018742309906638665584211
Soils-Ruderal plants           -0.013106339546159723186136
Soils-S. lycopersicum           0.005635970360478942398075
                                                        lwr
S. lycopersicum-Ruderal plants -0.0259610879506700364949268
Soils-Ruderal plants           -0.0197771653957649687360920
Soils-S. lycopersicum          -0.0006952026520835605585824
                                                       upr
S. lycopersicum-Ruderal plants -0.011523531862607294673495
Soils-Ruderal plants           -0.006435513696554476768819
Soils-S. lycopersicum           0.011967143373041444487370
                                                     p adj
S. lycopersicum-Ruderal plants 3.993645747435440296594e-08
Soils-Ruderal plants           2.678676650136324610685e-05
Soils-S. lycopersicum          9.137304112856081683702e-02

> TukeyHSD(abe)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(class_alpha, Class == "c_Betaproteobacteria"))

$Microbiome_type
                                                     diff
S. lycopersicum-Ruderal plants -0.01868715005202900569792
Soils-Ruderal plants           -0.04567396915789570521760
Soils-S. lycopersicum          -0.02698681910586670298913
                                                      lwr
S. lycopersicum-Ruderal plants -0.02911128147332506937994
Soils-Ruderal plants           -0.05530684122027461346871
Soils-S. lycopersicum          -0.03612922231197193689844
                                                       upr
S. lycopersicum-Ruderal plants -0.008263018630732942015893
Soils-Ruderal plants           -0.036041097095516796966486
Soils-S. lycopersicum          -0.017844415899761469079809
                                                     p adj
S. lycopersicum-Ruderal plants 1.317536605121016890507e-04
Soils-Ruderal plants           2.553512999999999807751e-15
Soils-S. lycopersicum          7.304666871376999997170e-10

> TukeyHSD(aga)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = Relative_abundance ~ Microbiome_type, data = filter(class_alpha, Class == "c_Gammaproteobacteria"))

$Microbiome_type
                                                      diff
S. lycopersicum-Ruderal plants -0.029132149244737007981465
Soils-Ruderal plants           -0.038541003550726246618119
Soils-S. lycopersicum          -0.009408854305989238636654
                                                      lwr
S. lycopersicum-Ruderal plants -0.04890671047025795664664
Soils-Ruderal plants           -0.05681454699300936272977
Soils-S. lycopersicum          -0.02675197908453785106131
                                                       upr
S. lycopersicum-Ruderal plants -0.009357588019216059316285
Soils-Ruderal plants           -0.020267460108443133975920
Soils-S. lycopersicum           0.007934270472559373788002
                                                     p adj
S. lycopersicum-Ruderal plants 1.954041696476394029958e-03
Soils-Ruderal plants           6.583786555269810492064e-06
Soils-S. lycopersicum          4.039503678317045709534e-01


```


```R
#Con la función HSD.test, se pueden definir los grupos con la opción, group=TRUE. Si se pone en false, 
#se puede mostrar el valor de P, los resultados son los mismos que al usar Tukey.HSD. 
#Estos grupos serán agregadas a la figura de abundancia relativa.
taal <- HSD.test(aal, "Microbiome_type", group=TRUE)
tade <- HSD.test(ade, "Microbiome_type", group=TRUE)
tabe <- HSD.test(abe, "Microbiome_type", group=TRUE)
taga <- HSD.test(aga, "Microbiome_type", group=TRUE)
```


```R
#Mostrar los grupos
> taal
$statistics
    MSerror  Df    Mean     CV
  0.0038044 103 0.24186 25.502

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r     Min     Max     Q25     Q50
Ruderal plants             0.18142 0.048902 27 0.10005 0.30563 0.15121 0.18908
S. lycopersicum            0.33877 0.082411 32 0.18707 0.53574 0.27718 0.33138
Soils                      0.21060 0.050892 47 0.10277 0.32870 0.16819 0.21249
                    Q75
Ruderal plants  0.20379
S. lycopersicum 0.39907
Soils           0.23705

$comparison
NULL

$groups
                Relative_abundance groups
S. lycopersicum            0.33877      a
Soils                      0.21060      b
Ruderal plants             0.18142      b

attr(,"class")
[1] "group"
> tabe 
$statistics
     MSerror  Df     Mean     CV
  0.00028138 103 0.036314 46.192

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance       std  r       Min      Max      Q25
Ruderal plants            0.062207 0.0221762 27 0.0280902 0.114009 0.045613
S. lycopersicum           0.043520 0.0211923 32 0.0130819 0.112924 0.031849
Soils                     0.016533 0.0070297 47 0.0010905 0.031319 0.011261
                     Q50      Q75
Ruderal plants  0.059806 0.071693
S. lycopersicum 0.037703 0.053771
Soils           0.017397 0.020921

$comparison
NULL

$groups
                Relative_abundance groups
Ruderal plants            0.062207      a
S. lycopersicum           0.043520      b
Soils                     0.016533      c

attr(,"class")
[1] "group"
> tade 
$statistics
     MSerror  Df     Mean    CV
  0.00013494 103 0.031687 36.66

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance       std  r       Min      Max      Q25
Ruderal plants            0.043156 0.0141824 27 0.0136169 0.066394 0.032815
S. lycopersicum           0.024414 0.0095357 32 0.0041725 0.043586 0.017437
Soils                     0.030050 0.0112775 47 0.0065559 0.054993 0.022763
                     Q50      Q75
Ruderal plants  0.040907 0.054220
S. lycopersicum 0.025403 0.030738
Soils           0.028274 0.038360

$comparison
NULL

$groups
                Relative_abundance groups
Ruderal plants            0.043156      a
Soils                     0.030050      b
S. lycopersicum           0.024414      b

attr(,"class")
[1] "group"
> taga
$statistics
    MSerror  Df     Mean     CV
  0.0010126 103 0.030314 104.97

$parameters
   test          name.t ntr StudentizedRange alpha
  Tukey Microbiome_type   3           3.3631  0.05

$means
                Relative_abundance      std  r       Min      Max      Q25
Ruderal plants            0.056197 0.060675 27 0.0131316 0.321395 0.026219
S. lycopersicum           0.027065 0.013398 32 0.0072323 0.073769 0.018578
Soils                     0.017656 0.008094 47 0.0040907 0.033182 0.010844
                     Q50      Q75
Ruderal plants  0.039891 0.052640
S. lycopersicum 0.024228 0.034162
Soils           0.018607 0.023415

$comparison
NULL

$groups
                Relative_abundance groups
Ruderal plants            0.056197      a
S. lycopersicum           0.027065      b
Soils                     0.017656      b

attr(,"class")
[1] "group"
```
