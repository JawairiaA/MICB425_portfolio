library("tidyverse")
library(tidyr)
library("magrittr")
library("cowplot")

RNA_10m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/10m/10m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_10m["Domain"]<-RNA_10m$Taxonomy
RNA_10m[8] <- lapply(RNA_10m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_10m_sep = separate(data=RNA_10m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_10m_sep["Depth"]<- "10_m"
R1<-select(RNA_10m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R1["Perc"]<-R1[2]*(1/sum(R1$Abundance))


RNA_100m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/100m/100m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_100m["Domain"]<-RNA_100m$Taxonomy
RNA_100m[8] <- lapply(RNA_100m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_100m_sep = separate(data=RNA_100m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_100m_sep["Depth"]<- "100_m"
R2<-select(RNA_100m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R2["Perc"]<-R2[2]*(1/sum(R2$Abundance))

RNA_120m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/120m/120m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_120m["Domain"]<-RNA_120m$Taxonomy
RNA_120m[8] <- lapply(RNA_120m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_120m_sep = separate(data=RNA_120m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_120m_sep["Depth"]<- "120_m"
R3<-select(RNA_120m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R3["Perc"]<-R3[2]*(1/sum(R3$Abundance))


RNA_135m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/135m/135m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_135m["Domain"]<-RNA_135m$Taxonomy
RNA_135m[8] <- lapply(RNA_135m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_135m_sep = separate(data=RNA_135m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_135m_sep["Depth"]<- "135_m"
R4<-select(RNA_135m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R4["Perc"]<-R4[2]*(1/sum(R4$Abundance))

RNA_150m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/150m/150m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_150m["Domain"]<-RNA_150m$Taxonomy
RNA_150m[8] <- lapply(RNA_150m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_150m_sep = separate(data=RNA_150m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_150m_sep["Depth"]<- "150_m"
R5<-select(RNA_150m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R5["Perc"]<-R5[2]*(1/sum(R5$Abundance))

RNA_165m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/165m/165m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_165m["Domain"]<-RNA_165m$Taxonomy
RNA_165m[8] <- lapply(RNA_165m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_165m_sep = separate(data=RNA_165m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_165m_sep["Depth"]<- "165_m"
R6<-select(RNA_165m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R6["Perc"]<-R6[2]*(1/sum(R6$Abundance))

RNA_200m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/200m/200m_RNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
RNA_200m["Domain"]<-RNA_200m$Taxonomy
RNA_200m[8] <- lapply(RNA_200m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
RNA_200m_sep = separate(data=RNA_200m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
RNA_200m_sep["Depth"]<- "200_m"
R7<-select(RNA_200m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
R7["Perc"]<-R7[2]*(1/sum(R7$Abundance))


R_all=rbind(R1, R2, R3, R4, R5, R6, R7)


#Order
p<-ggplot(data=R_all, aes(x=Depth, y= Perc, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")
p

P_unstandardized<-
  ggplot(data=R_all, aes(x=Depth, y= Abundance, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")
P_unstandardized

#Phylum
R_all%>%
  ggplot(aes(x=Depth, y= Perc, fill=Phylum))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")

#Family
R_all%>%
  ggplot(aes(x=Depth, y= Perc, fill=Family))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")

#Genus
R_all%>%
  ggplot(aes(x=Depth, y= Perc, fill=Genus))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")

____________

library("tidyverse")
library(tidyr)
library("magrittr")
library("cowplot")

DNA_10m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/10m/10m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_10m["Domain"]<-DNA_10m$Taxonomy
DNA_10m[8] <- lapply(DNA_10m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_10m_sep = separate(data=DNA_10m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_10m_sep["Depth"]<- "10_m"
D1<-select(DNA_10m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D1["Perc"]<-D1[2]*(1/sum(D1$Abundance))

DNA_100m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/100m/100m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_100m["Domain"]<-DNA_100m$Taxonomy
DNA_100m[8] <- lapply(DNA_100m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_100m_sep = separate(data=DNA_100m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_100m_sep["Depth"]<- "100_m"
D2<-select(DNA_100m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D2["Perc"]<-D2[2]*(1/sum(D2$Abundance))

DNA_120m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/120m/120m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_120m["Domain"]<-DNA_120m$Taxonomy
DNA_120m[8] <- lapply(DNA_120m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_120m_sep = separate(data=DNA_120m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_120m_sep["Depth"]<- "120_m"
D3<-select(DNA_120m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D3["Perc"]<-D3[2]*(1/sum(D3$Abundance))


DNA_135m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/135m/135m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_135m["Domain"]<-DNA_135m$Taxonomy
DNA_135m[8] <- lapply(DNA_135m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_135m_sep = separate(data=DNA_135m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_135m_sep["Depth"]<- "135_m"
D4<-select(DNA_135m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D4["Perc"]<-D4[2]*(1/sum(D4$Abundance))

DNA_150m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/150m/150m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_150m["Domain"]<-DNA_150m$Taxonomy
DNA_150m[8] <- lapply(DNA_150m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_150m_sep = separate(data=DNA_150m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_150m_sep["Depth"]<- "150_m"
D5<-select(DNA_150m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D5["Perc"]<-D5[2]*(1/sum(D5$Abundance))

DNA_165m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/165m/165m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_165m["Domain"]<-DNA_165m$Taxonomy
DNA_165m[8] <- lapply(DNA_165m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_165m_sep = separate(data=DNA_165m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_165m_sep["Depth"]<- "165_m"
D6<-select(DNA_165m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D6["Perc"]<-D6[2]*(1/sum(D6$Abundance))

DNA_200m<-read.delim("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/200m/200m_DNA/final_outputs/marker_contig_map.tsv", na.strings=c(""))
DNA_200m["Domain"]<-DNA_200m$Taxonomy
DNA_200m[8] <- lapply(DNA_200m[8], gsub, pattern = "Cystobacterineae;", replacement = "", fixed = TRUE)
DNA_200m_sep = separate(data=DNA_200m, col=Domain, into=c("Domain", "Phylum",	"Class","Order","Family",	"Genus",	"Species", "Subsp."), sep="\\;", na.strings=c(""))
DNA_200m_sep["Depth"]<- "200_m"
D7<-select(DNA_200m_sep, Depth, Abundance, Domain, Phylum, Class, Order, Family, Genus, Species, Subsp.)
D7["Perc"]<-D7[2]*(1/sum(D7$Abundance))


D_all=rbind(D1, D2, D3, D4, D5, D6, D7)

#Order
P2<-ggplot(data=D_all, aes(x=Depth, y= Perc, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")
P2


P2_unstandardized<-
  ggplot(data=D_all, aes(x=Depth, y= Abundance, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")
P2_unstandardized

#Phylum
D_all%>%
  ggplot(aes(x=Depth, y= Perc, fill=Phylum))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")

#Family
D_all%>%
  ggplot(aes(x=Depth, y= Perc, fill=Family))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")

#Genus
D_all%>%
  ggplot(aes(x=Depth, y= Perc, fill=Genus))+
  geom_bar(stat="identity")+
  labs(title="WHY", y="Relative Abundance")
#_____


RNA_unstandardized<-
  ggplot(data=R_all, aes(x=Depth, y= Abundance, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA expression: RPKM subsetted by Bacterial Orders across samples", y="Sum RPKM")+ 
  theme(legend.position="bottom",legend.text=element_text(size=11),legend.title=element_text(size=14))
RNA_unstandardized

RNA<-ggplot(data=R_all, aes(x=Depth, y= Perc, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA expression: Relative Abundance of Bacterial Orders across Depths", y="Relative Abundance")+
  theme(legend.position="bottom",legend.text=element_text(size=11),legend.title=element_text(size=14))
RNA

DNA_unstandardized<-
  ggplot(data=D_all, aes(x=Depth, y= Abundance, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA in genome: RPKM subsetted by Bacterial Orders across samples", y="Sum RPKM")+
  theme(legend.position="bottom",legend.text=element_text(size=11),legend.title=element_text(size=14))
DNA_unstandardized

DNA<-ggplot(data=D_all, aes(x=Depth, y= Perc, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA in genome: Relative Abundance of Bacterial Orders across Depths", y="Relative Abundance")+
  theme(legend.position="bottom",legend.text=element_text(size=11),legend.title=element_text(size=14))
DNA

#____
#Plots by replacing NA for unclassified: for Order

#RNA unstandardized
RNA_unstandardized<-R_all %>%
  mutate(Order = ifelse(is.na(Order), "unclassified", Order)) %>% 
  ggplot(aes(x=Depth, y= Abundance, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA expression: RPKM subsetted by Bacterial Orders across samples", y="Sum RPKM")
RNA_unstandardized

#RNA standardized
RNA<-R_all %>%
  mutate(Order = ifelse(is.na(Order), "unclassified", Order)) %>% 
  ggplot(aes(x=Depth, y= Perc, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA expression: Relative Abundance of Bacterial Orders across Depths", y="Relative Abundance")
RNA

#DNA unstandardized

DNA_unstandardized<- D_all %>%
  mutate(Order = ifelse(is.na(Order), "unclassified", Order)) %>% 
  ggplot(aes(x=Depth, y= Abundance, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA in genome: RPKM subsetted by Bacterial Orders across samples", y="Sum RPKM")
DNA_unstandardized

#DNA standardized
DNA<-D_all %>%
  mutate(Order = ifelse(is.na(Order), "unclassified", Order)) %>% 
  ggplot(aes(x=Depth, y= Perc, fill=Order))+
  geom_bar(stat="identity")+
  labs(title="napA in genome: Relative Abundance of Bacterial Orders across Depths", y="Relative Abundance")
DNA

#_____________


RNA_10<-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/10m/10m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.10 = Confident_Taxonomy, Abund.RNA.10 = Abundance, Query)
RNA_100<-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/100m/100m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.100 = Confident_Taxonomy, Abund.RNA.100 = Abundance, Query)
RNA_120 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/120m/120m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.120 = Confident_Taxonomy, Abund.RNA.120 = Abundance, Query)
RNA_135 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/135m/135m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.135= Confident_Taxonomy, Abund.RNA.135= Abundance, Query)
RNA_150 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/150m/150m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.150= Confident_Taxonomy, Abund.RNA.150= Abundance, Query)
RNA_165 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/165m/165m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.165= Confident_Taxonomy, Abund.RNA.165= Abundance, Query)
RNA_200 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/200m/200m_RNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.RNA.200= Confident_Taxonomy, Abund.RNA.200= Abundance, Query)





DNA_10<-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/10m/10m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.10 = Confident_Taxonomy, Abund.DNA.10 = Abundance, Query)
DNA_100<-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/100m/100m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.100 = Confident_Taxonomy, Abund.DNA.100 = Abundance, Query)
DNA_120 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/120m/120m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.120 = Confident_Taxonomy, Abund.DNA.120 = Abundance, Query)
DNA_135 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/135m/135m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.135= Confident_Taxonomy, Abund.DNA.135= Abundance, Query)
DNA_150 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/150m/150m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.150= Confident_Taxonomy, Abund.DNA.150= Abundance, Query)
DNA_165 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/165m/165m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.165= Confident_Taxonomy, Abund.DNA.165= Abundance, Query)
DNA_200 <-read_tsv("~/MICB425_portfolio/MICB 425 Project 2/gcloud raw outputs/200m/200m_DNA/final_outputs/marker_contig_map.tsv")%>%
  select(Tax.DNA.200= Confident_Taxonomy, Abund.DNA.200= Abundance, Query)


napA.all= DNA_10 %>%
  full_join(RNA_10, by="Query") %>%
  full_join(DNA_100, by="Query") %>%
  full_join(RNA_100, by="Query") %>%
  full_join(DNA_120, by="Query") %>%
  full_join(RNA_120, by="Query") %>%
  full_join(DNA_135, by="Query") %>%
  full_join(RNA_135, by="Query") %>%
  full_join(DNA_150, by="Query") %>%
  full_join(RNA_150, by="Query") %>%
  full_join(DNA_165, by="Query") %>%
  full_join(RNA_165, by="Query") %>%
  full_join(DNA_200, by="Query") %>%
  full_join(RNA_200, by="Query")%>%

  mutate(Taxonomy = ifelse(!is.na(Tax.RNA.10), Tax.RNA.10,
                         ifelse(!is.na(Tax.DNA.10), Tax.DNA.10,
                                ifelse(!is.na(Tax.RNA.100), Tax.RNA.100,
                                       ifelse(!is.na(Tax.DNA.100), Tax.DNA.100,
                                              ifelse(!is.na(Tax.RNA.120), Tax.RNA.120,
                                                ifelse(!is.na(Tax.DNA.120), Tax.DNA.120,
                                                ifelse(!is.na(Tax.RNA.135), Tax.RNA.135,
                                                     ifelse(!is.na(Tax.DNA.135), Tax.DNA.135,
                                                            ifelse(!is.na(Tax.RNA.150), Tax.RNA.150,
                                                                   ifelse(!is.na(Tax.DNA.150), Tax.DNA.150,
                                                                          ifelse(!is.na(Tax.RNA.165), Tax.RNA.165,
                                                                                 ifelse(!is.na(Tax.DNA.165), Tax.DNA.165,
                                                                                        ifelse(!is.na(Tax.RNA.200), Tax.RNA.200,
                                                                                               ifelse(!is.na(Tax.DNA.200), Tax.DNA.200,
                                              "unclassified"))))))))))))))) %>%
                                                select(-starts_with("Tax.")) %>% 
                                                gather("Key", "Abundance", starts_with("Abund")) %>% 
                                                separate(Key, c("Key","Type","Depth_m"), by = ".") %>% 
                                                select(Depth_m, Type, Abundance, Taxonomy, Query) %>% 
                                                mutate(Depth_m = as.numeric(Depth_m)) %>% 
                                                separate(Taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep="; ")

napA.all %>%
  filter(Type == "DNA") %>%
  mutate (Genus= ifelse(is.na(Genus), "unclassified", Genus))%>%
  ggplot(aes(x = "napA", y = Depth_m)) +
  # Use the size aesthetic to show abundance
  geom_point(aes(size = Abundance)) +
  # Reverse the why axis so depth increases going down
  scale_y_reverse(lim=c(200,10)) +
  labs(title = "Abundance of the mcrA gene (DNA) at different depths",
       x = "") +
  theme_classic()


napA.all %>% 
  # Change NAs to "unclassified" at the level you want to plot
  mutate(Genus = ifelse(is.na(Genus), "unclassified", Genus)) %>% 
  
  # Show both RNA and DNA using an x variable  
  ggplot(aes(x = Type, y = Depth_m)) +
  geom_point(aes(size = Abundance)) +
  scale_y_reverse(lim=c(200,10)) +
  labs(title = "Abundance of the napA gene (DNA vs. RNA) at different depths",
       x = "") +
  theme_classic()
                           
               
napA.all %>% 
  # Change NAs to "unclassified" at the level you want to plot
  mutate(Genus = ifelse(is.na(Family), "unclassified", Class)) %>% 
  
  ggplot(aes(x = Family, y = Depth_m)) +
  # Use an ifelse statement to make 0 values into NAs so that they don't show up on the plot
  # Use position_dodge to keep points from overlapping
  geom_point(aes(x= Family, y= Depth_m, size = ifelse(Abundance == 0, NA, Abundance), color = Type), position = position_dodge(0.5)) +
  scale_y_reverse(lim=c(200,10)) +
  labs(title = "Abundance of Class with napA (DNA vs. RNA) at different depths") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  # Rename legend
  scale_size_continuous(name = "Abundance")

library("phyloseq")

load("~/MICB425_portfolio/mothur_phyloseq.RData")

metadata = data.frame(mothur@sam_data)

plot1 = napA.all %>% 
  # Change NAs to "unclassified" at the level you want to plot
  mutate(Phylum = ifelse(is.na(Phylum), "unclassified", Phylum)) %>% 
  
  ggplot(aes(x = Phylum, y = Depth_m)) +
  geom_point(aes(size = ifelse(Abundance == 0, NA, Abundance), color = Type), position = position_dodge(0.5)) +
  scale_y_reverse(lim=c(200,10)) +
  labs(y = "") +
  theme_classic() +
  scale_size_continuous(name = "Abundance")+
  theme(axis.text.x = element_text(angle = 70, hjust = 1))

plot2=metadata %>% 
  select (Depth_m, NO2_uM, NO3_uM, NH4_uM, N2O_nM) %>%
  gather("Nutrients","uM", NO2_uM:N2O_nM)%>%
  arrange(Depth_m) %>% 
  ggplot() + geom_point(aes(y=Depth_m, x=uM)) +
  #geom_line(aes(y=Depth_m, x=uM)) + 
  geom_path(aes(y=Depth_m, x=uM,group = 1)) +
  scale_y_reverse(lim=c(200,10)) +
  facet_wrap(~Nutrients, scales="free_x", nrow=1) +
  theme(legend.position="none")+
  #unfortunately they didn't have N2
  labs(y = "Depth (m)",
       x = "Nitrogen Species")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

plot_grid(plot2, plot1, labels=c("A", "B"), rel_widths=c(1/3, 2/3), align = 'h', axis="b")


____________



Sum_reads_DNA=data.frame(
  Depth_m=c(10, 100, 120, 135, 150, 165, 200),
  Sum_RPKM=c(sum(D1$Abundance), sum(D2$Abundance), sum(D3$Abundance), sum(D4$Abundance),sum(D5$Abundance),sum(D6$Abundance), sum(D7$Abundance)))

Nitrogen=metadata %>% 
  select (NO2_uM, NO3_uM, NH4_uM, N2O_nM)

N_RPKM_DNA=cbind(Sum_reads_DNA, Nitrogen)

plot3=N_RPKM_DNA %>%  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) +
  geom_line(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) + 
  geom_line(aes(x=Depth_m, y=NO2_uM*5000, colour="NO2_uM")) + 
  geom_point(aes(x=Depth_m, y=NO2_uM*5000, colour="NO2_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(1/5000), name = "Nitrogen Species (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="napA gene abundance and NO2 (uM) across depth", y = "Sum RPKM" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.2, 0.9))

plot4=N_RPKM_DNA %>%  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) +
  geom_line(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) + 
  geom_line(aes(x=Depth_m, y=NO3_uM*50, colour="NO3_uM")) + 
  geom_point(aes(x=Depth_m, y=NO3_uM*50, colour="NO3_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(1/50), name = "Nitrogen Species (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="napA gene abundance and NO3 (uM) across depth", y = "Sum RPKM" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.2, 0.9))



Sum_reads_RNA=data.frame(
  Depth_m=c(10, 100, 120, 135, 150, 165, 200),
  Sum_RPKM=c(sum(R1$Abundance), sum(R2$Abundance), sum(R3$Abundance), sum(R4$Abundance),sum(R5$Abundance),sum(R6$Abundance), sum(R7$Abundance)))

Nitrogen=metadata %>% 
  select (NO2_uM, NO3_uM, NH4_uM, N2O_nM)

N_RPKM_RNA=cbind(Sum_reads_RNA, Nitrogen)

plot5=N_RPKM_RNA %>%  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) +
  geom_line(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) + 
  geom_line(aes(x=Depth_m, y=NO2_uM*10000, colour="NO2_uM")) + 
  geom_point(aes(x=Depth_m, y=NO2_uM*10000, colour="NO2_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(1/10000), name = "Nitrogen Species (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="napA gene expression and NO2 (uM) across depth", y = "Sum RPKM" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.2, 0.9))

plot6=N_RPKM_RNA %>%  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) +
  geom_line(aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM")) + 
  geom_line(aes(x=Depth_m, y=NO3_uM*100, colour="NO3_uM")) + 
  geom_point(aes(x=Depth_m, y=NO3_uM*100, colour="NO3_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(1/100), name = "Nitrogen Species (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="napA gene expression and NO3 (uM) across depth", y = "Sum RPKM" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.2, 0.9))

Nit=metadata %>% 
  select (Depth_m, NO2_uM, NO3_uM, NH4_uM, N2O_nM)

# both RNA and DNA with NO2 uM
plot7=ggplot() +
  geom_point(data= N_RPKM_RNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM RNA")) +
  geom_line(data= N_RPKM_RNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM RNA")) + 
  geom_point(data=N_RPKM_DNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM DNA")) +
  geom_line(data=N_RPKM_DNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM DNA")) + 
  geom_point(data=Nit, aes(x=Depth_m, y=NO2_uM*10000, colour="NO2_uM"))+
  geom_line(data=Nit, aes(x=Depth_m, y=NO2_uM*10000, colour="NO2_uM")) + 
  scale_y_continuous(sec.axis = sec_axis(~.*(1/10000), name = "Nitrogen Species (uM)")) +
  #scale_colour_manual(values = c("blue", "red", green))+
  labs(title="napA gene abundance and expression and NO2 (uM) across depth", y = "Sum RPKM" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.2, 0.9), plot.title=element_text(size=14))

# both RNA and DNA with NO3 uM
plot8=ggplot() +
  geom_point(data= N_RPKM_RNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM RNA")) +
  geom_line(data= N_RPKM_RNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM RNA")) + 
  geom_point(data=N_RPKM_DNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM DNA")) +
  geom_line(data=N_RPKM_DNA, aes(x=Depth_m, y=Sum_RPKM, colour= "SUM RPKM DNA")) + 
  geom_point(data=Nit, aes(x=Depth_m, y=NO3_uM*100, colour="NO3_uM"))+
  geom_line(data=Nit, aes(x=Depth_m, y=NO3_uM*100, colour="NO3_uM")) + 
  scale_y_continuous(sec.axis = sec_axis(~.*(1/100), name = "Nitrogen Species (uM)")) +
  #scale_colour_manual(values = c("blue", "red", green))+
  labs(title="napA gene abundance expression and NO3 (uM) across depth", y = "Sum RPKM" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.2, 0.9), plot.title=element_text(size=14))




plot_grid(plot7, plot8, labels=c("A", "B"), rel_widths=c(1/2, 1/2))



