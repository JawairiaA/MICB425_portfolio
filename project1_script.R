library(tidyverse)
library(ggplot2)

#Generating table of specific OTUs and their classifications
mothur_taxa=data.frame(mothur@tax_table)
mothur_otu=data.frame(mothur@otu_table)
mothur_taxa_t=t(mothur_taxa)
total=rbind(mothur_otu, mothur_taxa_t)
total_t=t(total)
total_t2=data.frame(total_t)
total_taxon=total_t2%>%
  filter(Genus== "Candidatus_Scalindua")


#Generating abundance vs depth plot for taxon
gp = subset_taxa(mothur, Genus== "Candidatus_Scalindua")
plot_bar(gp, fill="Species")

#Same for Qiime2 below:

qiime2_taxa=data.frame(qiime2@tax_table)
qiime2_otu=data.frame(qiime2@otu_table)
qiime2_otu_t=t(qiime2_otu)
qiime2_taxa_t=t(qiime2_taxa)
qiime2_total=rbind(qiime2_otu_t, qiime2_taxa_t)
qiime2_total_t=t(qiime2_total)
qiime2_total_t2=data.frame(qiime2_total_t)
qiime2_total_tax=qiime2_total_t2%>%filter(Genus== "D_5__Candidatus Scalindua")

gp2 = subset_taxa(qiime2, Genus== "D_5__Candidatus Scalindua")
plot_bar(gp2, fill="Genus")

-------------

library(tidyverse)
library(phyloseq)

load("mothur_phyloseq.RData")

set.seed(9376)
m.norm = rarefy_even_depth(mothur, sample.size=100000)
m.perc = transform_sample_counts(m.norm, function(x) 100 * x/sum(x))

m.alpha = estimate_richness(m.norm, measures = c("Chao1", "Shannon"))

m.meta.alpha = full_join(rownames_to_column(m.alpha), rownames_to_column(data.frame(m.perc@sam_data)), by = "rowname")

#generating Alpha-diversity + oxygen across depth plot
m.meta.alpha %>%  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Shannon, colour= " Shannon Diversity")) +
  geom_smooth(method='auto', aes(x=as.numeric(Depth_m), y=Shannon)) +
  labs(title="Alpha-diversity across depth", y="Shannon's diversity index", x="Depth (m)") +
  geom_line(aes(x=Depth_m, y=O2_uM/15, colour="O2_uM")) + 
  geom_point(aes(x=Depth_m, y=O2_uM/15, colour="O2_uM"))+
  scale_y_continuous(sec.axis = sec_axis(~.*(15), name = "O2 (uM)")) +
  scale_colour_manual(values = c("blue", "red"))+
  labs(title="Alpha-diversity and Oxygen across depth", y = "Shannon's diversity index" , x = "Depth (m)" , colour = "Parameter") +
  theme(legend.position = c(0.8, 0.9))

#gnerating alpha diversity across oxygen plot

m.meta.alpha %>% 
  ggplot() +
  geom_point(aes(x=O2_uM, y=Shannon)) +
  labs(title="Alpha-diversity across oxygen", y="Shannon's diversity index", x="Oxygen (uM)")

#generating oxic/anoxic vs alpha divrsity plot
m.meta.alpha %>% 
  mutate(O2_group = ifelse(O2_uM == 0, "anoxic", "oxic")) %>% 
  ggplot() +
  geom_boxplot(aes(x=O2_group, y=Shannon)) +
  labs(title="Alpha-diversity by oxic/anoxic", y="Shannon's diversity index", x="Oxygen")

write.table(total_taxon, file = "Taxon.csv", sep = ",", col.names = NA)

#Phyla across depth
m.perc %>% 
  plot_bar(fill="Class") + 
  geom_bar(aes(fill=Class), stat="identity") +
  labs(title="Class across samples")

m.perc %>% 
  
  plot_bar() + 
  geom_bar(aes(fill=Phylum), stat="identity") +
  facet_wrap(~Class, scales="free_y", nrow=10)+
  labs(title="Clases across samples; by Phylum")


m.norm %>% 
  subset_taxa(Genus=="Candidatus_Scalindua") %>% 
  tax_glom(taxrank = 'Genus') %>%
  psmelt() %>%
  
  lm(Abundance ~ Depth_m, .) %>% 
  summary()

#library(magrittr)
m.perc %>% 
  subset_taxa(Genus=="Candidatus_Scalindua") %>% 
  psmelt() %>% 
  group_by(Sample) %>% 
  summarize(Abundance_sum=sum(Abundance), Depth_m=mean(Depth_m)) %>% 
  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Abundance_sum)) +
  geom_smooth(method='lm', aes(x=as.numeric(Depth_m), y=Abundance_sum)) +
  labs(title="Abundance Candidatus Scalindua across depth")

#Linear model: Abundance Against Oxygen
m.perc %>% 
  subset_taxa(Genus=="Candidatus_Scalindua") %>% 
  psmelt() %>% 
  group_by(Sample) %>% 
  summarize(Abundance_sum=sum(Abundance), O2_uM=mean(O2_uM)) %>% 
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=Abundance_sum)) +
  geom_smooth(method='lm', aes(x=as.numeric(O2_uM), y=Abundance_sum)) +
  labs(title="Abundance Candidatus Scalindua across Oxygen conc.")

#Richness across taxon:
m.norm %>% 
  subset_taxa(Genus== "Candidatus_Scalindua")

m.norm %>% 
  subset_taxa(Genus== "Candidatus_Scalindua") %>%
  estimate_richness(measures = c("Observed"))

m.norm %>% 
  psmelt() %>% 
  filter(OTU=="Otu0242") %>% 
  lm(Abundance ~ Depth_m, .) %>% 
  summary()



m.perc %>% 
  subset_taxa(Genus == "Candidatus_Scalindua") %>% 
  psmelt() %>% 
  filter(OTU=="Otu1922") %>% 
  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Abundance)) +
  geom_smooth(method='lm', aes(x=Depth_m, y=Abundance)) +
  facet_wrap(~OTU, scales="free_y") +
  labs(title="Abundance of OTUs within Genus Candidatus Scalindua domain across depth")

m.perc %>% 
  subset_taxa(Genus == "Candidatus_Scalindua") %>% 
  psmelt() %>% 
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=Abundance)) +
  geom_smooth(method='lm', aes(x=O2_uM, y=Abundance)) +
  facet_wrap(~OTU, scales="free_y") +
  labs(title="Abundance of OTUs within Genus Candidatus Scalindua domain across Oxygen conc.")


m.norm %>% 
  psmelt() %>% 
  filter(OTU=="Otu0031") %>% 
  
  lm(Abundance ~ Depth_m, .) %>% 
  summary()
#OTU0031 p-value: 0.2156

#All otus: correcting for multiple comparisons
p.adjust(c(0.2156,0.2074, 0.2663, 0.2077, 0.5117,
           0.2077, 0.1701, 0.5322, 0.2077, 0.1191,
           0.4276, 0.885, 0.2077, 0.0164, 0.5322,
           0.7658, 0.2077, 0.2077, 0.2077, 0.2077,
           0.5322), method="fdr")


#Abundance of OTUs within Candidatus Scalindua domain across depth
m.perc %>%
  subset_taxa(Genus=="Candidatus_Scalindua") %>%
  psmelt() %>%
  
  ggplot() +
  geom_point(aes(x=Sample, y=OTU, size=Abundance, color=OTU)) +
  scale_size_continuous(range = c(0,5)) +
  labs(title="Abundance of OTUs within Candidatus Scalindua domain across depth")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Abundance of OTUs within Candidatus Scalindua domain across oxygen
m.perc %>% 
  subset_taxa(Genus=="Candidatus_Scalindua") %>%
  psmelt() %>% 
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=OTU, size=Abundance, color=OTU)) + 
  scale_size_continuous(range = c(0,5)) +
  labs(title="Abundance of OTUs within Candidatus Scalindua domain with O2 conc.")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#########

#Do the abundances of OTUs/ASVs within your taxon of interest change significantly with Ocygen conc.

m.norm %>% 
  psmelt() %>% 
  filter(OTU=="Otu0031") %>% 
  
  lm(Abundance ~ O2_uM, .) %>% 
  summary()

p.adjust(c(0.2637, 0.5091, 0.3793, 0.5095, 0.5558,
           0.5095, 0.5095, 0.5431, 0.5095, 0.5095,
           0.393, 0.7762, 0.5095, 0.5095, 0.5095,
           0.7507, 0.5095, 0.5095, 0.5095, 0.5095,
           0.5095), method="fdr")

#Abundance of bacterial Phyla across oxygen
m.perc %>%
  subset_taxa(Domain="Bacteria") %>%
  psmelt() %>%
  
  ggplot() +
  geom_point(aes(x=O2_uM, y=Phylum, size=Abundance, color=Phylum)) +
  scale_size_continuous(range = c(0,5)) +
  labs(title="Abundance of Bacterial Phyla with O2 conc.")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#Abundance of bacterial Phyla across depth
m.perc %>%
  subset_taxa(Domain="Bacteria") %>%
  psmelt() %>%
  
  ggplot() +
  geom_point(aes(x=Depth_m, y=Phylum, size=Abundance, color=Phylum)) +
  scale_size_continuous(range = c(0,5)) +
  labs(title="Abundance of Bacterial Phyla with Depth")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))






