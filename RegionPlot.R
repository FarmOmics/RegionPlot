library(data.table)
library(tidyverse)


data_file = 
region=c(28, 3430950, 3630950)
chr_name="chr"
snp="variant_id"
bp="pos"
p="pval_gi"
plink=
bfile=

LD.colours <- data.frame(LD = as.character(seq(from = 0, to = 1, by = 0.1)), Colour = c("#000080",rep(c("#000080", "#87CEFA", "#00FF00", "#FFA500", "#FF0000"), each = 2)), stringsAsFactors = FALSE)

### define plot region
chr=as.integer(region[1])
start=as.integer(region[2])
end=as.integer(region[3])


#### load data
mydata=fread(data_file, header=T) %>% as.data.frame
mydata %>% select(all_of(chr_name), all_of(snp), all_of(bp), all_of(p)) %>% as.data.frame -> mydata
#df %>% select(all_of(chr_name), all_of(snp), all_of(bp), all_of(p)) %>% as.data.frame -> mydata
names(mydata)=c("CHR", "SNP", "BP", "P")

#### calculate ld
ld_file=tempfile(pattern = "ld_file", "Temp")
cmd=paste0("/share/apps/plink-1.90/plink --bfile ", bfile, " --chr-set 39 --keep-allele-order --chr ", chr, " --r2 --out ", ld_file)
system(cmd)
ld=fread(paste0(ld_file, ".ld"), header=T) %>% as.data.frame


#### Get lead SNP from data file
arrange(mydata, P) %>% head(n=1) %>% pull(SNP) -> lead_snp


### subset ld matrix by lead_snp
filter(ld, SNP_A == lead_snp ) %>% select(SNP_B, R2) -> subset_ld
names(subset_ld)=c("SNP", "LD")
subset_ld$LD=round(subset_ld$LD, 1)
subset_ld=rbind(subset_ld, data.frame(SNP=lead_snp, LD="1"))

### merge d plotata
merge(merge(mydata, subset_ld, all.x=T, by = "SNP"), LD.colours, by="LD", all.x=T) -> mydata
mydata$LD=ifelse(is.na(mydata$LD), "0",mydata$LD)

null_color=subset(LD.colours, LD==0) %>% pull(Colour) %>% as.character
mydata$Colour=ifelse(is.na(mydata$Colour), null_color,mydata$Colour)


### make plot
ld_colors = as.vector(mydata$Colour)
names(ld_colors) = mydata$LD

mydata$LD=factor(mydata$LD, levels=c(sort(unique(mydata$LD))))

p=ggplot(mydata)+
  geom_point(mapping=aes(x = BP/1000000, y = -log10(P), color=LD), size=2)+
  scale_color_manual(values=c(ld_colors))+
  xlab("Position (Mb)")+
  ylab(expression("-log"[10]*italic(P)))+
  theme_classic(base_size = 20, base_line_size = 1)+
  theme(legend.title=element_blank(),
          axis.text.x = element_text(color="black", size=20),
          axis.text.y = element_text(color="black", size=20),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(1, 1, 1, 1))+
  guides(color = guide_legend(override.aes=list(shape = 15)))
ggsave("Example.pdf",p, width=10, height=3)






