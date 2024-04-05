require(ggplot2);require(reshape2);require(scales);require(ggpubr);require(tidyr)

s = read.csv('qr-optimal_astral2_estg_ests.csv')
head(s)
s$method =  factor(s$method) 
levels(s$method) = list("QR" = "QR", 
                        "QR-STAR" = "QR*", 
                        "Optimal" = "optimal")

s$height =  factor(s$height) 
levels(s$height) = list("500K" = "500000",
                        "2M" = "2000000", 
                        "10M" = "10000000")

ggplot(s, aes(x=as.factor(k), y=cd/396, fill = method))+
  geom_boxplot()+
  facet_grid(speciation~height)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_brewer(palette = "Greens")+
  scale_color_brewer(palette = "Greens",name="")+
  scale_x_discrete(name="Number of Genes")+
  coord_cartesian(clip="off",ylim=c(0,0.4))+
  scale_y_continuous(name="Rooted Species Tree Error (nCD)")+
  theme_bw()
ggsave("qrstar-ests-estg.pdf",width=10,height = 4)


q = read.csv('qr-avian-sim-10.csv')
head(q)
unique(q$rooted.tree.name)
q$rooted.tree.name =  factor(q$rooted.tree.name, levels=c("qr_exh.1.2.4.astraltree", "qr_le.1.2.4.astraltree", 
                                                          "mad.2.2-RAxML_result.astraltree.rooted", "rd.1.7.0-RAxML_result.astraltree", 
                                                          "MP.1.5-RAxML_result.astraltree", "MV.1.5-RAxML_result.astraltree")) 
levels(q$rooted.tree.name) = list("QR (Exh)" = "qr_exh.1.2.4.astraltree", 
                        "QR (LE)" = "qr_le.1.2.4.astraltree", 
                        "MAD" = "mad.2.2-RAxML_result.astraltree.rooted",
                        "Root Digger" = "rd.1.7.0-RAxML_result.astraltree",
                        "Midpoint" = "MP.1.5-RAxML_result.astraltree",
                        "Minvar" = "MV.1.5-RAxML_result.astraltree")
q$seq.length =  factor(q$seq.length) 
levels(q$seq.length) = list("true gene trees" = "true", 
                                  "1500bp" = "1500", 
                                  "1000bp" = "1000",
                                  "500bp" = "500",
                                  "250bp" = "250")


ggplot(data=q, aes(x=seq.length, y=cd/16, color=rooted.tree.name))+
  #geom_boxplot()+
  #stat_summary(position = position_dodge(width=0.86))+
  #stat_summary(geom="line")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(aes(group=rooted.tree.name),geom="line")+
  #stat_summary(geom="line")+
  stat_summary()+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_y_continuous(name="Rooted Species Tree Error (nCD)")+
  guides(fill=guide_legend(title="Rooting Method"))+
  theme_bw()+
  theme(legend.position = "bottom", legend.direction = "horizontal")+
  scale_x_discrete(name="Condition")
ggsave("qr-seq-len.pdf",width=5,height = 4)

q = read.csv('qr-avian-sim-n.csv')
head(q)
unique(q$rooted.tree.name)
q$rooted.tree.name =  factor(q$rooted.tree.name, levels=c("qr_exh.1.2.4.modeltree", "qr_le.1.2.4.modeltree", 
                                                          "mad.2.2-RAxML_result.modeltree.rooted", "rd.1.7.0-RAxML_result.modeltree", 
                                                          "MP.1.5-RAxML_result.modeltree", "MV.1.5-RAxML_result.modeltree")) 
levels(q$rooted.tree.name) = list("QR (Exh)" = "qr_exh.1.2.4.modeltree", 
                                  "QR (LE)" = "qr_le.1.2.4.modeltree", 
                                  "MAD" = "mad.2.2-RAxML_result.modeltree.rooted",
                                  "Root Digger" = "rd.1.7.0-RAxML_result.modeltree",
                                  "Midpoint" = "MP.1.5-RAxML_result.modeltree",
                                  "Minvar" = "MV.1.5-RAxML_result.modeltree")

ggplot(data=q, aes(x=as.factor(n), y=cd/(2*n-4), color=rooted.tree.name))+
  #geom_boxplot()+
  #stat_summary(position = position_dodge(width=0.86))+
  #stat_summary(geom="line")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(aes(group=rooted.tree.name),geom="line")+
  #stat_summary(geom="line")+
  stat_summary()+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_y_continuous(name="Rooting Error (nCD)")+
  guides(fill=guide_legend(title="Rooting Method"))+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  scale_x_discrete(name="Number of Species")
ggsave("qr-num-species.pdf",width=5,height = 3.8)


ggplot(aes(x=as.factor(n), y=cd/(2*n-4), color=rooted.tree.name),
       data=dcast(data=q,rooted.tree.name+n+replicate~'cd' ,value.var = "cd",fun.aggregate = mean))+
  geom_boxplot(outlier.alpha = 0.3,width=0.86)+
  stat_summary(position = position_dodge(width=0.86))+
  #stat_summary(geom="line")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  #stat_summary(aes(group=rooted.tree.name),geom="line")+
  #stat_summary(geom="line")+
  #stat_summary()+
  scale_fill_brewer(palette = "Dark2")+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_y_continuous(name="Rooting Error (nCD)")+
  guides(fill=guide_legend(title="Rooting Method"))+
  theme_bw()+
  theme(legend.position = "none", legend.direction = "horizontal")+
  scale_x_discrete(name="Number of Species")
ggsave("qr-num-species.pdf",width=5,height = 4)

