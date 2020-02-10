# Box plots of all conditions (for example patient samples (orig.ident)) for each major cell type (groups_bysize)

################################################################ PULP ###############################################################
    
    data_df         = cbind.data.frame(value= rep(1,length(all_healthy@meta.data$groups_bysize)), groups_bysize = all_healthy@meta.data$groups_bysize, condition = all_healthy@meta.data$orig.ident)

    df_summ         = aggregate( data_df$value, list( data_df$groups_bysize, data_df$condition), sum)

    names(df_summ)  = c('majorCT', 'Patient', 'counts')

    g  = ggplot(data = df_summ, aes(x=majorCT, y=counts)) + 
         geom_boxplot(outlier.size=0) +
         geom_point( size = 5, aes(color = Patient)) + 
         theme(axis.line = element_line(color="black", size = .75), axis.text=element_text(size=20,face="bold",angle = 45, hjust = 1),
         axis.title.x=element_blank(), axis.title.y=element_text(size=20,face="bold", vjust = 1))+
         theme(legend.text=element_text(size=22),legend.title=element_text(size=24))+
         labs(y= expression(paste('Number of Cells')))
         # +
         # guides(shape = guide_legend(override.aes = list(size = 5)))

   ggsave(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Box_byMCT.pdf"), width=14, height=10)

# LOG
    g  = ggplot(data = df_summ, aes(x=majorCT, y=counts)) + 
         geom_boxplot(outlier.size=0) +
         geom_point( size = 5, aes(color = Patient)) + 
        scale_y_continuous(trans='log10') +
         theme(axis.line = element_line(color="black", size = .75), axis.text=element_text(size=20,face="bold",angle = 45, hjust = 1),
         axis.title.x=element_blank(), axis.title.y=element_text(size=20,face="bold", vjust = 1))+
         theme(legend.text=element_text(size=22),legend.title=element_text(size=24))+
         labs(y= expression(paste('Number of Cells (','log'['10'],')')))
         # +
         # guides(shape = guide_legend(override.aes = list(size = 5)))

   ggsave(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPulp", "Box_byMCT_log.pdf"), width=14, height=10)


############################################################### PERIO ###############################################################

    data_df         = cbind.data.frame(value= rep(1,length(all_perio@meta.data$groups_bysize)), groups_bysize = all_perio@meta.data$groups_bysize, condition = all_perio@meta.data$orig.ident)

    df_summ         = aggregate(data_df$value, list(data_df$groups_bysize, data_df$condition), sum)

    names(df_summ)  = c('majorCT', 'Patient', 'counts')

    g  = ggplot(data = df_summ, aes(x=majorCT, y=counts)) + 
         geom_boxplot(outlier.size=0) +
         geom_point(aes(color = Patient), size = 5) + 
         theme(axis.line = element_line(color="black", size = .75), axis.text=element_text(size=20,face="bold",angle = 45, hjust = 1),
         axis.title.x=element_blank(), axis.title.y=element_text(size=20,face="bold", vjust = 1))+
         theme(legend.text=element_text(size=22),legend.title=element_text(size=24))+
         labs(y= expression(paste('Number of Cells')))
         # +
         # guides(shape = guide_legend(override.aes = list(size = 5)))

   ggsave(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Box_byMCT.pdf"), width=14, height=10)

# LOG
    g  = ggplot(data = df_summ, aes(x=majorCT, y=counts)) + 
         geom_boxplot(outlier.size=0) +
         geom_point( size = 5, aes(color = Patient)) + 
        scale_y_continuous(trans='log10') +
         theme(axis.line = element_line(color="black", size = .75), axis.text=element_text(size=20,face="bold",angle = 45, hjust = 1),
         axis.title.x=element_blank(), axis.title.y=element_text(size=20,face="bold", vjust = 1))+
         theme(legend.text=element_text(size=22),legend.title=element_text(size=24))+
         labs(y= expression(paste('Number of Cells (','log'['10'],')')))
         # +
         # guides(shape = guide_legend(override.aes = list(size = 5)))

   ggsave(file = file.path("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/AllPerio", "Box_byMCT_log.pdf"), width=14, height=10)

############################################################### PERIO/PULP Proportions ###############################################################
melted_celltype_condition2         = cbind.data.frame(cellType =merged_harmony@meta.data$groups_bysize, condition=merged_harmony@meta.data$condition)
melted_celltype_by_condition2      = as.data.frame.matrix(table(melted_celltype_condition2))
cellType_normalized_by_condition2  = as.data.frame(t(melted_celltype_by_condition2)/rowSums(t(melted_celltype_by_condition2)))

datm3 <- cellType_normalized_by_condition2 %>%
mutate(ind = factor(row_number())) %>%
gather(variable, value, -ind)

ggplot(datm3, aes(x = variable, y = value, fill = ind)) +
geom_bar(position = position_fill(), stat = "identity")+
scale_y_continuous(labels = scales::percent_format())+
labs(x="Cell type", y="Proportion of sample",
title="Proportion of cell type from sample (normalized)", fill = "Sample") +
scale_fill_manual(labels = c("Perio", "Pulp"), values = c("blue", "red")) +
theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/Normalized_proportions_perio_pulp.pdf")



