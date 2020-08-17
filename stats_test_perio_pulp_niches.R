# wilcoxon ranked sum test for proportion of cell types across perio/pulp for 5-patient groups. 
# Statistically tricky since: 1. Small sample sizes and 2. some samples are "dependent" (perio/pulp from same patient...)
library(RColorBrewer)

load("/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/r_output/Review.Rdata")

outputsDir 		=  "/IMCR_shares/Moorlab/Common/Tooth_project/R_analysis/ldvr_analyses/"

colorPalette 	= colorRampPalette(brewer.pal(8, 'Set1'))

pd 				= position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)

############################################# 
# Generate condition meta.data from orig.ident if not already present in Seurat object
# red_orig.ident 		 				= gsub('[[:digit:]]+', '', merged_harmony@meta.data$orig.ident)
# red_orig.ident2 	 					= gsub('v', '', red_orig.odent)
# merged_harmony@meta.data$condition  = red_orig.ident2

  	cellType 							= merged_harmony@meta.data$orig.ident
  	majorCT 							= merged_harmony@meta.data$groups_bysize
  	cellCategories 						= merged_harmony@meta.data$condition

    metas_test                          = table(cbind.data.frame(cellType, majorCT))

    patient_and_label                   = unique(cbind.data.frame(cellType, cellCategories))

    metas_patient_and_label             = merge( metas_test, patient_and_label, by="cellType") 

    names(metas_patient_and_label)      = c('cellType', 'majorCT',  "Freq" ,'cellCategory' )

	# Normalization is relative to file size:

    files_tot_pop                    = aggregate(metas_patient_and_label$Freq, list(metas_patient_and_label$cellType), sum)
    names(files_tot_pop)             = c('cellType','tot_pop')
    metas_summ_tot_pop               = merge(metas_patient_and_label, files_tot_pop, by="cellType") 
    metas_summ_tot_pop$relativ_freq  = metas_summ_tot_pop$Freq/metas_summ_tot_pop$tot_pop

    # metas_colors                     = colorPalette(length(unique(merged_harmony@meta.data$condition))) # random color squeme
     metas_colors                     = c(perio = 'blue', pulp='red')

    # # First plot proportions in Perio/pulp as adjacent boxplots - NORMALIZED BY SAMPLE SIZE
    g  = ggplot(data = metas_summ_tot_pop, aes(x=majorCT, y=relativ_freq, fill=cellCategory)) + 
         geom_boxplot(outlier.size=0) +
         scale_fill_manual(values = metas_colors)+
         geom_point(aes(color = cellCategory), position = pd, show.legend = FALSE) +
         scale_color_manual(values = metas_colors)+
         # scale_y_continuous(trans='log10') +
         labs(fill='Sample type')+ 
         theme(axis.line = element_line(color="black", size = .75), axis.text=element_text(size=20, face="bold", angle = 45, hjust = 1),
         axis.title.x=element_blank(), axis.title.y=element_text(size=22,face="bold"))+
         theme(legend.text=element_text(size=22),legend.title=element_text(size=24))+ 
         labs(y= 'Frequency in sample')

         # labs(y= expression(paste('Frequency of major cell type (','log'['10'],')')))
	sample_name = 'tooth'
    ggsave(file.path(outputsDir, paste0('RelativMajorCellTypesby_', sample_name, '_boxplotlin.png')), width = 22, height = 10)

	# Wilcoxon rank sum test
	NmajorCT                         = sort(unique(metas_summ_tot_pop$majorCT))
    cellCategories_pairs             = combn(unique(metas_summ_tot_pop$cellCategory),2)

    if (length(unique(metas_summ_tot_pop$cellCategory))==2){
            ncc                              = 1
    }else{
            ncc                              = 1:ncol(cellCategories_pairs) 
    }

	 for (cc in ncc) {

        if (length(unique(metas_summ_tot_pop$cellCategory))==2){
            temp_cc          = as.character(cellCategories_pairs)
        }else{
            temp_cc          = cellCategories_pairs[,cc]
        }

        temp_subset1     = subset(metas_summ_tot_pop, cellCategory == temp_cc[1], select = c('majorCT','relativ_freq'))
        temp_subset2     = subset(metas_summ_tot_pop, cellCategory == temp_cc[2], select = c('majorCT','relativ_freq'))

        # This is the rank sum test, not signed rank! wilcoxon on samples for which there were more than zero cells from a given cluster. This excludes zeros. This is not paired.
        # temp_wilcoxon    = lapply(1:length(NmajorCT), function(x) if( length(subset(temp_subset1, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq) > 0 & length(subset(temp_subset2, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq) > 0 )
        # wilcox.test(subset(temp_subset1, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq,  subset(temp_subset2, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq )     )

        # temp_wilcoxon    = lapply(1:length(NmajorCT), function(x) if( length(subset(temp_subset1, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq) > 0 & length(subset(temp_subset2, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq) > 0 )
        #     wilcox.test(subset(temp_subset1, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq,  subset(temp_subset2, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq )     )
 
        ###################### NOW INCLUDE PATIENT SAMPLES WHERE THERE WERE ZEROS AND PERFORM A PAIRED WILCOXON TEST  ######################
        temp_wilcoxon    = lapply(1:length(NmajorCT), function(x) wilcox.test(subset(temp_subset1, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq,  subset(temp_subset2, majorCT == NmajorCT[x], select = c('relativ_freq'))$relativ_freq, paired = TRUE,)     )

        # Extract p-vals for each major cell type to print 
        temp_ps         = lapply(1:length(temp_wilcoxon), function(x) temp_wilcoxon[[x]]$p.value)

        adjusted_ps     = lapply(1:length(temp_wilcoxon), function(x) p.adjust(temp_wilcoxon[[x]]$p.value, method = 'bonferroni', n = length(temp_ps)))

    
        write.table(cbind(as.character(NmajorCT), temp_ps) , file.path(outputsDir,paste0(sample_name,'wilcoxon_',paste0(temp_cc, collapse='_'),'.csv')),
                                    col.names=TRUE, row.names=FALSE, sep=", ")

        # write.table(cbind(as.character(NmajorCT), adjusted_ps ) , file.path(outputsDir,paste0(sample_name,'p_adjusted_BH_wilcoxon_',paste0(temp_cc, collapse='_'),'.csv')),
        #                             col.names=TRUE, row.names=FALSE, sep=", ")

        write.table(cbind(as.character(NmajorCT), adjusted_ps ) , file.path(outputsDir,paste0(sample_name,'p_adjusted_Bonferroni_wilcoxon_',paste0(temp_cc, collapse='_'),'.csv')),
                                    col.names=TRUE, row.names=FALSE, sep=", ")


    }


