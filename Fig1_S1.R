rm(list = ls()) ; graphics.off()
#------------------------------------------------------------------------------------#
# package prepare
#------------------------------------------------------------------------------------#
library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(export)
library(RColorBrewer)
library(plyr)
library(survival)
#library(survminer)
load("/data/Fig1_S1_data.Rdata")
load("/data/CBCGA.Extended_MergedData_V2.5_220722.Rdata")
#------------------------------------------------------------------------------------#
# Fig1
#------------------------------------------------------------------------------------#
#### 1. data prepare ------------
## clinical 
CBCGAClin = CBCGA.Extended_Cohort.Info[CBCGA.Extended_Cohort.Info$CBCGA_Cohort == "Yes",]
CBCGAClin$TMB = TMB[rownames(CBCGAClin)]

CBCGAClin$Lymph_node_status = LN[rownames(CBCGAClin)]
CBCGAClin$Lymph_node_status[CBCGAClin$Lymph_node_status > 0] = 1

CBCGAClin$Ki67 = Ki67[rownames(CBCGAClin)]

CBCGAClin = CBCGAClin[colnames(CBCGA_merge_mat),]

## DNA --- CBCGA_merge_mat

## RNA
exp.fpkm.TT = CBCGA.Extended_RNA_913_FPKM_symbol_proteinCode # havenot already log2-transformed
exp.fpkm.TT = exp.fpkm.TT[,str_detect(colnames(exp.fpkm.TT),"_T")]
colnames(exp.fpkm.TT) = substring(colnames(exp.fpkm.TT),1,4)
tmp = as.data.frame( matrix(nrow = nrow(exp.fpkm.TT), 
                            ncol = length(setdiff(sample_ord, colnames(exp.fpkm.TT)) ) ) )
rownames(tmp) = rownames(exp.fpkm.TT) ; colnames(tmp) = setdiff(sample_ord, colnames(exp.fpkm.TT))

exp.fpkm.TT = exp.fpkm.TT[, intersect(sample_ord , colnames(exp.fpkm.TT) )]
exp.fpkm.TT.addNA.log = log2( cbind(exp.fpkm.TT,tmp) +1 )

## Protein
protein_TT_log2 = CBCGA.Extended_PRO_normalized # have already scaled and log2-transformed
protein_TT_log2 = protein_TT_log2[,str_detect(colnames(protein_TT_log2),"_T")]
colnames(protein_TT_log2) = substring(colnames(protein_TT_log2),1,4)

tmp = as.data.frame( matrix(nrow = nrow(protein_TT_log2), 
                            ncol = length(setdiff(sample_ord, colnames(protein_TT_log2)) ) ) )
rownames(tmp) = rownames(protein_TT_log2) ; colnames(tmp) = setdiff(sample_ord, colnames(protein_TT_log2))
protein_TT_log2 =  protein_TT_log2[, intersect(sample_ord , colnames(protein_TT_log2) )]
protein_TT_log2_addNA = cbind(protein_TT_log2,tmp) 

## polar
polar_TT_log2 = CBCGA.Extended_pol
# polar_TT_log2$class = CBCGA_pol_anno[rownames(polar_TT_log2),"Metabolite_class"]
# polar_TT_log2 = aggregate(polar_TT_log2[1:(ncol(polar_TT_log2)-1)],list(polar_TT_log2$class),mean)
# rownames(polar_TT_log2) = polar_TT_log2$Group.1  ; polar_TT_log2 = polar_TT_log2[,-1]
# polar_TT_log2 = polar_TT_log2[c( "Amino acid","Carbohydrates","Nucleotide","Peptide","Vitamins and Cofactors","Xenobiotics" ) , ]
polar_TT_log2 = polar_TT_log2[ , str_detect(colnames(polar_TT_log2),"_T") ]
colnames(polar_TT_log2) = substring(colnames(polar_TT_log2),1,4)

tmp = as.data.frame( matrix(nrow = nrow(polar_TT_log2), 
                            ncol = length(setdiff(sample_ord, colnames(polar_TT_log2)) ) ) )
rownames(tmp) = rownames(polar_TT_log2) ; colnames(tmp) = setdiff(sample_ord, colnames(polar_TT_log2))
polar_TT_log2 =  polar_TT_log2[, intersect(sample_ord , colnames(polar_TT_log2) )]
polar_TT_log2_addNA = cbind(polar_TT_log2,tmp) 

## lipid
lipid_TT_log2 = CBCGA.Extended_lip
# lipid_TT_log2$class = CBCGA_lip_anno[rownames(lipid_TT_log2),"Lipid.super.class"]
# lipid_TT_log2 = aggregate(lipid_TT_log2[1:(ncol(lipid_TT_log2)-1)],list(lipid_TT_log2$class),mean)
# rownames(lipid_TT_log2) = lipid_TT_log2$Group.1  ; lipid_TT_log2 = lipid_TT_log2[,-1]
lipid_TT_log2 = lipid_TT_log2[ , str_detect(colnames(lipid_TT_log2),"_T") ]
colnames(lipid_TT_log2) = substring(colnames(lipid_TT_log2),1,4)

tmp = as.data.frame( matrix(nrow = nrow(lipid_TT_log2), 
                            ncol = length(setdiff(sample_ord, colnames(lipid_TT_log2)) ) ) )
rownames(tmp) = rownames(lipid_TT_log2) ; colnames(tmp) = setdiff(sample_ord, colnames(lipid_TT_log2))
lipid_TT_log2 =  lipid_TT_log2[, intersect(sample_ord , colnames(lipid_TT_log2) )]
lipid_TT_log2_addNA = cbind(lipid_TT_log2,tmp) 

##
CBCGARNA_mat <- exp.fpkm.TT.addNA.log[RNA4plot, ]
CBCGARNA_mat <- CBCGARNA_mat[, colnames(CBCGA_merge_mat)]
write.csv(CBCGARNA_mat,file = "/results/TableS1F.csv")
CBCGARNA_mat <- as.matrix(CBCGARNA_mat)

CBCGARNA_mat_colname <- colnames(CBCGARNA_mat)
CBCGARNA_mat <- apply(CBCGARNA_mat, 1, scale) %>% t()
colnames(CBCGARNA_mat) <- CBCGARNA_mat_colname

CBCGAProtein_mat <- protein_TT_log2_addNA[Protein4plot, ]
CBCGAProtein_mat <- CBCGAProtein_mat[, colnames(CBCGA_merge_mat)]
write.csv(CBCGAProtein_mat,file = "/results/TableS1G.csv")
CBCGAProtein_mat <- as.matrix(CBCGAProtein_mat)

CBCGAPolar_mat <- polar_TT_log2_addNA[ picked.Polar[order(names(picked.Polar))] , colnames(CBCGA_merge_mat)]
write.csv(CBCGAPolar_mat,file = "/results/TableS1H.csv")
CBCGAPolar_mat_colname <- colnames(CBCGAPolar_mat)
CBCGAPolar_mat <- apply(CBCGAPolar_mat, 1, scale) %>% t()
colnames(CBCGAPolar_mat) <- CBCGAPolar_mat_colname
CBCGAPolar_mat <- as.matrix(CBCGAPolar_mat)

CBCGALipid_mat <- lipid_TT_log2_addNA[ picked.Lipid[order(names(picked.Lipid))], colnames(CBCGA_merge_mat)]
write.csv(CBCGALipid_mat,file = "/results/TableS1I.csv")
CBCGALipid_mat_colname <- colnames(CBCGALipid_mat)
CBCGALipid_mat <- apply(CBCGALipid_mat, 1, scale) %>% t()
colnames(CBCGALipid_mat) <- CBCGALipid_mat_colname
CBCGALipid_mat <- as.matrix(CBCGALipid_mat)


#### 2. parameter prepare -----
col = c("Splicing" = "#6b4694", "Missense" = "#4d93cd", "Nonsense" = "#e63a3a", "Frameshift" = '#fab71b', "Inframe" = "#ecd71e",
        "Nonstop" = '#ab5b9e', "Translation_Start_Site" = '#018b38', "No" = '#eaeaea', 'Amplification' = '#BC102B', 'Gain' = '#DF989E',
        'Loss' = '#AEC4D6', 'Deletion' = '#5385AC')
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # small purple
  Splicing = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.67,
              gp = gpar(fill = col["Splicing"], col = NA))
  },
  # small green
  Missense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
              gp = gpar(fill = col["Missense"], col = NA))
  },
  # small black
  Nonsense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
              gp = gpar(fill = col["Nonsense"], col = NA))
  },
  # small black
  Frameshift = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
              gp = gpar(fill = col["Frameshift"], col = NA))
  },
  # small black
  Inframe = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
              gp = gpar(fill = col["Inframe"], col = NA))
  },
  # small black
  Nonstop = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
              gp = gpar(fill = col["Nonstop"], col = NA))
  },
  # small black
  Translation_Start_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h*0.33,
              gp = gpar(fill = col["Translation_Start_Site"], col = NA))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
              gp = gpar(fill = col["Amplification"], col = NA))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
              gp = gpar(fill = col["Gain"], col = NA))
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
              gp = gpar(fill = col["Loss"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
              gp = gpar(fill = col["Deletion"], col = NA))
  },
  No = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1.5, "pt"), h-unit(2.5, "pt"),
              gp = gpar(fill = col["No"], col = NA))
  }
)

column_title <- ""
heatmap_legend_param <- list(title = "", at = c("Splicing", "Missense", 'Nonsense', 'Frameshift', 'Inframe', "Nonstop", "Translation_Start_Site",
                                                'Amplification', 'Gain', 'Loss', 'Deletion', 'No'),
                             labels = c("Splicing", "Missense", 'Nonsense', 'Frameshift', 'Inframe', "Nonstop", "Translation_Start_Site",
                                        'Amplification', 'Gain', 'Loss', 'Deletion', 'No'))
col_fun_age <- colorRamp2(c(24, 53, 90), c('#ecf3e4', '#9ec27d', '#489205'))
col_fun_ki67 <- colorRamp2(c(3, 30, 98), c('#EBEFF6', '#99AFD2', '#3C74AE'))
col_fun_hrd <- colorRamp2(c(0, 16, 75), c('#FDF9E7', '#EFE089', '#D7C801'))

TMB1 <- as.numeric(CBCGAClin$TMB)
TMB1[which(TMB1 > 6)] <- NA
TMB2 <- as.numeric(CBCGAClin$TMB)
TMB2[which(TMB2 <= 6)] <- NA

pos <- 5 + 0.1
neg <- -5 + 0.1
poscut <- 50
negcut <- 50
mybreaks1 <- -(seq(-sqrt(-neg), 0, sqrt(-neg)/negcut)) ^2
mybreaks2 <- (seq(0, sqrt(pos), sqrt(pos)/poscut)[-1]) ^2
mybreaks  <- c(mybreaks1, mybreaks2)

all(rownames(picked.Pro) == rownames(CBCGAProtein_mat))
all(rownames(picked.RNA) == rownames(CBCGARNA_mat))
picked.RNA$path <- factor(picked.RNA$path, levels = c('ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE', 'ERBB2_SIGNALING_PATHWAY',
                                                      'SIGNALING_BY_ERBB2_ECD_MUTANTS', 'NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                                                      'POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE', 'T_CELL_RECEPTOR_SIGNALING_PATHWAY',
                                                      'MYC_TARGETS', 'CELL_CYCLE', 'P53_SIGNALING_PATHWAY', 'NEGATIVE_REGULATION_OF_DNA_REPAIR'))
picked.RNA$name <- rownames(picked.RNA)
picked.RNA <- picked.RNA[order(picked.RNA$path), ]
picked.RNA <- picked.RNA[, -2, drop = F]
picked.Pro$path <- factor(picked.Pro$path, levels = c('ESTROGEN_RESPONSE_EARLY', 'ESTROGEN_RESPONSE_LATE', 'ERBB2_SIGNALING_PATHWAY',
                                                      'SIGNALING_BY_ERBB2_ECD_MUTANTS', 'NATURAL_KILLER_CELL_MEDIATED_IMMUNITY',
                                                      'POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE', 'T_CELL_RECEPTOR_SIGNALING_PATHWAY',
                                                      'MYC_TARGETS', 'CELL_CYCLE', 'P53_SIGNALING_PATHWAY', 'NEGATIVE_REGULATION_OF_DNA_REPAIR'))
picked.Pro$name <- rownames(picked.Pro)
picked.Pro <- picked.Pro[order(picked.Pro$path), ]
picked.Pro <- picked.Pro[, -2, drop = F]

right_annotation_RNA <- rowAnnotation(df = picked.RNA,
                                      show_annotation_name = F,
                                      col = list(path = c("T_CELL_RECEPTOR_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[1],
                                                          "POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE" = brewer.pal(12, "Set3")[2] ,
                                                          "NATURAL_KILLER_CELL_MEDIATED_IMMUNITY" = brewer.pal(12, "Set3")[3],
                                                          "INNATE_IMMUNE_RESPONSE_ACTIVATING_SIGNAL_TRANSDUCTION" = brewer.pal(12, "Set3")[4],
                                                          "SIGNALING_BY_ERBB2_ECD_MUTANTS" = brewer.pal(12, "Set3")[5],
                                                          "ERBB2_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[6],
                                                          "ESTROGEN_RESPONSE_EARLY" = brewer.pal(12, "Set3")[7],
                                                          "ESTROGEN_RESPONSE_LATE" = brewer.pal(12, "Set3")[8],
                                                          "MYC_TARGETS" = brewer.pal(12, "Set3")[9],
                                                          "CELL_CYCLE" = brewer.pal(12, "Set3")[10],
                                                          "P53_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[11],
                                                          "NEGATIVE_REGULATION_OF_DNA_REPAIR" = brewer.pal(12, "Set3")[12]))
)
right_annotation_Pro <- rowAnnotation(df = picked.Pro,
                                      show_annotation_name = F,
                                      col = list(path = c("T_CELL_RECEPTOR_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[1],
                                                          "POSITIVE_REGULATION_OF_INNATE_IMMUNE_RESPONSE" = brewer.pal(12, "Set3")[2] ,
                                                          "NATURAL_KILLER_CELL_MEDIATED_IMMUNITY" = brewer.pal(12, "Set3")[3],
                                                          "INNATE_IMMUNE_RESPONSE_ACTIVATING_SIGNAL_TRANSDUCTION" = brewer.pal(12, "Set3")[4],
                                                          "SIGNALING_BY_ERBB2_ECD_MUTANTS" = brewer.pal(12, "Set3")[5],
                                                          "ERBB2_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[6],
                                                          "ESTROGEN_RESPONSE_EARLY" = brewer.pal(12, "Set3")[7],
                                                          "ESTROGEN_RESPONSE_LATE" = brewer.pal(12, "Set3")[8],
                                                          "MYC_TARGETS" = brewer.pal(12, "Set3")[9],
                                                          "CELL_CYCLE" = brewer.pal(12, "Set3")[10],
                                                          "P53_SIGNALING_PATHWAY" = brewer.pal(12, "Set3")[11],
                                                          "NEGATIVE_REGULATION_OF_DNA_REPAIR" = brewer.pal(12, "Set3")[12]))
                                      
)

CBCGARNA_mat <- CBCGARNA_mat[rownames(picked.RNA), ]
CBCGAProtein_mat <- CBCGAProtein_mat[rownames(picked.Pro), ]

tmp = data.frame(row.names = c(rownames(CBCGAPolar_mat)),
                 class = c(rep("Polar metabolite",nrow(CBCGAPolar_mat))))
right_annotation_polar = rowAnnotation(df = tmp,
                                       show_annotation_name = F,
                                       col = list(class = c("Polar metabolite" = "#FFC000"))
                                       
)

tmp = data.frame(row.names = c(rownames(CBCGALipid_mat)),
                 class = c(rep("Lipids",nrow(CBCGALipid_mat))))
right_annotation_lipid = rowAnnotation(df = tmp,
                                       show_annotation_name = F,
                                       col = list(class = c("Lipids" = "#4472C6"))
                                       
)


## MutSig
mutsig_ord = c("SBS5","SBS1",
               "SBS13","SBS2",
               "SBS8","SBS18","SBS17a","SBS17b",
               "SBS20","SBS26","SBS30","SBS3","SBS6",
               "SBSNA")
mutsig_col = c('SBS5' ='#0084C8' , 'SBS1' = '#A6CEE3', 
               'SBS13' =  '#009100' ,'SBS2' = '#9ADE00' ,
               'SBS8' ='#6A3D9A'  ,
               'SBS17a' ='#33C6BB', 'SBS17b' =  '#99E3DD', 'SBS18' = '#008A80',
               'SBS3' = '#FF6600','SBS6' = '#DC0000', 'SBS20' = '#FFFF3E','SBS26' = '#FFC022','SBS30' = '#FF9900', 
               'SBSNA' = '#eaeaea')

MutSig = MutSig[mutsig_ord,]
mutsig_col = mutsig_col[mutsig_ord]


#### draw -------
CBCGA_landscape <- oncoPrint(CBCGA_merge_mat, alter_fun = alter_fun, col = col, show_pct = T, pct_side = 'right',
                             column_title = '', heatmap_legend_param = heatmap_legend_param, column_order = sample_ord,
                             row_names_side = "left", show_column_names = F, remove_empty_columns = F,
                             row_names_gp = gpar(fontsize = 9), column_title_gp = gpar(fontsize = 0),
                             
                             top_annotation = HeatmapAnnotation(TMB2 = anno_points(TMB2, size = unit(1, 'mm'), border = T, ylim = c(6, 20), height = unit(0.75, 'cm'),
                                                                                   axis_param = list(at = seq(8, 20, by = 6), labels = seq(8, 20, by = 6))),
                                                                TMB1 = anno_points(TMB1, size = unit(1, 'mm'), border = T, ylim = c(0, 5.95), height = unit(1.5, 'cm'),
                                                                                   axis_param = list(at = c(0, 2, 4, 6), labels = c(0, 2, 4, 6))),
                                                                IHC_Subtype = CBCGAClin$Clinical_Subtype,
                                                                Intrinsic_subtype = CBCGAClin$PAM50_classifier,
                                                                Age.at.surgery = CBCGAClin$Age,
                                                                Lymph_node_status = CBCGAClin$Lymph_node_status,
                                                                Ki67 = as.numeric(CBCGAClin$Ki67),
                                                                HRD = as.numeric(CBCGAClin$HRD),
                                                                Relapse = factor(as.character(CBCGAClin$RFS_status)),
                                                                COSMIC_sig = anno_barplot(t(MutSig), height = unit(1.5, 'cm'),
                                                                                          gp = gpar(fill =mutsig_col,  col = NA)),
                                                                # COSMIC_sig = anno_barplot(t(MutSig), height = unit(1.5, 'cm'),
                                                                #                            gp = gpar(fill = mutsig_col,col = NA)),
                                                                annotation_name_side = "left",
                                                                annotation_name_gp = gpar(fontsize = 9),
                                                                col = list(IHC_Subtype = c('HR+HER2-' = '#0085c4', 'HR+HER2+' = '#7ab801',
                                                                                           'HR-HER2+' = '#f2af01', 'TNBC' = '#dc5035'),
                                                                           Intrinsic_subtype = c('LumA' = '#0077c9', 'LumB' = '#74d2e8',
                                                                                                 'Her2' = '#7552cd', 'Basal' = '#e4002c',
                                                                                                 'Normal' = '#cecece',"Unknown" = "grey"),
                                                                           Age.at.surgery = col_fun_age,
                                                                           Lymph_node_status = c('0' = 'white', '1' = 'black'),
                                                                           Ki67 = col_fun_ki67,
                                                                           HRD = col_fun_hrd,
                                                                           Relapse = c('0' = 'white', '1' = 'black','NA' = "grey"))),
                             column_split = factor(CBCGAClin$Clinical_Subtype, levels = c('HR+HER2-', 'HR+HER2+', 'HR-HER2+', 'TNBC')),
                             row_split = factor(c(rep('Somatic', 13), rep('Germline', 5), rep('AMP', 9), rep('DEL', 8)), 
                                                levels = c('Somatic', 'Germline', 'AMP', 'DEL')),
                             gap = unit(0.1, 'inch'))  %v%
  Heatmap(CBCGARNA_mat, show_row_names = F, show_column_names = F, height = unit(3, 'cm'), cluster_columns = F, cluster_rows = F, na_col = '#eaeaea',
          col = colorRamp2(mybreaks, c(colorRampPalette(c("green", "#000000"))(51), colorRampPalette(c("#000000", "red"))(50))), right_annotation = right_annotation_RNA) %v%
  Heatmap(CBCGAProtein_mat, show_row_names = F, show_column_names = F, height = unit(3, 'cm'), cluster_columns = F, cluster_rows = F, na_col = '#eaeaea',
          col = colorRamp2(mybreaks, c(colorRampPalette(c("#3A53A4", "#FFFFFF"))(51), colorRampPalette(c( "#FFFFFF", "#E21E26"))(50))), right_annotation = right_annotation_Pro) %v%
  
  Heatmap(CBCGAPolar_mat, show_row_names = F, show_column_names = F, height = unit(1, 'cm'), cluster_columns = F, cluster_rows = F, na_col = '#eaeaea',
          show_row_dend = F,right_annotation = right_annotation_polar,
          col = colorRamp2(mybreaks, c(colorRampPalette(c("#26B2E3", "#000000"))(51), colorRampPalette(c( "#000000", "#F4EB08"))(50))) ) %v%
  Heatmap(CBCGALipid_mat, show_row_names = F, show_column_names = F, height = unit(1, 'cm'), cluster_columns = F, cluster_rows = F, na_col = '#eaeaea',
          show_row_dend = F,right_annotation = right_annotation_lipid,
          col = colorRamp2(mybreaks, c(colorRampPalette(c("#26B2E3", "#000000"))(51), colorRampPalette(c( "#000000", "#F4EB08"))(50))) )

graph2pdf(CBCGA_landscape,file = "/results/Fig1.pdf",width = 18, height = 15)

#------------------------------------------------------------------------------------#
# FigS1A
#------------------------------------------------------------------------------------#
clinic_CBCGA = CBCGA_Cohort.Info
clinic_CBCGA$Clinical_data = "Yes"
clinic_CBCGA$Clinical_data = factor(clinic_CBCGA$Clinical_data,levels = c("Yes","No"))
clinic_CBCGA$`Exome sequencing (paired)` = factor(clinic_CBCGA$`Exome sequencing (paired)`,levels = c("Yes","No"))
clinic_CBCGA$`OncoScan array` = factor(clinic_CBCGA$`OncoScan array`,levels = c("Yes","No"))
clinic_CBCGA$`RNA sequencing` = factor(clinic_CBCGA$`RNA sequencing`,levels = c("Yes","No"))
clinic_CBCGA$Proteomics = factor(clinic_CBCGA$Proteomics,levels = c("Yes","No"))
clinic_CBCGA$Metabolomics = factor(clinic_CBCGA$Metabolomics,levels = c("Yes","No"))
clinic_CBCGA$Pathomics = factor(clinic_CBCGA$Pathomics,levels = c("Yes","No"))
clinic_CBCGA$Radiomics = factor(clinic_CBCGA$Radiomics,levels = c("Yes","No"))
#clinic_CBCGA$Imageomic = factor(clinic_CBCGA$Imageomic, levels = c("Both","Digital_pathology","Radiomics","No"))

ord = rownames(clinic_CBCGA)[order(clinic_CBCGA$Clinical_data,clinic_CBCGA$`RNA sequencing`, 
                                   clinic_CBCGA$`Exome sequencing (paired)`,clinic_CBCGA$`OncoScan array`,
                                   clinic_CBCGA$Proteomics,
                                   clinic_CBCGA$Metabolomics,
                                   #clinic_CBCGA$Imageomic,
                                   clinic_CBCGA$Pathomics,clinic_CBCGA$Radiomics, 
                                   decreasing = F)]

clinic_CBCGA = clinic_CBCGA[ord,]
##### simulated data --------
sim = as.data.frame( matrix(rep(1, nrow(clinic_CBCGA)), nrow = 1, ncol = nrow(clinic_CBCGA)))
colnames(sim) = ord

###### draw -------
p = Heatmap(sim , 
            #alter_fun = alter_fun,col = col,
            #alter_fun_is_vectorized = F,
            #show_pct = T, pct_side = 'right',
            column_title = '',
            row_names_side = "left", show_column_names = F,# remove_empty_columns = F,
            row_names_gp = gpar(fontsize = 9), 
            column_title_gp = gpar(fontsize = 10),
            #column_split = factor(clinic_CBCGA$Clinical_Subtype, levels = c('HR+HER2-', 'HR+HER2+', 'HR-HER2+','TNBC')),
            column_order = colnames(sim), #row_order = cnvlist_ord,
            gap = unit(0.5,"inch"),
            top_annotation = HeatmapAnnotation(Clinical_data = clinic_CBCGA$Clinical_data,
                                               RNA = clinic_CBCGA$`RNA sequencing`,
                                               Mutation = clinic_CBCGA$`Exome sequencing (paired)`,
                                               CNA = clinic_CBCGA$`OncoScan array`, 
                                               Proteomics = clinic_CBCGA$Proteomics,
                                               Metabolomics = clinic_CBCGA$Metabolomics,
                                               #Imageomic = clinic_CBCGA$Imageomic,
                                               Digital_pathology = clinic_CBCGA$Pathomics,
                                               Radiomics = clinic_CBCGA$Radiomics,
                                               annotation_name_side = "left",
                                               annotation_name_gp = gpar(fontsize = 9),
                                               col = list(Clinical_data = c('Yes' = '#93C6C1', 'No' = 'white'), 
                                                          Mutation = c('Yes' = '#4F85C1', 'No' = 'white'), 
                                                          CNA = c('Yes' = '#91B720', 'No' = 'white'), 
                                                          RNA = c('Yes' = '#DFAF1E', 'No' = 'white'), 
                                                          Proteomics = c('Yes' = '#BD5135', 'No' = 'white'), 
                                                          Metabolomics = c('Yes' = '#AE1426', 'No' = 'white'), 
                                                          # Imageomic = c('Radiomics' = '#6682AA','Digital_pathology' = '#9B2A59',
                                                          #              'Both' = '#5F2783','No' = 'white')
                                                          Radiomics = c('Yes' = '#9B2A59', 'No' = 'white'),
                                                          Digital_pathology = c('Yes' = '#5F2783', 'No' = 'white')
                                               )
            )
) 
graph2pdf(p, file = "/results/FigS1A.pdf",height = 10,width = 11.15)


#------------------------------------------------------------------------------------#
# FigS1B
#------------------------------------------------------------------------------------#
clinic_CBCGA = CBCGA_Cohort.Info
clinic_CBCGA = clinic_CBCGA[!is.na(clinic_CBCGA$`Intrinsic subtype (PAM50)`),]
cluster.merge4plot = as.data.frame(table(clinic_CBCGA$`Intrinsic subtype (AIMS)`, clinic_CBCGA$`Intrinsic subtype (PAM50)`))
colnames(cluster.merge4plot)[1:2] = c("AIMS","PAM50")

####
cluster.merge4plot2 = ddply(cluster.merge4plot,"AIMS", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$PAM50 = as.character(cluster.merge4plot2$PAM50)
cluster.merge4plot2$AIMS = as.character(cluster.merge4plot2$AIMS)

cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Her2"] = "HER2-enriched"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumB"] = "Luminal B"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Basal"] = "Basal-like"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Normal"] = "Normal-like"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumA"] = "Luminal A"

cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "Her2"] = "HER2-enriched"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "LumB"] = "Luminal B"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "Basal"] = "Basal-like"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "Normal"] = "Normal-like"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "LumA"] = "Luminal A"


cluster.merge4plot2$PAM50 = factor(cluster.merge4plot2$PAM50, levels = c("Normal-like","Luminal B","Luminal A","HER2-enriched","Basal-like"))
cluster.merge4plot2$AIMS = factor(cluster.merge4plot2$AIMS, levels = rev(c("Normal-like","Luminal B","Luminal A","HER2-enriched","Basal-like")))

color = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
names(color) = c("Luminal A","Luminal B","HER2-enriched","Basal-like","Normal-like")

tmp = fisher.test(table(clinic_CBCGA$`Intrinsic subtype (AIMS)`, clinic_CBCGA$`Intrinsic subtype (PAM50)`),simulate.p.value=TRUE)
p = ggplot(cluster.merge4plot2, aes(x=AIMS, y = percent, fill = PAM50 )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],3)) ) 
graph2pdf(p, file = "./results/FigS1B_part1.pdf")

######
cluster.merge4plot2 = ddply(cluster.merge4plot,"PAM50", transform,
                            percent = Freq / sum(Freq) *100)
cluster.merge4plot2$PAM50 = as.character(cluster.merge4plot2$PAM50)
cluster.merge4plot2$AIMS = as.character(cluster.merge4plot2$AIMS)

cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Her2"] = "HER2-enriched"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumB"] = "Luminal B"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Basal"] = "Basal-like"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "Normal"] = "Normal-like"
cluster.merge4plot2$PAM50[cluster.merge4plot2$PAM50 == "LumA"] = "Luminal A"

cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "Her2"] = "HER2-enriched"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "LumB"] = "Luminal B"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "Basal"] = "Basal-like"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "Normal"] = "Normal-like"
cluster.merge4plot2$AIMS[cluster.merge4plot2$AIMS == "LumA"] = "Luminal A"


cluster.merge4plot2$PAM50 = factor(cluster.merge4plot2$PAM50, levels = rev(c("Normal-like","Luminal B","Luminal A","HER2-enriched","Basal-like")) )
cluster.merge4plot2$AIMS = factor(cluster.merge4plot2$AIMS, levels = c("Normal-like","Luminal B","Luminal A","HER2-enriched","Basal-like")) 

color = c("#1D76BC","#76CFE6","#6E59A6","#E11D2E","#CDCFD0")
names(color) = c("Luminal A","Luminal B","HER2-enriched","Basal-like","Normal-like")

tmp = fisher.test(table(clinic_CBCGA$`Intrinsic subtype (AIMS)`, clinic_CBCGA$`Intrinsic subtype (PAM50)`),simulate.p.value=TRUE)
p = ggplot(cluster.merge4plot2, aes(x=PAM50, y = percent, fill = AIMS )) + geom_bar(stat = "identity") +
  scale_fill_manual(values = color)+
  theme(panel.border = element_rect(fill=NA, colour = "black", size=0.5), 
        axis.ticks = element_line(size = 0.5),
        axis.text.x=element_text(size = 10),
        axis.text.y=element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank() 
  ) +
  ggtitle(paste0("p = ",round(tmp[["p.value"]],3)) ) 
graph2pdf(p, file = "./results/FigS1B_part2.pdf")



#------------------------------------------------------------------------------------#
# FigS1D-E
#------------------------------------------------------------------------------------#
clinic_CBCGA = CBCGA_Cohort.Info
colnames(clinic_CBCGA)[colnames(clinic_CBCGA) == "Intrinsic subtype (PAM50)"] = "PAM50"
colnames(clinic_CBCGA)[colnames(clinic_CBCGA) == "Clinical subtype"] = "Clinical_subtype"

colnames(clinic_CBCGA)[colnames(clinic_CBCGA) == "RFS time (month)"] = "RFS_months"
colnames(clinic_CBCGA)[colnames(clinic_CBCGA) == "RFS status"] = "RFS_status"

clinic_CBCGA$PAM50 = factor(clinic_CBCGA$PAM50, levels = c("LumA", "LumB","Her2","Basal","Normal"))
clinic_CBCGA$Clinical_subtype = factor(clinic_CBCGA$Clinical_subtype, levels = c("HR+HER2-", "HR+HER2+","HR-HER2+","TNBC"))

km = survfit(Surv(RFS_months,RFS_status) ~ PAM50,data = clinic_CBCGA)
#p = ggsurvplot(km,pval = T,pval.method = T,palette =  c("#1B77BE",  "#7ACBE2",  "#6354A2", "#E4052D","#CECECF"),
#           break.x.by = 12,break.y.by = 0.1, ylim = c(0.4,1),xlim = c(0,108) ) 
#ggsave(p$plot,filename = "FigS1D.pdf",height = 6,width = 6) 

km = survfit(Surv(RFS_months,RFS_status) ~ Clinical_subtype,data = clinic_CBCGA)
#p = ggsurvplot(km,pval = T,pval.method = T,palette =  c("#2085C4",  "#78B829",  "#F1AE15", "#DC4F35"),
#           break.x.by = 12,break.y.by = 0.1, ylim = c(0.4,1),xlim = c(0,108) ) 
#ggsave(p$plot,filename = "/results/FigS1E.pdf",height = 6,width = 6) 


