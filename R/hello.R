# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

hello <- function() {
  print("Hello, world!")
}

PCA_GD <- function(data_path, metadata_path,output_filename = "geno_data",sep = ",",pruned_data = "no",maf = 0.01,geno = 0.02,indep_pairwise = c(50,5,0.2),hwe = 1e-8,max_alleles = 2){

  source("/data/sata_data1/ab1/devashish/PhD_research/Rscripts/ggplot_theme_Publication-2.R")
  df <- read.csv(metadata_path,sep = ",")

  if (pruned_data == "no"){
  # performed maf,geno, LD pruning, 50 SNP window, sliding window of 5 and r2 greater than 0.2
  #data_path <- "/data/sata_data1/ab1/devashish/PhD_research/datasets/GTEx_final_genotype_gene_expression/Genotype_data/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz"
    plink_cmd <- paste0("plink2 --vcf ",data_path," --chr 1-22 --maf ",maf," --geno ",geno,
                        " --indep-pairwise ",indep_pairwise[1]," ",indep_pairwise[2]," ",indep_pairwise[3],
                        "  --hwe ",hwe," --max-alleles ",max_alleles," --rm-dup exclude-mismatch --snps-only --make-bed --out ",output_filename,"_QC")

    system(plink_cmd)

    plink_cmd_fin <- paste0("plink2 --bfile ",output_filename,"_QC --extract ",output_filename,"_QC.prune.in --make-bed --out ",output_filename,"_QC_PCA_PRUNED")

    system(plink_cmd_fin)

    indiv_data  = read.table(paste0(output_filename,"_QC_PCA_PRUNED.fam"),header = F)

    plink_pca <- paste0("plink2 --bfile ",output_filename,"_QC_PCA_PRUNED --pca ",nrow(indiv_data)," --out ",output_filename,"_QC_PCA_RESULT")
    system(plink_pca)

  }else{
    indiv_data  = read.table(data_path,header = F)
    plink_pca <- paste0("plink2 --bfile ",data_path,"_QC_PCA_PRUNED --pca ",nrow(indiv_data)," --out ",output_filename,"_QC_PCA_RESULT")
    system(plink_pca)

  }
  #pca_eig_vec_2 <- merge(pca_eig_vec,df,by.x = "sample.id",by.y = "SUBJID")

  #write.table(pca_eig_vec_2,"pcadata_GTEx_V9.txt",col.names = T,row.names = F,quote = F,sep = "\t")

}

PCA_plot <- function(pcadata_path,metadata_path,sep = ",",x = 1,y = 2,color_var = "population",shape_var = NULL,shape_vec = NULL,figfile = "PCA_plot"){

  library(dplyr)
  library(ggplot2)
  library(randomcoloR)

  df <- read.csv(metadata_path,sep = sep)
  pca_eig_vec <- read.table(paste0(pcadata_path,".eigenvec"))
  pca_eig_vec <- pca_eig_vec[,-1]
  colnames(pca_eig_vec) <- c("sample.id",paste0("EV",1:nrow(df)))
  pca_eig_vec_2 <- merge(pca_eig_vec,df,by = "sample.id",all.x = T)

  EV_x <- ifelse(x > 0, pca_eig_vec[, paste0("EV",x)], -pca_eig_vec[, paste0("EV",abs(x))])
  EV_y <- ifelse(y > 0, pca_eig_vec[, paste0("EV",y)], -pca_eig_vec[, paste0("EV",abs(y))])

  pca_eig_val <- scan(paste0(pcadata,".eigenval"))
  pca_eig_val_pt <- (pca_eig_val/sum(pca_eig_val))*100
  pca_eig_val_pt <- round(pca_eig_val_pt,2)

  palette <- distinctColorPalette(length(unique(pca_eig_vec_2[,shape_var])))

  if (!is.null(shape_var)){
    pca_eig_vec_2[,shape_var] = factor(pca_eig_vec_2[,shape_var])
    shape_df = pca_eig_vec_2[,shape_var]
  }

  shape_var <- shape_vec[levels(pca_eig_vec_2[,shape_var])]

  png(paste0(figfile,".png"),width=4500, height=3500, res = 300)

  p <- ggplot(pca_eig_vec_2,aes(x=EV_x ,y=EV_y,fill=pca_eig_vec_3[,color_var],shape = shape_df)) +
    geom_point(size = 2) +
    xlab(paste("PC",x,"(", pca_eig_val_pt[x],"%)",sep=" ")) +
    ylab(paste("PC",y,"(", pca_eig_val_pt[y],"%)",sep=" ")) +
    theme(legend.title=element_text(size=20), legend.text=element_text(size=18),axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14)) +
    theme_bw() +
    scale_fill_manual(values = palette,guide = guide_legend(override.aes = list(shape = 21)))+
    scale_shape_manual(values = shape_vec)+
    labs(fill=color_var,shape=shape_var)

  p <- p + theme_Publication()

  print(p)

  dev.off()

  library(dplyr)

  PC_summary_df <-  pca_eig_vec_2 %>%
    group_by_at(shape_var,color_var) %>%
    summarise_at(.vars = c(paste0("EV",x),paste0("EV",y)), .funs = mean)

  PC_summary_df <- data.frame(PC_summary_df)

  EV_x1 <- ifelse(x > 0, PC_summary_df[, paste0("EV",x)], -PC_summary_df[, paste0("EV",abs(x))])
  EV_y1 <- ifelse(y > 0, PC_summary_df[, paste0("EV",y)], -PC_summary_df[, paste0("EV",abs(y))])

  png(paste0(fig_file,"_summary.png"),width=4500, height=3500, res = 300)

  p <- ggplot(PC_summary_df,aes(x=EV_x1 ,y=EV_y1,fill=pca_eig_vec_3[,color_var],shape = shape_df)) +
    geom_point(size = 5) +
    xlab(paste("PC",x,"(", pca_eig_val_pt[x],"%)",sep=" ")) +
    ylab(paste("PC",y,"(", pca_eig_val_pt[y],"%)",sep=" ")) +
    theme(legend.title=element_text(size=20), legend.text=element_text(size=18),axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14)) +
    theme_bw() +
    scale_fill_manual(values = palette,guide = guide_legend(override.aes = list(shape = 21)))+
    scale_shape_manual(values = shape_vec)+
    labs(fill="linguistic group",shape="ethnic group")
  p <- p + theme_Publication()

  print(p)

  dev.off()

}

PCA_outliers <- function(metadata_path,pcadata_path,color_var,x = 1,y=2,shape_var = NULL, shape_vec = NULL){

  df <- read.csv(metadata_path,sep = sep)
  pca_eig_vec <- read.table(paste0(pcadata_path,".eigenvec"))
  pca_eig_vec <- pca_eig_vec[,-1]
  colnames(pca_eig_vec) <- c("sample.id",paste0("EV",1:nrow(df)))
  EV_x <- ifelse(x > 0, pca_eig_vec[, paste0("EV",x)], -pca_eig_vec[, paste0("EV",abs(x))])
  EV_y <- ifelse(y > 0, pca_eig_vec[, paste0("EV",y)], -pca_eig_vec[, paste0("EV",abs(y))])

  pca_eig_vec_2 <- merge(pca_eig_vec,df,by = "sample.id",all.x = T)

  pca_eig_val <- scan(paste0(pcadata,".eigenval"))
  pca_eig_val_pt <- (pca_eig_val/sum(pca_eig_val))*100
  pca_eig_val_pt <- round(pca_eig_val_pt,2)

  findf <- data.frame()

  for (pop in unique(df[,color_var])){
    df2 <- df[df[,color_var] == pop,]
    df2$D2 <- mahalanobis(df2[,paste0("EV",1:3)],center = colMeans(df2[,paste0("EV",1:3)]),cov=cov(df2[,paste0("EV",1:3)]))
    df2 <- df2[order(df2$D2,decreasing=T),]
    df2$pvalue <- pchisq(df2$D2,lower.tail = F, df=3)
    findf <- rbind(findf,df2)
  }

  findf[is.na(findf$D2),c("D2","pvalue")] <- 1
  write.table(findf,"PCA_data_outlier_pvalue.txt",col.names=T,row.names=F,sep = "\t",quote =F)

  dfsub <- findf[findf$pvalue <= 1e-6,]
  write.table(dfsub[,c(2,2)],"outlier_indiv_info.txt",sep="\t",col.names=F,row.names=F,quote = F)


  pdf("pop_specific_plots_D2_pvalue.pdf", width = 12, height = 12)

  for (pop in unique(findf[,color_var])){
    df2 <- findf[findf[,color_var] == pop,]
    p <- ggplot(findf,aes(x=EV1 ,y=EV2,shape = shape_var)) +
      geom_point(alpha = 0.05,fill = "blue",size = 2) +
      geom_point(data = df2, aes(x = PC1, y = PC2, shape = shape_var,size = -log10(pvalue), fill = ifelse(-log10(pvalue) > 6, "red", "black")))+
      xlab(paste("PC1(", pca_eig_val_pt[1],"%)",sep=" ")) +
      ylab(paste("PC2(", pca_eig_val_pt[2],"%)",sep=" ")) +
      theme(legend.title=element_text(size=20), legend.text=element_text(size=18),axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14),legend.position = NA) +
      theme_bw() +
      scale_fill_manual(values = c("black", "red"),guide = "none")+
      scale_shape_manual(values = shape_vec)+
      labs(shape=shape_var,title = paste0("Population: ",pop,",",shape_var,":",paste(unique(df2$Institute),collapse=",")))


    p <- p + theme_Publication()

    print(p)

  }
  dev.off()

}
