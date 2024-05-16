
PCA_GD <- function(data_path, metadata_path,geno_format = "pfile",
                   output_filename = "geno_data",pruned_data = "no",
                   maf = 0.01,geno = 0.02,indep_pairwise = c(50,5,0.2),hwe = 1e-8,
                   max_alleles = 2){

  df <- read.csv(metadata_path,check.names = F,sep = "\t")

  if (pruned_data == "no"){
  # performed maf,geno, LD pruning, 50 SNP window, sliding window of 5 and r2 greater than 0.2
    plink_cmd <- paste0("plink2 --",geno_format," ",data_path," --chr 1-22 --maf ",maf," --geno ",geno,
                        " --indep-pairwise ",indep_pairwise[1]," ",indep_pairwise[2]," ",indep_pairwise[3],
                        "  --hwe ",hwe," --max-alleles ",max_alleles," --rm-dup exclude-mismatch --snps-only --make-bed --out ",output_filename,"_QC")

    system(plink_cmd)

    plink_cmd_fin <- paste0("plink2 --bfile ",output_filename,"_QC --extract ",output_filename,"_QC.prune.in --make-bed --out ",output_filename,"_QC_PCA_PRUNED")

    system(plink_cmd_fin)

    indiv_data  = read.table(paste0(output_filename,"_QC_PCA_PRUNED.fam"),header = F)

    plink_pca <- paste0("plink2 --bfile ",output_filename,"_QC_PCA_PRUNED --pca ",nrow(indiv_data)," --out ",output_filename,"_QC_PCA_RESULT")
    system(plink_pca)

  }else{
    indiv_data  = read.table(metadata_path,header = T)
    plink_pca <- paste0("plink2 --bfile ",data_path," --pca ",nrow(indiv_data)," --out ",output_filename,"_QC_PCA_RESULT")
    system(plink_pca)

  }

}

PCA_plot <- function(pcadata_path,metadata_path,
                     sampleid_var,color_var,x = 1,y = 2,
                     shape_var = NULL,color_vec = NULL,shape_vec= NULL,
                     figfile = "PCA_plot",...){

  theme_Publication <- function(base_size=14, base_family="sans") {
    suppressWarnings(library(grid))
    suppressWarnings(library(ggthemes))
    (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour="black"),
              axis.ticks = element_line(),
              panel.grid.major = element_line(colour="#f0f0f0"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.box = "vetical",
              legend.key.size= unit(0.5, "cm"),
              #legend.margin = unit(0, "cm"),
              legend.title = element_text(face="italic"),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")
      ))

  }
  col_sel <- function(x,df){
    if(x > 0){
      return(df[,paste0("EV",x)])
    } else {
      return(-df[,paste0("EV",abs(x))])
    }
  }

  suppressPackageStartupMessages(library(dplyr))
  suppressWarnings(library(ggplot2))
  suppressWarnings(library(randomcoloR))

  df <- read.csv(metadata_path,check.names = F,sep = "\t")
  pca_eig_vec <- read.table(paste0(pcadata_path,".eigenvec"))
  pca_eig_vec <- pca_eig_vec[,-1]
  colnames(pca_eig_vec) <- c(sampleid_var,paste0("EV",1:(ncol(pca_eig_vec)-1)))
  pca_eig_vec_2 <- merge(pca_eig_vec,df,by = sampleid_var)

  EV_x <- col_sel(x,pca_eig_vec_2)
  EV_y <- col_sel(y,pca_eig_vec_2)

  pca_eig_val <- scan(paste0(pcadata_path,".eigenval"))
  pca_eig_val_pt <- (pca_eig_val/sum(pca_eig_val))*100
  pca_eig_val_pt <- round(pca_eig_val_pt,2)

  if (is.null(color_vec)){
    palette <- distinctColorPalette(length(unique(pca_eig_vec_2[,color_var])))
  }else{
    palette <- color_vec
  }
  if (!is.null(shape_var)){
    pca_eig_vec_2[,shape_var] = factor(pca_eig_vec_2[,shape_var])
    shape_df = pca_eig_vec_2[,shape_var]

    if(is.null(shape_vec)){
      shape_vec <- c(21:25,sample(7:14,(length(unique(pca_eig_vec_2[,shape_var]))-5)))
      names(shape_vec) <- levels(pca_eig_vec_2[,shape_var])
    }else{
      shape_vec <- shape_vec[levels(pca_eig_vec_2[,shape_var])]
      names(shape_vec) <- levels(pca_eig_vec_2[,shape_var])
    }
    col_shape_df <- pca_eig_vec_2[,c(shape_var,color_var)]
    col_shape_df <- col_shape_df[!duplicated(col_shape_df[,color_var]),]
    rownames(col_shape_df) <- NULL
    col_shape_df$shape_vec <- shape_vec[as.character(col_shape_df[,shape_var])]
  }

  if (!is.null(shape_var) & !is.null(color_var)){
    png(paste0(figfile,".png"),width=4500, height=3500, res = 300)

    p <- ggplot(pca_eig_vec_2,aes(x=EV_x ,y=EV_y,fill=pca_eig_vec_2[,color_var],shape = shape_df)) +
    geom_point(...) +
    xlab(paste("PC",x,"(", pca_eig_val_pt[x],"%)",sep=" ")) +
    ylab(paste("PC",y,"(", pca_eig_val_pt[y],"%)",sep=" ")) +
    theme(legend.title=element_text(size=20), legend.text=element_text(size=18),axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14)) +
    theme_bw() +
    scale_shape_manual(values = shape_vec)+
    scale_fill_manual(values = palette,guide = guide_legend(override.aes = list(shape = col_shape_df$shape_vec)))+
    labs(fill=color_var,shape=shape_var,title =paste0(figfile))

    p <- p + theme_Publication()
    print(p)
    dev.off()

  } else if(is.null(shape_var) & !is.null(color_var)){
    png(paste0(figfile,".png"),width=4500, height=3500, res = 300)

    p <- ggplot(pca_eig_vec_2,aes(x=EV_x ,y=EV_y,fill=pca_eig_vec_2[,color_var])) +
      geom_point(shape = 21,...) +
      xlab(paste("PC",x,"(", pca_eig_val_pt[x],"%)",sep=" ")) +
      ylab(paste("PC",y,"(", pca_eig_val_pt[y],"%)",sep=" ")) +
      theme(legend.title=element_text(size=20), legend.text=element_text(size=18),
            axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14)) +
      theme_bw() +
      scale_fill_manual(values = palette,guide = guide_legend(override.aes = list(shape = 21)))+
      labs(fill=color_var,title =paste0(figfile))

    p <- p + theme_Publication()
    print(p)
    dev.off()
  }

  if (!is.null(shape_var) & !is.null(color_var)){
  PC_summary_df <-  data.frame(pca_eig_vec_2 %>%
                                group_by(!!sym(shape_var), !!sym(color_var)) %>%
                                 summarise(across(.cols = c(paste0("EV", abs(x)), paste0("EV", abs(y))), .fns = mean)))
  } else if(is.null(shape_var) & !is.null(color_var)){
    PC_summary_df <-  data.frame(pca_eig_vec_2 %>%
                                   group_by(!!sym(color_var)) %>%
                                   summarise(across(.cols = c(paste0("EV", abs(x)), paste0("EV", abs(y))), .fns = mean)))
  }
  PC_summary_df <- data.frame(PC_summary_df)

  EV_x1 <- col_sel(x,PC_summary_df)
  EV_y1 <- col_sel(y,PC_summary_df)

  if (!is.null(shape_var)){
    PC_summary_df[,shape_var] = factor(PC_summary_df[,shape_var])
    shape_df1 = PC_summary_df[,shape_var]
  }

  if (!is.null(shape_var) & !is.null(color_var)){

  png(paste0(figfile,"_summary.png"),width=4500, height=3500, res = 300)

  p1 <- ggplot(PC_summary_df,aes(x=EV_x1 ,y=EV_y1,fill=PC_summary_df[,color_var],shape = shape_df1)) +
    geom_point(...) +
    xlab(paste("PC",x,"(", pca_eig_val_pt[x],"%)",sep=" ")) +
    ylab(paste("PC",y,"(", pca_eig_val_pt[y],"%)",sep=" ")) +
    theme(legend.title=element_text(size=20), legend.text=element_text(size=18),
          axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14)) +
    theme_bw() +
    scale_fill_manual(values = palette,guide = guide_legend(override.aes = list(shape = col_shape_df$shape_vec)))+
    scale_shape_manual(values = shape_vec)+
    labs(fill=color_var,shape=shape_var,title =paste0(figfile,"_summary"))
  p1 <- p1 + theme_Publication()

  print(p1)
  dev.off()

  } else if(is.null(shape_var) & !is.null(color_var)){
    png(paste0(figfile,"_summary.png"),width=4500, height=3500, res = 300)

    p1 <- ggplot(PC_summary_df,aes(x=EV_x1 ,y=EV_y1,fill=PC_summary_df[,color_var])) +
      geom_point(shape = 21,...) +
      xlab(paste("PC",x,"(", pca_eig_val_pt[x],"%)",sep=" ")) +
      ylab(paste("PC",y,"(", pca_eig_val_pt[y],"%)",sep=" ")) +
      theme(legend.title=element_text(size=20), legend.text=element_text(size=18),
            axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14)) +
      theme_bw() +
      scale_fill_manual(values = palette,guide = guide_legend(override.aes = list(shape = 21)))+
      labs(fill=color_var,title =paste0(figfile,"_summary"))

    p1 <- p1 + theme_Publication()
    print(p1)
    dev.off()
  }
  return(list(p,p1))
}

PCA_outliers <- function(pcadata_path,metadata_path,color_var,
                         sampleid_var,p_threshold = 1e-4,
                         ncomp = 3,x = 1,y=2,shape_var = NULL,...){
  theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
      + theme(plot.title = element_text(face = "bold",
                                        size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
              text = element_text(),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.border = element_rect(colour = NA),
              axis.title = element_text(face = "bold",size = rel(1)),
              axis.title.y = element_text(angle=90,vjust =2),
              axis.title.x = element_text(vjust = -0.2),
              axis.text = element_text(),
              axis.line.x = element_line(colour="black"),
              axis.line.y = element_line(colour="black"),
              axis.ticks = element_line(),
              panel.grid.major = element_line(colour="#f0f0f0"),
              panel.grid.minor = element_blank(),
              legend.key = element_rect(colour = NA),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.box = "vetical",
              legend.key.size= unit(0.5, "cm"),
              #legend.margin = unit(0, "cm"),
              legend.title = element_text(face="italic"),
              plot.margin=unit(c(10,5,5,5),"mm"),
              strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
              strip.text = element_text(face="bold")
      ))

  }

    col_sel <- function(x,df){
    if(x > 0){
      return(df[,paste0("EV",x)])
    } else {
      return(-df[,paste0("EV",abs(x))])
    }
  }

    df <- read.csv(metadata_path,check.names = F,sep = "\t")
  pca_eig_vec <- read.table(paste0(pcadata_path,".eigenvec"))
  pca_eig_vec <- pca_eig_vec[,-1]
  colnames(pca_eig_vec) <- c(sampleid_var,paste0("EV",1:nrow(df)))

  pca_eig_vec_2 <- merge(pca_eig_vec,df,by = sampleid_var)
  EV_x <- col_sel(x,pca_eig_vec_2)
  EV_y <- col_sel(y,pca_eig_vec_2)

  pca_eig_val <- scan(paste0(pcadata_path,".eigenval"))
  pca_eig_val_pt <- (pca_eig_val/sum(pca_eig_val))*100
  pca_eig_val_pt <- round(pca_eig_val_pt,2)

  findf <- data.frame()

  for (pop in unique(pca_eig_vec_2[,color_var])){
    df2 <- pca_eig_vec_2[pca_eig_vec_2[,color_var] == pop,]
    df2$D2 <- mahalanobis(df2[,paste0("EV",1:ncomp)],
                          center = colMeans(df2[,paste0("EV",1:ncomp)]),cov=cov(df2[,paste0("EV",1:3)]))
    df2 <- df2[order(df2$D2,decreasing=T),]
    df2$pvalue <- pchisq(df2$D2,lower.tail = F, df=ncomp)
    findf <- rbind(findf,df2)
  }

  findf[is.na(findf$D2),c("D2","pvalue")] <- 1
  findf <- findf[,c(sampleid_var,color_var,paste0("EV",1:ncomp),"D2","pvalue")]
  write.table(findf,"PCA_data_outlier_pvalue.txt",col.names=T,row.names=F,sep = "\t",quote =F)

  print(head(findf))

  pdf("pop_specific_plots_D2_pvalue.pdf", width = 12, height = 12)

  for (pop in unique(findf[,color_var])){
    df2 <- findf[findf[,color_var] == pop,]
    p <- ggplot(findf,aes(x=EV_x ,y=EV_y)) +
      geom_point(alpha = 0.05,fill = "blue",shape = 21,...) +
      geom_point(data = df2, shape = 21,aes(x = EV1, y = EV2,size = -log10(pvalue),
                                 fill = ifelse(-log10(pvalue) > -log10(p_threshold), "red", "black")),...)+
      xlab(paste("PC1(", pca_eig_val_pt[1],"%)",sep=" ")) +
      ylab(paste("PC2(", pca_eig_val_pt[2],"%)",sep=" ")) +
      theme(legend.title=element_text(size=20),
            legend.text=element_text(size=18),axis.text.x = element_text(size =14),axis.text.y = element_text(size = 14),legend.position = NA) +
      theme_bw() +
      scale_fill_manual(values = c("black", "red"),guide = "none")

    if (is.null(shape_var)){
     p <- p+ labs(shape=shape_var,title = paste0(color_var,": ",pop))
    }else{
      p <- p+ labs(shape=shape_var,
                  title = paste0(color_var,": ",pop,",",shape_var,":",paste(unique(df2[,shape_var]),collapse=",")))
    }

    p <- p + theme_Publication()
    print(p)
   }
  dev.off()

}
