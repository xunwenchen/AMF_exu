plot_box <- function(df, x, y, z){
  ggplot2::theme_set(theme_bw())
  
  plot <- ggplot2::ggplot(df, aes(x = {{x}}, y = {{y}}), fill = {{z}})+
    geom_boxplot()+
    # geom_jitter(alpha = 0.5, width = 0.2)+
    stat_summary(fun=mean, geom="point", 
                 shape=20, size=5, color="indianred")
  
  return(plot)
}


plot_w_unifrac <- function(phylo_ob, by){
  library(phyloseq)
  ordu <-  ordinate(phylo_ob, "PCoA", "unifrac", weighted = TRUE)
  plot_ordination(phylo_ob, ordu, color=by, shape=by)+
    geom_point(size=5, alpha=0.5)+
    scale_colour_brewer(type="qual", palette="Set1")+
    ggtitle("PCoA on weighted UniFrac distance")
  
}


plot_uw_unifrac <- function(phylo_ob, by){
  library(phyloseq)
  ordu_uw = ordinate(phylo_ob, "PCoA", "unifrac", weighted = FALSE)
  plot_ordination(phylo_ob, ordu_uw, color=by, shape=by)+
  geom_point(size=5, alpha=0.5)+
  scale_colour_brewer(type="qual", palette="Set1")+
  ggtitle("PCoA on unweighted UniFrac distance")
}


# # The following not working
# adonis <- function(phylo_ob, by){
#   dist <- phyloseq::distance(phylo_ob, method ='wunifrac')
#   library(pairwiseAdonis)
#   result <- pairwise.adonis(dist, sample_data(phylo_ob)$by, perm = 999)
#   
#   
# }


# # The following not working
# plot_uw_unifrac2 <- function(phylo_ob, var1, x, var2, y, by){
#   library(phyloseq)
#   subset1 <- subset_samples(phylo_ob, {{var1}} == x)
#   subset2 <- subset_samples(subset1, {{var2}} == y)
#   
#   ordu_uw = ordinate(subset2, "PCoA", "unifrac", weighted = FALSE)
#   plot_ordination(subset2, ordu_uw, color=by, shape=by)+
#     geom_point(size=5, alpha=0.5)+
#     scale_colour_brewer(type="qual", palette="Set1")+
#     ggtitle("PCoA on unweighted UniFrac distance")
# }


## Create function
## function to plot within-group beta diversity distance
beta_boxplot <- function(physeq, method = "bray", group) {
  
  # physeq: phyloseq-class object
  # method: beta-diversity metric. Default "bray", i.e., Bray-Curtis dissimilarity 
  # group: factorial variable to group
  
  ## Packages
  require("phyloseq") # v.1.30.0
  require("ggplot2") # v.3.3.2
  
  ## Identify the correspondence: group and samples
  group2samp <- list() # list to save the correspondence between group <--> samples
  group_list <- as.factor(get_variable(sample_data(physeq), group)) # list of group elements
  for (groups in levels(group_list)) { # loop over the no. of group levels
    target_group <- which(group_list == groups) # vct pos of the curr group variable 
    group2samp[[ groups ]] <- sample_names(physeq)[target_group] # matching samples: based on vct pos
  }  
  
  ## Calculate beta-diversity
  beta_div_dist <- phyloseq::distance(physeq = physeq, method = method)
  beta_div_dist <- as(beta_div_dist, "matrix")
  
  
  
  
  
  ## Coerce distance mtx into a tidy data frame by group
  dist_df <- data.frame() # save results in df 
  counter <- 1 
  for (groups in names(group2samp)) { # loop over group fct levels 
    sub_dist <- beta_div_dist[ group2samp[[groups]], group2samp[[groups]] ] # subset dist mtx to curr samples
    #print(sub_dist)
    no_samp_col <- ncol(sub_dist) # n cols: curr sub dist
    no_samp_row <- nrow(sub_dist) # n rows: curr sub dist
    for ( cols in seq(no_samp_col) ) { # loop over cols: curr sub_dist
      if ( cols > 1 ) {
        for ( rows in seq((cols-1)) ) { # loop over rows: curr sub_dist 
          ## Save results
          dist_df[ counter, "sample_pair" ] <- paste0( colnames(sub_dist)[cols], "-",  
                                                       rownames(sub_dist)[rows] ) # sample pair
          dist_df[ counter, "group" ] <- groups # group  
          dist_df[ counter, "beta_div_method" ] <- method # method
          dist_df[ counter, "beta_div_value" ] <- sub_dist[rows, cols] # beta-diversity for the sample pair     
          counter = counter + 1
        }
      }
    }
  }
  
  ## Create a ggplot2 boxplot
  plot_boxplot <- ggplot(data = dist_df, aes(x = group, y = beta_div_value, color = group)) + 
    geom_boxplot() + #geom_boxplot(outlier.shape=NA)
    geom_jitter(alpha = 0.5, width = 0.2) + 
    theme_bw() + 
    xlab(group) + ylab(method) 
    # theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  
  ## Save df and boxplot into a list 
  # list_Out <- list("data" = dist_df, "plot" = plot_boxplot) 
  # 
  # return(list_Out)
}


# function to plot between-group beta diversity comparison

beta_boxplot_btw <- function(phylo_ob, dist_method, pair1, pair2){
  
  # calculate distance matrix
  dist <- phyloseq::distance(phylo_ob, method = dist_method)
  dist <- as.matrix(dist)
  
  # assign '0' to upper triangle of the dist matrix
  dist[upper.tri(dist)] <- 0
  
  # using melt() in 'reshape2' package to convert dist matrix to dataframe in common form and delete previously assigned '0'.
  library(reshape2)
  dist_melt <- subset(melt(dist), value!=0)
  
  # assign Sample_ID from metadata file to dist dataframe - dist_melt
  # Var1 and Var2 are column names of dist_melt, Sample_ID is a column name in metadata (load metadata first if have not)
  metadata.in.fun <- read_xlsx('data/metadata.xlsx')
  library(dplyr)
  temp <- dplyr::left_join(dist_melt, metadata.in.fun, by = c ('Var1'='Sample_ID'))
  temp2 <- dplyr::left_join(temp, metadata.in.fun, by = c ('Var2'='Sample_ID'))
  
  # generate new columns, and exclude unwanted columns
  temp3 <- temp2 %>% dplyr::transmute(distance = value, Compartment = Compartment.x, Pair = paste0(Group.x, sep = "_VS_", Group.y))
  
  # see what pairs can be used for comparison
  groups_to_select <- temp3 %>% group_by(Pair) %>% summarise()
  print('select two pairs for comparison:')
  print(groups_to_select)
  # plot 
  plot <- temp3 %>% 
    filter(Pair == pair1 | Pair == pair2) %>% 
    ggplot(aes(x=Pair, y = distance, color = Pair))+
    geom_boxplot()+
    geom_jitter(alpha = 0.5, width = 0.2)
  return(plot)
  
}


log2fc <- function(phylo_ob, Group, trt_group, ref_group){
  library("DESeq2")
  
  # since an error shows "In DESeqDataSet(se, design = design, ignoreRank): some variables in design formula are characters, converting to factors". The variables in design formula ("~ Group" in our case) are converted to factors using the following as.factor() function
  sample_data(phylo_ob)$Group <- as.factor(sample_data(phylo_ob)$Group)

  
  # The following two lines actually do all the complicated DESeq2 work. The function phyloseq_to_deseq2 converts the phyloseq-format microbiome data (pl_hps in our case) into a DESeqDataSet with dispersions estimated, using the experimental design formula, also shown (the ~Group term in our case). The DESeq function does the rest of the testing, in this case with default testing framework, but we can actually use alternatives. 
  deseq2_temp <- phyloseq_to_deseq2(phylo_ob, ~ Group)
  
  deseq2_temp2  <-  DESeq(deseq2_temp, test="Wald", fitType="parametric")
  
  
  # The following results function call creates a table of the results of the tests. Very fast. The hard work was already stored with the rest of the DESeq2-related data in our latest version of the deseq2_temp2 object (see above). I then order by the adjusted p-value, removing the entries with an NA value (seems we did not remove them). The rest of this example is just formatting the results table with taxonomic information for nice(ish) display in the HTML output.
  
  res <- results(deseq2_temp2, 
                 contrast=c(Group,trt_group,ref_group), 
                 cooksCutoff = FALSE)
  res # use this to check whether it is M+Cd vs M, or M vs M+Cd:
  # Group M.Cd vs M 
  # Wald test p-value: Group M.Cd vs M 
  # DataFrame with 20861 rows and 6 columns
  
  
  
  
  alpha = 0.01
  sigtab = res[which(res$padj < alpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(phylo_ob)[rownames(sigtab), ], "matrix"))
  head(sigtab)
  dim(sigtab)
  
  # Let's look at the OTUs that were significantly different between the two groups (M vs M+Cd in our case? I think so) The following makes a nice ggplot2 summary of the results.
  
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  # Phylum level
  # x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
  
  # class level
  
  
  # Order level
  # x = tapply(sigtab$log2FoldChange, sigtab$Order, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Order = factor(as.character(sigtab$Order), levels=names(x))
  
  # Family level
  x = tapply(sigtab$log2FoldChange, sigtab$Family, function(x) max(x))
  x = sort(x, TRUE)
  sigtab$Family = factor(as.character(sigtab$Family), levels=names(x))
  
  
  # Genus level
  # x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
  # x = sort(x, TRUE)
  # sigtab$Genus = factor(as.character(sigtab$Genus), levels=names(x))
  
  
  ggplot(sigtab, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=3) 
  # + 
  #   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
}


flat_mt <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    Source = rownames(cormat)[row(cormat)[ut]],
    Target = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}



metabo_otu_cor_p_tab <- function(df_otu, df_metabo, which_sample, to_which_sample, cor_p_value, cor_R_value){
  # For each treatment in a compartment, only those OTUs detected in all 8 samples were retained. These more-stringent OTU tables were used to correlate with the metabolomic profiles. 
  # omit rows contain 0 otu number
  # df2 <- df_otu[!rowSums(df_otu==0),]
  df2 <- df_otu
  # omit rows with group-max value smaller than 5
  df2$rowMax <- apply(df2, 1, max, na.rm=TRUE) # find row max value assign to column 'rowMax'. apply(df,1,...) --- '1' means manipulation is performed on rows. '2' is on columns
  
  df3 <- df2[df2$rowMax >= 15, ]     # select only rows with rowMax more than x zotu counts (the larger the more stringent)
  df4 <- select(df3, -rowMax) # drop the rowMax column
  
  # remove rows with >=5 occurrence of 0 zotu count
  df4$count.0 <- apply(df4, 1, function(x) length(which(x==0))) # count the no. of occurrence of 0
  df5 <- df4[df4$count.0 <= 0, ]  # retain only those equal or less than x (the smaller the more stringent in this case)
  df6 <- select(df5, -count.0) # drop the count.0 column
  
  # transpose and convert to dataframe
  df7 <- as.data.frame(t(as.data.frame(df6)))
  
  
  # delete compound name which is not needed later
  metabo_df2 <- df_metabo[,-2]
  # transpose
  metabo_df3 <- t(metabo_df2)
  # take the first row and assign it as the colname
  colnames(metabo_df3) <- metabo_df3[1, ]
  # remove first row which is not needed (repeated)
  metabo_df4 <- metabo_df3[-1, ]
  # select a treatment group in the same compartment (e.g., Rhizosphere and Ctrl). Each compartment and each treatment has a network
  metabo_df5 <- as.data.frame(metabo_df4[c(which_sample:to_which_sample),])
  # check if the df is numeric, if not, need to conver to numeric as follow
  metabo_df6 <- as.data.frame(lapply(metabo_df5, as.numeric))
  
  # construct correlation matrix. This matrix is large. Maybe not needed
  # temp_cor_mt <- cor(df3, metabo_df6, 'pearson', use = "complete.obs")
  
  # construct R and P value tables
  library(Hmisc)
  # cor_p_mt is a very large file
  cor_p_mt <- rcorr(as.matrix(df7), as.matrix(metabo_df6), type = 'spearman')
  
  # flattern the tables using the function flat_mt() (own function)
  
  # using the function
  # flat_cor_p_tab is a large file
  flat_cor_p_tab <- flat_mt(cor_p_mt$r, cor_p_mt$P)
  
  
  library(stringr)
  # select only zotu paired with compound (not zotu vs zotu, or com vs com)
  zotu_com_tab <- flat_cor_p_tab %>% 
    filter(str_detect(Source, 'Zotu') & str_detect(Target, 'Com') | str_detect(Source, 'Com') & str_detect(Target, 'Zotu'))
  
  zotu_com_tab_sig <- zotu_com_tab[abs(zotu_com_tab$cor) >= cor_R_value &
                                     zotu_com_tab$p <= cor_p_value, ]
  
  return(zotu_com_tab_sig)
}


construct_node_file <- function(edge_raw, TAX, metabo_df){
  # edge_raw is the melted/flattern zotu-com correlation with p and r value; TAX is the phylo object of tax info. metabo_df is the pos and neg combined metabo file with all treatments and samples. 
  node_raw <- 
    as.data.frame(
      unique(
        c(edge_raw$Source, edge_raw$Target)))
  
  # assign "Name" as column name for later merging.
  colnames(node_raw)[1]  <- "Name" 
  
  # prepare for zotu_tax_tab file
  zotu_tax_tab <- as.data.frame(TAX)
  zotu_tax_tab <- rownames_to_column(zotu_tax_tab, "Name") # here used 'Name' is for later merging with node file 
  
  
  # prepare for metabo file
  metabo_Name_com_name <- metabo_df %>% select(com_id, com_name) %>% rename(Name = com_id)
  
  # first merge
  
  node_temp1 <- merge(node_raw, zotu_tax_tab, by = 'Name', all.x = TRUE)
  
  # last merge
  
  node_temp2 <- merge(node_temp1, metabo_Name_com_name, by = 'Name', all.x = TRUE)
  
  
  # assign a column 'ID' and order before 'Name'. It is needed for later VLookup in excel to create Source and Target numbers, to be used in gephi. 
  
  node_final <- node_temp2 %>% mutate(ID = 1:n(), .before = 'Name')
  
  
  return(node_final)
  
}





node_merge_otu <- function(node_file_address){
  
  node <- read_csv(node_file_address)
  zotu_tax_tab <- as.data.frame(TAX)
  zotu_tax_tab2 <- rownames_to_column(zotu_tax_tab, "Name") # here used 'Name' is for later merging with node file
  node2 <- merge(node, zotu_tax_tab2, by = 'Name', all.x = TRUE)
  node3 <- node2 %>% mutate(ID = 1:n(), .before = 'Name')
  return(node3)
  
}



plot_points_with_line <- function(df, x, y){
  
  ggplot2::theme_set(theme_bw())
  
  plot <- ggplot2::ggplot(df, aes(x = {{x}}, y = {{y}}), group=1)+
    # add group = 1 here is needed to create line connecting the points
    geom_line(color = 'grey', size = 1.5)+
    geom_point(shape=21, color="black", fill="#69b3a2", size=5)
  
  return(plot)
  
}




plot_pi_zi <- function(df){
  
  ggplot(df, aes(Pi, Zi))+
    geom_point(size = 2, alpha = 0.4)+
    geom_hline(yintercept=2.5, linetype="dashed", color = "darkgrey")+
    geom_vline(xintercept=0.62, linetype="dashed", color = "darkgrey")+
    xlab('Among-module connectivity (Pi)')+
    ylab('Within-module connectivity (Zi)')+
    xlim(0, 0.8)+
    ylim(0, 4)
  # set pi = 0.62, zi = 2.5
}




plot_node_degree <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
    # + coord_flip()
  
}



plot_node_degree_2_trt <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                           all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
    # + coord_flip()
  
}



plot_node_degree0 <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}



plot_node_degree_2_trt0 <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Accum_degree = sum(node.degree), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Accum_degree, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+ 
    coord_flip()
  
}



plot_btw0 <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}

plot_btw_2_trt0 <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+ 
    coord_flip()
  
}



plot_btw <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}



plot_btw_2_trt <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Betweenness = sum(node.betw), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Betweenness, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}





plot_module <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}


plot_module_2_trt <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    # geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)
  # + coord_flip()
  
}


plot_module0 <- function(df1, df2, df3, df4, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Ctrl') 
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'Cd') 
  c <- df3 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  d <- df4 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b, c, d)
  df_merge <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)
  df_merge$trt <- factor(df_merge$trt, levels = c('Ctrl', 'Cd', 'M', 'M+Cd'))
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}


plot_module_2_trt0 <- function(df1, df2, tax_level){
  
  a <- df1 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M')
  b <- df2 %>% group_by({{tax_level}}) %>% 
    summarise(Number_of_module = sum(No.module), trt = 'M+Cd')
  
  df_list <- list(a, b)
  
  
  df_merge <- Reduce(function(x, y) merge(x, y,
                                          all = T), df_list)
  
  ggplot(df_merge, aes({{tax_level}}, Number_of_module, color = trt, group = trt))+
    geom_point(size = 4, alpha = 0.4)+
    geom_line()+
    scale_x_discrete(limits=rev)+
    coord_flip()
  
}
