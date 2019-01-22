# This script is to do venn diagram to show how the analysis overlap

library(VennDiagram)

# define a function to do the venn diagram since most of the code is the same(currently max 3 things to compare)
draw_venn_diagram <- function(diff_list, fill, size, title){
  venn = venn.diagram(x = diff_list, 
                       filename = NULL,
                       height = 2000,
                       width = 2000, fill = fill,
                       cat.default.pos = "text", 
                       cat.cex = size,
                       main = title);
  png(str_interp('../figs/${title}.png'), width = 4, height = 4, units = 'in', res = 300)
  grid.draw(venn)
  dev.off()
}

# plot a venn diagram to show this overlap relationship (use all of the data, unfiltered result)

diff_list <- list(RF14_cntrl_sig = res14$Locus,
                  RF89_cntrl_sig = res89$Locus)
fill <- c("orange", "light blue")
size  <- rep(0.5, 2)


# did the pooling all KD together myself using IRfinder
pooledKD <- read_tsv('../data/myPooledKD/RFKD_vs_ctrl.txt',comment = '#',col_names = T) %>%
  mutate(Locus =paste(Chr,paste(Start, End, sep = '-'), sep = ':' ) )



# plot a venn diagram to show this overlap relationship
diff_list <- list(RF14_cntrl = res_all_14$Locus,
                  RF89_cntrl = res_all_89$Locus,
                  pooled_KD_all = pooledKD$Locus)
fill <- c("orange", "light blue", 'pink')
size  <- rep(0.5, 3)
venn <- venn.diagram(x = diff_list, 
                     filename = NULL,
                     height = 2000,
                     width = 2000, fill = fill,
                     cat.default.pos = "text", 
                     cat.cex = size,
                     main = "Overlapped between the raw result of the locus between the pooled and separate");
png('../figs/Overlapped between the raw result of the locus between the pooled and separate.png', width = 4, height = 4, units = 'in', res = 300)
grid.draw(venn)
dev.off()

# do the venn diagram of the filtered result
pooledKD_fil <- filterRawResult('../data/myPooledKD/RFKD_vs_ctrl.txt')

# venn
diff_list <- list(RF14_cntrl = res14$Locus,
                  RF89_cntrl = res89$Locus,
                  pooled_KD_fil = pooledKD_fil$Locus)
fill <- c("orange", "light blue", 'pink')
size  <- rep(0.5, 3)
venn <- venn.diagram(x = diff_list, 
                     filename = NULL,
                     height = 2000,
                     width = 2000, fill = fill,
                     cat.default.pos = "text", 
                     cat.cex = size,
                     main = "OL between the filtered result of the locus between the pooled and separate");
png('../figs/OL between the filtered result of the locus between the pooled and separate.png', width = 4, height = 4, units = 'in', res = 300)
grid.draw(venn)
dev.off()


res_all_14 <- raw_df_list[[1]]
res_all_89 <- raw_df_list[[2]]
res14T <- true_df_list[[1]]
res89T <- true_df_list[[2]]
res14F <- false_df_list[[1]]
res89F <- false_df_list[[2]]

# trying to make a venn diagram of what exactly I am comparing
diff_list <- list( 
  RF14_cntrl_SIG_T = res14T$Locus,
  RF89_cntrl_SIG_T = res89T$Locus,
  RF14_cntrl_SIG_F = res14F$Locus,
  RF89_cntrl_SIG_F = res89F$Locus
  
)

fill <- c("#999999", "#E69F00", "#56B4E9", "#009E73")
size  <- rep(0.5, length(fill))
venn <- venn.diagram(x = diff_list, 
                     filename = NULL,
                     height = 2000,
                     width = 2000, fill = fill,
                     cat.default.pos = "text", 
                     cat.cex = size,
                     main = "Overview of sig among the raw results");
png('../figs/Overview of sig among the raw results_all.png', width = 4, height = 4, units = 'in', res = 300)
grid.draw(venn)
dev.off()




