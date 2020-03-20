source('/workspace/workspace_pipelines/mace_paper_2020/00_packages_functions.R', echo=TRUE)
load.pigz("mace_paper_2020.RData")

#max plot in the style of plot_cells_alt
gene_list<-c("CD200R1", "DOK1", "DOK2")
agg_df<-plot_cells(cds_trimmed, genes = gene_list ,norm_method = "size_only")[["data"]]
agg_df1<-cumulative_max_expr(extracted_df = agg_df, gene_list = gene_list)

max_plot<-ggplot() + 
  geom_point(data = agg_df1[is.na(agg_df1$max_val), ],
             aes(data_dim_1, data_dim_2), 
             size = 1, 
             stroke = 0.25, 
             color = "grey80",
             shape = 1,
             alpha = 1) + 
  geom_point(data = na.omit(agg_df1, cols = "max_val"), 
             aes(x = data_dim_1, y = data_dim_2,color = max_val), 
             size = 1, 
             stroke = 0,
             na.rm = TRUE,
             alpha = 1,
             shape = 16) + 
  viridis::scale_color_viridis(option = "viridis",
                               name = "max expr", 
                               na.value = "transparent", 
                               end = 0.8, 
                               alpha = 1)+
  labs(x = "UMAP 1", y = "UMAP 2")
max_plot
save_plot(max_plot, filename = "plots_out/max_plot.pdf", base_width = 4, base_height = 3)

# plot each individually
plot_cells_alt(cds_trimmed, genes = "DOK2", cell_size = 1) %>% save_plot(filename = "plots_out/dok2.pdf", base_width = 4, base_height = 3)
plot_cells_alt(cds_trimmed, genes = "DOK1", cell_size = 1) %>% save_plot(filename = "plots_out/dok1.pdf", base_width = 4, base_height = 3)
plot_cells_alt(cds_trimmed, genes = "CD200R1", cell_size = 1) %>% save_plot(filename = "plots_out/cd200r1.pdf", base_width = 4, base_height = 3)             


save.image.pigz("mace_paper_2020.RData",n.cores = 39)
