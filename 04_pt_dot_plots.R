source('/workspace/workspace_pipelines/mace_paper_2020/00_packages_functions.R', echo=TRUE)
load.pigz("mace_paper_2020.RData")

all_cells_data<-tbl_df(plot_cells(cds = cds_trimmed)[["data"]])

#plot all cells by cluster
all_cells_plot<-ggplot(data = all_cells_data, aes(x = data_dim_1, y = data_dim_2, fill = cluster_assignment, color = cluster_assignment))
all_cells_plot<-all_cells_plot+
  geom_point(size = 0.5, alpha = 0.2, shape = 21)+
  theme(legend.position = "none")+
  labs(x = "UMAP 1", y = "UMAP 2")
all_cells_plot
save_plot(all_cells_plot, filename = "plots_out/all_cells_plot.pdf", base_width = 3.3, base_height = 3)

# plot overlay by patient
## overlay for pt_497
background_for_497<-all_cells_data %>%
  filter(pt != "pt_497")
foreground_for_497<-all_cells_data %>%
  filter(pt == "pt_497")

## overlay for pt_505
background_for_505<-all_cells_data %>%
  filter(pt != "pt_505")
foreground_for_505<-all_cells_data %>%
  filter(pt == "pt_505")

## overlay for pt_520
background_for_520<-all_cells_data %>%
  filter(pt != "pt_520")
foreground_for_520<-all_cells_data %>%
  filter(pt == "pt_520")

#overlay plot function
overlay_plot<-function(foreground_data, background_data) {
  plot<-ggplot()+
    geom_point(data = background_data, aes(x = data_dim_1, y = data_dim_2), color = "grey80", shape = 1, alpha = 0.2, size = 0.5)+
    geom_point(data = foreground_data, aes(x = data_dim_1, y = data_dim_2, color = cluster_assignment, fill = cluster_assignment), shape = 21, alpha = 0.2, size = 0.5)+
    theme(legend.position = "none")+
    labs(x = "UMAP 1", y = "UMAP 2")
  return(plot)
  
}

save_plot(overlay_plot(foreground_data = foreground_for_497, background_data = background_for_497), 
          filename = "plots_out/pt_497.pdf", 
          base_width = 3.3, base_height = 3)

save_plot(overlay_plot(foreground_data = foreground_for_505, background_data = background_for_505), 
          filename = "plots_out/pt_505.pdf", 
          base_width = 3.3, base_height = 3)

save_plot(overlay_plot(foreground_data = foreground_for_520, background_data = background_for_520), 
          filename = "plots_out/pt_520.pdf", 
          base_width = 3.3, base_height = 3)
 
save.image.pigz("mace_paper_2020.RData",n.cores = 39)
