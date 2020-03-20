load.pigz("mace_paper_2020.RData")

demo_genes<-c("CD14", "IL7R", "CD8A", "KLRB1", "CD79A", "FCGR3A", "PF4", "IL3RA","PCNA")

plot_cells_alt(cds_trimmed, genes = demo_genes,label_cell_groups = F) %>% 
  save_plot(filename = "plots_out/marker_genes.pdf", base_width = 5, base_height = 4)

save.image.pigz("mace_paper_2020.RData",n.cores = 39)
