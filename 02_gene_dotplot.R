load.pigz("mace_paper_2020.RData")

demo_genes<-c("CD14", "IL7R", "CD8A", "KLRB1", "CD79A", "FCGR3A", "PF4", "IL3RA","PCNA")

#dot plot
gene_dot_plot<-plot_genes_by_group(cds = cds_trimmed, 
                                   markers = demo_genes, 
                                   group_cells_by = "cluster_assignment", 
                                   ordering_type = "maximal_on_diag")
gene_dot_plot<-gene_dot_plot + labs(x = NULL, y = NULL)
gene_dot_plot
save_plot(gene_dot_plot, filename = "plots_out/gene_dot_plot.pdf", base_width = 6, base_height = 5)

save.image.pigz("mace_paper_2020.RData",n.cores = 39)
