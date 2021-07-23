source('/workspace/workspace_pipelines/mace_paper_2020/00_packages_functions.R', echo=TRUE)

#load the data
cds_497<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/mace_aug2019/output_mace_aug2019/mace_aug2019_497", barcode_filtered = TRUE)
cds_505<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/mace_aug2019/output_mace_aug2019/mace_aug2019_505", barcode_filtered = TRUE)
cds_513<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/mace_aug2019/output_mace_aug2019/mace_aug2019_513", barcode_filtered = TRUE)
cds_520<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/mace_aug2019/output_mace_aug2019/mace_aug2019_520", barcode_filtered = TRUE)

colData(cds_497)$pt<-"pt_497"
colData(cds_505)$pt<-"pt_505"
colData(cds_513)$pt<-"pt_513"#exclude this one to be consistent with submitted paper
colData(cds_520)$pt<-"pt_520"

#combine the cds's we will use
cds_list<-list(cds_497,cds_505,cds_520)
cds<-combine_cds(cds_list = cds_list, keep_all_genes = TRUE)

#trim off uninformative genes
cds_trimmed<-cds[substr(rowData(cds)$gene_short_name,1,2)!="RP",]

## Normalize and pre-process the data
cds_trimmed<-preprocess_cds(cds_trimmed, num_dim = 100)

## Reduce dimensionality and visualize cells
cds_trimmed<-reduce_dimension(cds_trimmed, cores = 39)

# Group cells into clusters
cds_trimmed<-cluster_cells(cds_trimmed,k = 20)
plot_cells(cds_trimmed, group_cells_by = "cluster", color_cells_by = "cluster", group_label_size = 5)

# manually assign celltypes to clusters
cds_trimmed$cluster<-monocle3::clusters(cds_trimmed)
cds_trimmed$cluster_assignment<-dplyr::recode(cds_trimmed$cluster, "1" = "MDSC",
                              "3" = "Naive/Memory T",
                              "2" = "CD8",
                              "4" = "NK",
                              "5" = "B",
                              "6" = "high-MT T",
                              "7" = "CD16+ Monocytes",
                              "8" = "Granulocyte",
                              "9" = "Platelet",
                              "10" = "pDC",
                              "11" = "Proliferating T")

# identify marker genes and output data
marker_test_res <- top_markers(cds_trimmed, group_cells_by="cluster_assignment", reference_cells=1000, cores=39,genes_to_test_per_group = 100)
marker_test_res %>% dplyr::rename(cluster = cell_group) %>% arrange(cluster) %>% write_csv(path = "data_out/cluster_top_markers.csv")

# general group stats
## cell counts by cluster
tbl_df(colData(cds_trimmed)) %>% 
  group_by(cluster_assignment) %>% 
  summarise(n = n()) %>% 
  mutate(running_total = cumsum(n)) %>%
  write_csv(path = "data_out/cluster_assignment_stats.csv")

## cell counts by patient
tbl_df(colData(cds_trimmed)) %>% 
  group_by(pt) %>% 
  summarise(n = n()) %>%
  write_csv(path = "data_out/patient_stats.csv")

save.image.pigz(file = "mace_paper_2020.RData", n.cores = 39)
