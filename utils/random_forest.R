suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("tidymodels")
  library("glue")
})  

# folder
result_folder = "../results/Seurat/callpeaks_GFPsorted/"

# Seurat objects
g4 = readRDS(file = "../results/Seurat/final/sorted_brain/res0.8/outputs/Seurat_object.Rds")
g4_counts = as.matrix(g4@assays$GA@data)

# MOL markers (FindMarker analysis, Marques et al. scRNA-Seq)
markers = read_tsv("../results/Seurat/final/sorted_brain/res0.1/integration/outputs/scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv")
markers = markers[markers$p_val < 0.05 &
                    markers$avg_log2FC > 0.5, ]
# oligodendrocyte/MOL markers (Zeisel et al., Marques et al.)
zeisel_markers = read_tsv("../data/Zeisel_et_al_2015/Zeisel_et_al_markers.txt", col_names = FALSE)
marques_markers = read_tsv("../data/GSE75330/marker_genes.txt", col_names = FALSE)
zeisel_markers = zeisel_markers %>% dplyr::filter(X2 == "Oligodendrocyte")
marques_markers = marques_markers %>% dplyr::filter(X2 == "MOL") 
literature_markers = unique(c(marques_markers$X1, zeisel_markers$X1))

# prediction scores of Seurat 
pred_score = readRDS("../results/Seurat/final/sorted_brain/res0.8/integration/outputs/g4_cell_label_preds.Rds")
pred_score = t(pred_score@data)
ids = rownames(pred_score)
pred_score = as_tibble(pred_score)

# highly predicted MOL cells
pred_score_mol = pred_score %>% 
  dplyr::select(-max) %>% 
  mutate(cell_id = ids) %>% 
  dplyr::filter(MOL > 0.75)
mol_ids = pred_score_mol %>% pull(cell_id)

# highly predicted OPC cells
pred_score_opc = pred_score %>% 
  dplyr::select(-max) %>% 
  mutate(cell_id = ids) %>% 
  dplyr::filter(OPC > 0.25)
opc_ids = pred_score_opc %>% pull(cell_id)

# scRNA-Seq markers
mol_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "MOL"), "MOL", cluster)) %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "NFOL"), "NFOL", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "MOL")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(100, wt = avg_log2FC) # top 100 marker
all_mol_markers = mol_markers$gene
mol_markers = unique(mol_markers$gene)[1:20]

opc_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "OPC"), "OPC", cluster)) %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "NFOL"), "NFOL", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "OPC")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(100, wt = avg_log2FC) # top 100 marker
opc_markers = unique(opc_markers$gene)[1:20]

## random forest (tidymodels with randomForest)
# gene activity scores (G4) of set of genes
# target: MOL type or non-MOL type
random_forest_mol = function(predictors,
                             curve_filename
) {
  
  set.seed(42)
  
  g4 = readRDS(file = "../results/Seurat/final/sorted_brain/res0.8/outputs/Seurat_object.Rds")
  g4_counts = as.matrix(g4@assays$GA@data)
  
  # GA scores of highly predicted MOL cells
  g4_counts = t(g4_counts)
  ids = rownames(g4_counts)
  g4_counts = as_tibble(g4_counts)
  g4_counts = g4_counts %>% 
    mutate(cell_id = ids) %>% 
    mutate(cell_type = ifelse(cell_id %in% mol_ids, "MOL", "non-MOL")) %>% 
    dplyr::select(cell_type, everything())
  
  g4_counts = g4_counts[, c("cell_type", intersect(predictors, colnames(g4_counts)))]
  
  # splitting
  trees_split = initial_split(g4_counts, strata = cell_type, prop = 0.6)
  
  # recipe
  recipe = training(trees_split) %>%
    recipe(cell_type ~ .) %>%
    #step_corr(all_predictors()) %>%
    step_center(all_predictors(), -all_outcomes()) %>%
    step_scale(all_predictors(), -all_outcomes()) %>%
    prep()
  testing = recipe %>%
    bake(testing(trees_split)) 
  training = juice(recipe)
  
  # for cross validation
  folds = vfold_cv(training, strata = cell_type, v = 5)
  
  ## random forest (tidymodels with ranger) 
  rf_spec = rand_forest(
    mtry = tune(), trees = tune(), min_n = tune()
  ) %>%
    set_engine(
      "ranger", num.threads = 7, importance = "impurity"
    ) %>%
    set_mode("classification")
  
  # define workflow
  rf = workflow() %>% 
    add_recipe(recipe) %>%
    add_model(rf_spec)
  
  # space-filling designs grid
  rf_grid = 
    grid_latin_hypercube(
      min_n(), 
      mtry(range = c(2, 18)), 
      trees(), 
      size = 80)
  
  # tuning
  tune_res = 
    rf %>% 
    tune_grid(
      resamples = folds, grid = rf_grid, 
      metrics = metric_set(roc_auc)
    )
  
  min_n = show_best(x = tune_res, metric = "roc_auc", n = 1)$min_n
  mtry = show_best(x = tune_res, metric = "roc_auc", n = 1)$mtry
  trees = show_best(x = tune_res, metric = "roc_auc", n = 1)$trees
  
  ## random forest (tidymodels with randomForest)
  rf = rand_forest(trees = trees, min_n = min_n, mtry = mtry, mode = "classification") %>%
    set_engine("randomForest") %>%
    fit(cell_type ~ ., data = training)
  
  # rf_out = rf %>%
  #   predict(testing) %>%
  #   bind_cols(testing) %>%
  #   glimpse()
  # 
  # rf %>%
  #   predict(testing) %>%
  #   bind_cols(testing) %>%
  #   metrics(truth = cell_type, estimate = .pred_class)
  
  probs = rf %>%
    predict(testing, type = "prob") %>%
    bind_cols(testing)
  
  roc = probs %>%
    roc_curve(cell_type, ".pred_MOL") %>%
    autoplot + ggtitle("")
  print(roc)
  
  pr = probs %>%
    pr_curve(cell_type, ".pred_MOL") %>%
    autoplot + ggtitle("")
  print(pr)
  
  ggsave(
    glue("{result_folder}{curve_filename}-Rf_ROC_curve-MOL.pdf"),
    plot = roc,
    width = 6,
    height = 3,
    device = "pdf"
  )
  
  ggsave(
    glue("{result_folder}{curve_filename}-Rf_PR_curve-MOL.pdf"),
    plot = pr,
    width = 6,
    height = 3,
    device = "pdf"
  )
  
  predict = predict(rf, testing, type = "prob") %>%
    bind_cols(predict(rf, testing)) %>%
    bind_cols(dplyr::select(testing, cell_type)) %>%
    metrics(cell_type, ".pred_MOL", estimate = .pred_class)
  roc = round(predict$.estimate[4], 2)
  
  print(glue("MOL : {as.character(roc)}"))
  return(roc)
  
}

mol_roc = random_forest_mol(curve_filename = "highly_pred_MOLs", predictors = mol_markers)
opc_roc = random_forest_mol(curve_filename = "highly_pred_MOLs-OPC_markers-", predictors = opc_markers)
random_genes = random_forest_mol(curve_filename = "random_genes", predictors = sample(rownames(g4_counts), size = 20))

# gene activity scores (G4) of set of genes
# target: OPC type or non-OPC type
random_forest_opc = function(predictors,
                             curve_filename
) {
  
  set.seed(42)
  
  g4 = readRDS(file = "../results/Seurat/final/sorted_brain/res0.8/outputs/Seurat_object.Rds")
  g4_counts = as.matrix(g4@assays$GA@data)
  
  # GA scores of highly predicted OPC cells
  g4_counts = t(g4_counts)
  ids = rownames(g4_counts)
  g4_counts = as_tibble(g4_counts)
  g4_counts = g4_counts %>% 
    mutate(cell_id = ids) %>% 
    mutate(cell_type = ifelse(cell_id %in% opc_ids, "OPC", "non-OPC")) %>% 
    dplyr::select(cell_type, everything())
  
  g4_counts = g4_counts[, c("cell_type", intersect(predictors, colnames(g4_counts)))]
  
  # splitting
  trees_split = initial_split(g4_counts, strata = cell_type, prop = 0.6)
  
  # recipe
  recipe = training(trees_split) %>%
    recipe(cell_type ~ .) %>%
    #step_corr(all_predictors()) %>%
    step_center(all_predictors(), -all_outcomes()) %>%
    step_scale(all_predictors(), -all_outcomes()) %>%
    prep()
  testing = recipe %>%
    bake(testing(trees_split)) 
  training = juice(recipe)
  
  # for cross validation
  folds = vfold_cv(training, strata = cell_type, v = 5)
  
  ## random forest (tidymodels with ranger) 
  rf_spec = rand_forest(
    mtry = tune(), trees = tune(), min_n = tune()
  ) %>%
    set_engine(
      "ranger", num.threads = 7, importance = "impurity"
    ) %>%
    set_mode("classification")
  
  # define workflow
  rf = workflow() %>% 
    add_recipe(recipe) %>%
    add_model(rf_spec)
  
  # space-filling designs grid
  rf_grid <- 
    grid_latin_hypercube(
      min_n(), 
      mtry(range = c(2, 18)), 
      trees(), 
      size = 80)
  
  # tuning
  tune_res = 
    rf %>% 
    tune_grid(
      resamples = folds, grid = rf_grid, 
      metrics = metric_set(roc_auc)
    )
  
  min_n = show_best(x = tune_res, metric = "roc_auc", n = 1)$min_n
  mtry = show_best(x = tune_res, metric = "roc_auc", n = 1)$mtry
  trees = show_best(x = tune_res, metric = "roc_auc", n = 1)$trees
  
  ## random forest (tidymodels with randomForest)
  rf = rand_forest(trees = trees, min_n = min_n, mtry = mtry, mode = "classification") %>%
    set_engine("randomForest") %>%
    fit(cell_type ~ ., data = training)
  
  # rf_out = rf %>%
  #   predict(testing) %>%
  #   bind_cols(testing) %>%
  #   glimpse()
  
  # rf %>%
  #   predict(testing) %>%
  #   bind_cols(testing) %>%
  #   metrics(truth = cell_type, estimate = .pred_class)
  
  probs = rf %>%
    predict(testing, type = "prob") %>%
    bind_cols(testing)
  
  roc = probs %>%
    roc_curve(cell_type, ".pred_OPC") %>%
    autoplot + ggtitle("")
  print(roc)
  
  ggsave(
    glue("{result_folder}{curve_filename}-Rf_ROC_curve_OPC.pdf"),
    plot = roc,
    width = 6,
    height = 3,
    device = "pdf"
  )
  
  predict = predict(rf, testing, type = "prob") %>%
    bind_cols(predict(rf, testing)) %>%
    bind_cols(dplyr::select(testing, cell_type)) %>%
    metrics(cell_type, ".pred_OPC", estimate = .pred_class)
  roc = round(predict$.estimate[4], 2)
  
  print(glue("OPC : {as.character(roc)}"))
  return(roc)
  
}

mol_roc_on_opc_cells = random_forest_opc(curve_filename = "highly_pred_OPCs", predictors = mol_markers)
opc_roc_on_opc_cells = random_forest_opc(curve_filename = "highly_pred_OPCs-OPC_markers-", predictors = opc_markers)
random_genes = random_forest_opc(curve_filename = "random_genes", predictors = sample(rownames(g4_counts), size = 20))
