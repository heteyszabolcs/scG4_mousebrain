suppressPackageStartupMessages({
  library("tidyverse")
  library("data.table")
  library("tidymodels")
  library("glue")
  library("themis")
})  

# folders
result_folder = "../results/Seurat/callpeaks_GFPsorted/"
seurat_folder = "../results/Seurat/final/sorted_brain/res0.8/outputs/"
seurat_integration_output = "../results/Seurat/final/sorted_brain/res0.8/integration/outputs/"

# Seurat objects
g4 = readRDS(file = glue("{seurat_folder}Seurat_object.Rds"))
g4_counts = as.matrix(g4@assays$GA@data)

# MOL markers (FindMarker analysis, Marques et al. scRNA-Seq)
markers = read_tsv(glue("{seurat_integration_output}scRNA-Seq_Marques_et_al-FindAllMarkers_output.tsv"),  show_col_types = FALSE)
markers = markers[markers$p_val < 0.05 &
                    markers$avg_log2FC > 0.5, ]
# oligodendrocyte/MOL markers (Zeisel et al., Marques et al.)
zeisel_markers = read_tsv("../data/Zeisel_et_al_2015/Zeisel_et_al_markers.txt", col_names = FALSE, show_col_types = FALSE)
marques_markers = read_tsv("../data/GSE75330/marker_genes.txt", col_names = FALSE, show_col_types = FALSE)
zeisel_markers = zeisel_markers %>% dplyr::filter(X2 == "Oligodendrocyte")
marques_markers = marques_markers %>% dplyr::filter(X2 == "MOL") 
literature_markers = unique(c(marques_markers$X1, zeisel_markers$X1))

# prediction scores of Seurat 
pred_score = readRDS(glue("{seurat_integration_output}g4_cell_label_preds.Rds")) 
pred_score = t(pred_score@data)
ids = rownames(pred_score)
pred_score = as_tibble(pred_score)

# predictions of scBridge integration
pred_bridge = fread("../results/scBridge/output/scbridge_predictions.csv", header = TRUE)
mol_ids_bridge = pred_bridge %>% filter(Prediction == "MOL") %>% pull(V1)

# highly predicted MOL cells
pred_score_mol = pred_score %>% 
  dplyr::select(-max) %>% 
  mutate(cell_id = ids) %>% 
  dplyr::filter(MOL > 0.75)
mol_ids = pred_score_mol %>% pull(cell_id)

length(mol_ids_bridge)
length(mol_ids)
length(intersect(mol_ids_bridge, mol_ids))

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
  top_n(200, wt = avg_log2FC) # top 200 marker
all_mol_markers = mol_markers$gene
mol_markers = unique(mol_markers$gene)[1:20]

opc_markers = markers %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "OPC"), "OPC", cluster)) %>% 
  mutate(cluster = ifelse(str_detect(markers$cluster, "NFOL"), "NFOL", cluster)) %>% 
  dplyr::filter(str_detect(cluster, "OPC")) %>% 
  arrange(desc(avg_log2FC)) %>% 
  top_n(200, wt = avg_log2FC) # top 100 marker
all_opc_markers = opc_markers$gene
opc_markers = unique(opc_markers$gene)[1:20]

## random forest (tidymodels with randomForest), xgboost, logreg
# gene activity scores (G4) of set of genes
# target: MOL type or non-MOL type
random_forest_mol = function(predictors,
                             curve_filename,
                             mol_ids) {
  set.seed(42)
  
  g4 = readRDS(file = glue("{seurat_folder}Seurat_object.Rds"))
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
    # step_corr(all_predictors()) %>%
    # step_center(all_predictors(), -all_outcomes()) %>%
    # step_scale(all_predictors(), -all_outcomes()) %>%
    step_smote(cell_type, over_ratio = 1) %>%
    prep()
  testing = recipe %>%
    bake(testing(trees_split))
  training = juice(recipe)
  
  # for cross validation
  folds = vfold_cv(training, strata = cell_type, v = 5)
  
  ## random forest (tidymodels with ranger)
  rf_spec = rand_forest(mtry = tune(),
                        trees = tune(),
                        min_n = tune()) %>%
    set_engine("randomForest") %>%
    set_mode("classification")
  
  # define workflow
  rf = workflow() %>%
    add_recipe(recipe) %>%
    add_model(rf_spec)
  
  # tuning
  # space-filling designs grid
  rf_grid =
    grid_latin_hypercube(min_n(),
                         mtry(range = c(2, 18)),
                         trees(),
                         size = 80)
  
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
  rf = rand_forest(
    trees = trees,
    min_n = min_n,
    mtry = mtry,
    mode = "classification"
  ) %>%
    set_engine("randomForest") %>%
    fit(cell_type ~ ., data = training)
  
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
  
  print(
    glue(
      "MOL : {as.character(roc)}"
    )
  )
  return(roc)
  
}

# extreme gradient boosting
xgb_mol = function(predictors, mol_ids) {
  
  set.seed(42)
  
  g4 = readRDS(file = glue("{seurat_folder}Seurat_object.Rds"))
  g4_counts = as.matrix(g4@assays$GA@data)
  
  # GA scores of the selected genes
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
    # step_corr(all_predictors()) %>%
    # step_center(all_predictors(), -all_outcomes()) %>%
    # step_scale(all_predictors(), -all_outcomes()) %>%
    step_smote(cell_type, over_ratio = 1) %>% 
    prep()
  testing = recipe %>%
    bake(testing(trees_split)) 
  training = juice(recipe)
  
  # for cross validation
  folds = vfold_cv(training, strata = cell_type, v = 5)
  
  ## extreme gradient boosting
  library("xgboost")
  
  xgb_spec =
    boost_tree(
      tree_depth = tune(),
      min_n = tune(),
      loss_reduction = tune(),
      sample_size = tune(),
      mtry = tune(),
      learn_rate = tune()
    ) %>%
    set_engine("xgboost") %>%
    set_mode("classification") 
  
  xgb = workflow() %>% 
    add_recipe(recipe) %>%
    add_model(xgb_spec)
  
  xgb_grid = grid_latin_hypercube(
    tree_depth(),
    min_n(),
    loss_reduction(),
    sample_size = sample_prop(),
    finalize(mtry(), training),
    learn_rate(),
    size = 80
  )
  
  tune_res = 
    xgb %>% 
    tune_grid(
      resamples = folds, grid = xgb_grid, 
      metrics = metric_set(roc_auc)
    )
  best_params = select_best(tune_res, "roc_auc")
  final_xgb = finalize_workflow(xgb, best_params)
  
  # fit model to training set, evaluates on the test
  last_fit_xgb = last_fit(
    final_xgb,
    split = trees_split,
    metrics = metric_set(recall, precision, f_meas,
                         accuracy, kap,
                         roc_auc, sens, spec)
  )
  return(print(return(last_fit_xgb %>%
                 collect_metrics())))

}

# logistic regression 
logreg_mol = function(predictors,
                      curve_filename,
                      mol_ids) {
  set.seed(42)
  
  g4 = readRDS(file = glue("{seurat_folder}Seurat_object.Rds"))
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
  
  split = initial_split(g4_counts, strata = cell_type, prop = 0.6)
  
  # Create training data
  train = split %>% training() %>% mutate(cell_type = as.factor(cell_type))
  # Create testing data
  test = split %>% testing() %>% mutate(cell_type = as.factor(cell_type))
  
  logreg = logistic_reg() %>%
    # Set the engine
    set_engine("glm") %>%
    # Set the mode
    set_mode("classification") %>%
    # Fit the model
    fit(cell_type ~ ., data = train)
  
  tidy(logreg, exponentiate = TRUE)  %>%
    filter(p.value < 0.05)
  
  # prediction on test
  pred_class = predict(logreg,
                       new_data = test,
                       type = "class")
  pred_proba = predict(logreg,
                       new_data = test,
                       type = "prob")
  results = test %>%
    select(cell_type) %>%
    bind_cols(pred_class, pred_proba)
  
  # model evaluation
  print(conf_mat(results, truth = cell_type,
                 estimate = .pred_class))
  print(accuracy(results, truth = cell_type,
                 estimate = .pred_class))
  print(sens(results, truth = cell_type,
             estimate = .pred_class))
  print(spec(results, truth = cell_type,
             estimate = .pred_class))
  print(precision(results, truth = cell_type,
                  estimate = .pred_class))
  print(f_meas(results, truth = cell_type,
               estimate = .pred_class))
  print(mcc(results, truth = cell_type,
            estimate = .pred_class))
  all_metrics = metric_set(accuracy, sens, spec, precision, recall, f_meas, kap, mcc)
  all_metrics = all_metrics(results,
                            truth = cell_type,
                            estimate = .pred_class)
  
  print(roc_auc(results,
                truth = cell_type,
                ".pred_MOL"))
  
  
  roc = results %>%
    roc_curve(truth = cell_type, ".pred_MOL") %>%
    autoplot + ggtitle("")
  
  ggsave(
    glue("{result_folder}{curve_filename}-logreg_ROC_curve-MOL.pdf"),
    plot = roc,
    width = 6,
    height = 3,
    device = "pdf"
  )
  
  return(all_metrics)
}
  
# run on MOL markers
mol_roc = random_forest_mol(curve_filename = "scBridge_MOLs", 
                            predictors = mol_markers, 
                            mol_ids = mol_ids_bridge)

mol_roc = random_forest_mol(curve_filename = "scBridge_MOLs-upG4_MOL", 
                            predictors = c("Otol1", "Mmp20", "Serpinb3b", "Trim12a", "Glra3",
                                           "Glra3"), 
                            mol_ids = intersect(mol_ids_bridge, mol_ids))

mol_xgb = xgb_mol(predictors = all_mol_markers, mol_ids = mol_ids_bridge)

mol_logreg = logreg_mol(predictors = all_mol_markers, mol_ids = mol_ids_bridge, curve_filename = "scBridge_MOLs")
mol_logreg_negcontrol = logreg_mol(predictors = sample(rownames(g4_counts), size = 200), mol_ids = mol_ids_bridge, curve_filename = "scBridge_random")

# run on OPC markers
opc_roc = random_forest_mol(curve_filename = "highly_pred_MOLs-OPC_markers-", predictors = opc_markers)
random_genes = random_forest_mol(curve_filename = "random_genes", predictors = sample(rownames(g4_counts), size = 20))

# gene activity scores (G4) of set of genes
# target: OPC type or non-OPC type
random_forest_opc = function(predictors,
                             curve_filename
)  {
  set.seed(42)
  
  g4 = readRDS(file = glue("{seurat_folder}Seurat_object.Rds"))
  g4_counts = as.matrix(g4@assays$GA@data)
  
  # GA scores of highly predicted OPC cells
  g4_counts = t(g4_counts)
  ids = rownames(g4_counts)
  g4_counts = as_tibble(g4_counts)
  g4_counts = g4_counts %>%
    mutate(cell_id = ids) %>%
    mutate(cell_type = ifelse(cell_id %in% mol_ids, "OPC", "non-OPC")) %>%
    dplyr::select(cell_type, everything())
  
  g4_counts = g4_counts[, c("cell_type", intersect(predictors, colnames(g4_counts)))]
  
  # splitting
  trees_split = initial_split(g4_counts, strata = cell_type, prop = 0.6)
  
  # recipe
  recipe = training(trees_split) %>%
    recipe(cell_type ~ .) %>%
    # step_corr(all_predictors()) %>%
    # step_center(all_predictors(), -all_outcomes()) %>%
    # step_scale(all_predictors(), -all_outcomes()) %>%
    step_smote(cell_type, over_ratio = 1) %>%
    prep()
  testing = recipe %>%
    bake(testing(trees_split))
  training = juice(recipe)
  
  # for cross validation
  folds = vfold_cv(training, strata = cell_type, v = 5)
  
  ## random forest (tidymodels with ranger)
  rf_spec = rand_forest(mtry = tune(),
                        trees = tune(),
                        min_n = tune()) %>%
    set_engine("randomForest") %>%
    set_mode("classification")
  
  # define workflow
  rf = workflow() %>%
    add_recipe(recipe) %>%
    add_model(rf_spec)
  
  # tuning
  # space-filling designs grid
  rf_grid =
    grid_latin_hypercube(min_n(),
                         mtry(range = c(2, 18)),
                         trees(),
                         size = 80)
  
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
  rf = rand_forest(
    trees = trees,
    min_n = min_n,
    mtry = mtry,
    mode = "classification"
  ) %>%
    set_engine("randomForest") %>%
    fit(cell_type ~ ., data = training)
  
  probs = rf %>%
    predict(testing, type = "prob") %>%
    bind_cols(testing)
  
  roc = probs %>%
    roc_curve(cell_type, ".pred_OPC") %>%
    autoplot + ggtitle("")
  print(roc)
  
  pr = probs %>%
    pr_curve(cell_type, ".pred_OPC") %>%
    autoplot + ggtitle("")
  
  ggsave(
    glue("{result_folder}{curve_filename}-Rf_ROC_curve-OPC.pdf"),
    plot = roc,
    width = 6,
    height = 3,
    device = "pdf"
  )
  
  ggsave(
    glue("{result_folder}{curve_filename}-Rf_PR_curve-OPC.pdf"),
    plot = pr,
    width = 6,
    height = 3,
    device = "pdf"
  )
  
  predict = predict(rf, testing, type = "prob") %>%
    bind_cols(predict(rf, testing)) %>%
    bind_cols(dplyr::select(testing, cell_type)) %>%
    metrics(cell_type, ".pred_OPC", estimate = .pred_class)
  roc = round(predict$.estimate[4], 2)
  
  print(
    glue(
      "OPC : {as.character(roc)}"
    )
  )
  return(roc)
  
}

mol_roc_on_opc_cells = random_forest_opc(curve_filename = "highly_pred_OPCs", predictors = mol_markers)
opc_roc_on_opc_cells = random_forest_opc(curve_filename = "highly_pred_OPCs-OPC_markers-", predictors = opc_markers)
random_genes = random_forest_opc(curve_filename = "random_genes", predictors = sample(rownames(g4_counts), size = 20))
