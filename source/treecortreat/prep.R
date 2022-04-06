library(tidyverse)

sample_meta = readRDS("/dcl02/hongkai/data/covid/data/current/meta.rds")

sample_meta = sample_meta %>%
  mutate(sample = `Library Sample Code`,
         study = gsub('(.*)-(.*)-(.*)', '\\3', sample),
         age = Age,
         sex = Sex,
         severity = type,
         batch = study) %>%
  select(sample, study, age, sex, severity, batch) %>%
  unique()

print(head(sample_meta))
saveRDS(sample_meta, "/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/sample_meta.rds")


clu = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/cluster.rds")
celltype = readRDS("/dcl02/hongkai/data/covid/data/current/pbmc/celltype.rds")
cell_meta = data.frame(celltype = clu) %>%
  mutate(barcode = names(clu),
         sample = gsub('(.*):(.*)', '\\1', barcode))

new_label = data.frame(old = names(celltype),
                       full = celltype) %>%
  mutate(new = gsub('(.*):(.*)', '\\2', full)) %>%
  mutate(large = case_when(
    new == 'T/NK cell' ~ 'T',
    new == 'Monocyte' ~ 'Mono',
    new == 'B cell'~ 'B',
    T ~ new
  )) %>%
  mutate(new = paste0(large,'_', old))

print(head(new_label))
saveRDS(new_label, "/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/new_label.rds")


cell_meta = cell_meta %>%
  mutate(old = celltype) %>%
  inner_join(new_label %>% select(old, new)) %>%
  mutate(celltype = new) %>%
  select(-new)
print(head(cell_meta))
saveRDS(cell_meta, "/dcl02/hongkai/data/rli/covid/wrapup/treecortreat/data/cell_meta.rds")