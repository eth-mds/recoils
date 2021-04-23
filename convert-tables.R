library(arrow)
library(fst)
library(ricu)

list.files(file.path(data_dir(), "covid19", "single_timestamp"))

if (!dir.exists(file.path(data_dir(), "covid19"))) 
  dir.create(file.path(data_dir(), "covid19"))

if (!dir.exists(file.path(data_dir(), "covid19"))) 
  dir.create(file.path(data_dir(), "covid19"))

convert_names <- c(
  'episodes', 'admissions',
  'medications', 'intubations', 'range_measurements',
  'diagnoses', 'parameters', 'comorbidities', 'outcomes'
)

convert_names <- c(
  convert_names,
  gsub(".parquet", "", paste0("single_timestamp/", list.files("single_timestamp")))
)

#convert_names <- grep("part_", list.files(), value = T)
#convert_names <- gsub(".parquet", "", convert_names)

for (tab_name in convert_names) {
  
  if (file.exists(paste0(tab_name, ".parquet"))) {
    
    tbl <- read_parquet(paste0(tab_name, ".parquet"))
    file.remove(paste0(tab_name, ".parquet"))
    
    write_fst(tbl, path = file.path(data_dir(), "covid19", paste0(tab_name, ".fst")))
    
  }
  
  print(tab_name)
  
}
