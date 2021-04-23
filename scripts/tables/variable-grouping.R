
cfg <- config("feature2")
sys <- sapply(cfg, `[[`, "category")

u_sys <- unique(sys)

sort_sys <- lapply(u_sys, 
                   function(s) sapply(cfg[which(sys == s)], `[[`, "full_name"))

names(sort_sys) <- u_sys

groups <- Reduce(
  c,
  lapply(sort_sys, function(x) paste(sapply(x, str_to_title), collapse = ", "))
)

res <- data.frame(System = u_sys, Features = groups)

res

my_doc <- read_docx()

my_doc <- my_doc %>%
  body_add_table(res, style = "table_template")

print(my_doc, target = "~/covid-19/naps-res/SupplementaryTable1.docx")
