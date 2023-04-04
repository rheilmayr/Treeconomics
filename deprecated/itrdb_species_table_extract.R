# load required libraries -------------------------------------------------
library(pdftools)
library(naniar)
library(tidyverse)
library(rstudioapi)    
require(dplyr)

# setting up --------------------------------------------------------------
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print(getwd())


# function to read subset of pages
read_single_page <- function(pdf, page){
  tmp <- tempfile()
  on.exit(unlink(tmp))
  tempfile <- pdftools::pdf_subset(pdf, tmp, pages = page)
  pdftools::pdf_text(tmp)
}

# set link to pdf file
pdf_file = 'https://github.com/jasonjb82/general/raw/master/species/Grissino-Mayer%201993%20dendro%20ssp%20list.pdf'

# PARSING -----------------------------------------------------------------
# parse first page
i = 5
df1 <- read_single_page(pdf_file, page = i) %>%
  strsplit("\n") %>%
  as_tibble(.name_repair = make.names) %>%
  slice(16:100) %>%
  mutate(
    CDI = str_sub(X,0,9), 
    code = str_sub(X,10,25),
    text = str_sub(X,26,100)) %>%
  # remove original string
  select(-X) %>%
  slice(2:100) %>%
  # remove white spaces around values
  mutate_all(str_trim) %>%
  #filter(!(CDI=="CD")) %>%
  as.data.frame() %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  fill(CDI,code) %>%
  mutate(page_no = i) %>%
  replace(is.na(.), "") %>%
  mutate(CDI = ifelse(is.na(CDI),0,CDI))

# parse last page
i = 26
df2 <- read_single_page(pdf_file, page = i) %>%
  strsplit("\n") %>%
  as_tibble(.name_repair = make.names) %>%
  slice(3:11) %>%
  mutate(
    CDI = str_sub(X,0,9), 
    code = str_sub(X,10,25),
    text = str_sub(X,26,100)) %>%
  # remove original string
  select(-X) %>%
  slice(2:100) %>%
  # remove white spaces around values
  mutate_all(str_trim) %>%
  #filter(!(CDI=="CD")) %>%
  as.data.frame() %>%
  mutate_if(is.character, list(~na_if(.,""))) %>%
  fill(CDI,code) %>%
  mutate(page_no = i) %>%
  replace(is.na(.), "") %>%
  mutate(CDI = ifelse(is.na(CDI),0,CDI))


# parse rest of the pages
combined_df <- data.frame(CDI=character(),
                          dode=character(),
                          text=character())

for (i in 6:25) {
  
  df3 = read_single_page(pdf_file, page = i) %>%
    strsplit("\n") %>%
    as_tibble(.name_repair = make.names) %>%
    slice(3:100) %>%
    mutate(
      CDI = str_sub(X,0,8), 
      code = str_sub(X,14,19),
      text = str_sub(X,20,100),
      page_no = i,
      text = str_replace(text, "I ",""),
      code = str_replace(code," ","")) %>%
    # remove original string
    select(-X) %>%
    slice(2:100) %>%
    # remove white spaces around values
    mutate_all(str_trim) %>%
    #filter(!(CDI=="CD")) %>%
    as.data.frame() %>%
    mutate_if(is.character, list(~na_if(.,""))) %>%
    mutate(CDI = ifelse(!is.na(code) & is.na(CDI),"",CDI)) %>%
    mutate(code = str_pad(code,4, side = 'right', pad = 'I')) %>%
    fill(CDI,code) 
  
  
  combined_df <- rbind(combined_df,df3)
}

# Merge dataframes
merge_df <- rbind(df1,combined_df,df2) %>%
  mutate(page_no = as.integer(page_no))

# Clean up dataframe
comb_df <- merge_df %>%
  group_by(CDI,code) %>% 
  mutate(text = str_replace(text,'[*]','')) %>%
  mutate(text = str_replace(text,'[.]','')) %>%
  summarise_at(vars(-page_no),funs(paste(., collapse = " ; "))) %>%
  separate(text, c("name", "common_name","common_name1","common_name2","common_name3"), " ; ") %>%
  mutate(common_name = ifelse(!is.na(common_name1),paste0(common_name," ",common_name1),common_name)) %>%
  mutate(common_name = ifelse(!is.na(common_name2),paste0(common_name," ",common_name1," ",common_name2),common_name)) %>%
  mutate(common_name = ifelse(!is.na(common_name2),paste0(common_name," ",common_name1," ",common_name2," ",common_name3),common_name)) %>%
  select(-common_name1,-common_name2,-common_name3) %>%
  #mutate(Genus = stri_extract_first(Name, regex="\\w+")) %>%
  separate(name, into = c("genus", "species","authority"), sep = "\\s",extra="merge") %>%
  mutate(genus = ifelse(grepl('SP$',code),paste0(genus," ",species),genus),
         species = ifelse(grepl('SP$',code),NA,species)) %>%
  mutate(species = ifelse(grepl('^L.',authority),paste0(species," ",authority),species),
         authority = ifelse(grepl('^L.',authority),NA,authority)) %>%
  mutate(authority = str_replace(authority,"[']",'')) %>%
  mutate(genus = ifelse(code == "MIXI",paste0(genus," ",species),genus)) %>%
  mutate(species = ifelse(code == "MIXI",NA,species)) %>%
  mutate(genus = ifelse(code == "CONI",paste0(genus," ",species),genus)) %>%
  mutate(genus = str_replace(genus,'[,]','')) %>% 
  mutate(species = ifelse(code == "CONI",NA,species)) %>%
  mutate(common_name = str_squish(common_name)) %>%
  arrange(code)

# export to csv
write.csv(comb_df,"species_table.csv",row.names = FALSE)