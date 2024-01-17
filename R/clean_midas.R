library(seqinr)
library(stringr)
library(stringi)
library(tidyverse)

midas_seq <- read.fasta('datbase/SINTAX fa file MiDAS 3.6.fa',
                        forceDNAtolower = F, as.string = T)
midas_names <- names(midas_seq)
tmp <-  midas_names %>% str_split( pattern = ";", n = 3,simplify = T)
raw_names <- str_split_fixed(tmp[, 2], pattern = ",", 7) %>% 
  data.frame %>% mutate_all(str_remove, pattern = "([^:]+)") %>% 
  mutate_all(str_remove, pattern = ":") %>% 
  unite(col = "clean_taxonomy", X1:X7, sep = ";", remove = F) %>% 
  mutate(clean_taxonomy = paste0(clean_taxonomy, ";") )
names(midas_seq) <- raw_names$clean_taxonomy
write.fasta(midas_seq, names = names(midas_seq),
            file.out = "datbase/midas_dada2.fa")
