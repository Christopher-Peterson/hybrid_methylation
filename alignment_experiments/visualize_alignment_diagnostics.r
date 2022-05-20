# Visualize the mapping results of one of the alignments
library(tidyverse)

files = dir("alignment_experiments/.scratch/", full.names = TRUE)
reports = files %>% set_names(basename(files)) %>% map(read_lines)

get_field = function(x, pattern) {
  # browser()
  x %>% str_subset(fixed(pattern)) %>% 
    str_remove(fixed(pattern)) %>% 
    str_remove(fixed(":")) %>% 
    trimws() %>% type.convert()
}

patterns = tribble(~"Name",      ~"pattern",
        "N_seq",     "Sequences analysed in total:",
        "N_unique", "Number of alignments with a unique best hit from the different alignments:",
        "N_unaligned", "Sequences with no alignments under any condition:",
        "N_not_unique", "Sequences did not map uniquely:",
        "N_discarded", "Sequences which were discarded because genomic sequence could not be extracted:",
        "meth_CpG",          "C methylated in CpG context:",
        "meth_CHG",        "C methylated in CHG context:",
        "meth_CHH",        "C methylated in CHH context:",
        "meth_CN_or_CHN",        "C methylated in Unknown context (CN or CHN):",
 ) %>%  with(pattern %>% set_names(Name))

parse_report = function(report) {
  map(patterns, ~get_field(report, .x)) %>% 
    as_tibble()
}

# parse_report(reports[[1]])
results = imap_dfr(reports, ~parse_report(.x) %>% mutate(file = .y))

smry = results %>% 
  transmute(file=file, unique_alignment = N_unique/N_seq,
                      no_alignment = N_unaligned/N_seq,
                      not_mapped_uniquely = N_not_unique/N_seq,
                      discarded = N_discarded/N_seq)
meth_smry = results %>% select(file, starts_with("meth"))
library(cowplot)
fig = smry %>% pivot_longer(-file, names_to = "stat", values_to = "proportion") %>% 
  filter(stat != "discarded") %>% 
  mutate(stat = recode(stat, unique_alignment = "Unique Best",
                       no_alignment = "None",
                       not_mapped_uniquely = "Not Mapped Uniquely")) %>%
  ggplot(aes(x = file, y = proportion, fill = stat)) + 
  geom_col() + 
  theme_cowplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank()) +
  ggtitle("Single End Bismark Alignment, Relaxed") + 
  scale_fill_brewer("Alignment", palette = "Dark2") + 
  theme(legend.text = element_text(size = 8))
    # facet_wrap(~stat) +

ggsave("alignment_experiments/figures/se_relaxed_alignment.png", fig, width = 5, height = 4, dpi = 300)

meth_context_fig = meth_smry %>% pivot_longer(-file, names_to = "Context", names_prefix = "meth_", values_to = "percent") %>% 
  mutate(percent = as.character(percent) %>% str_remove("%") %>% as.numeric()) %>% 
  mutate(Context = recode(Context, CN_or_CHN = "CN or CHN")) %>% 
  ggplot(aes(x = file, y = percent, fill = Context)) + 
  geom_col() + 
  theme_cowplot() +
  theme(legend.position = "bottom",
        axis.text.x = element_blank()) +
  ggtitle("Single End Bismark Alignment, Relaxed") + 
  scale_fill_brewer("Methylation\nContext", palette = "Set2") + 
  theme(legend.text = element_text(size = 8))
ggsave("alignment_experiments/figures/se_relaxed_context.png", meth_context_fig, width = 5, height = 3.5, dpi = 300)
