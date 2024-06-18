
design<-model.matrix(~0 + immuno_only_wide$Age + immuno_only_wide$Sex + immuno_only_wide$LDH_elevated + immuno_only_wide$staging_groups + immuno_only_wide$OS_months)
colnames(design) <- c('Age', 'Sex','LDHelevatedNA','LDHelevated_No','LDHelevated_yes','staging_m1cd','OS_months')

fit <- immuno_only_wide %>%
  select(where(is.numeric)) %>% 
  select(-Age,-Sex,-thaw_freeze_cycles,-progress_days,-progress_months,-death,-OS_months) %>% 
  t() %>% 
  lmFit(design=design)

cont <- makeContrasts(OS_months, levels = colnames(design))

# apply contrast
contrast_fit<-contrasts.fit(fit, cont)


# apply empirical Bayes smoothing to the SE
ebays_fit<-eBayes(contrast_fit)

# summary
print(summary(decideTests(ebays_fit)))

# extract DE results
DE_OS_months_immuno_only <-topTable(ebays_fit, n=ncol(immuno_only_wide), adjust.method="fdr", confint=TRUE)

# label as up or down regulated based on values 
DE_OS_months_immuno_only <- DE_OS_months_immuno_only %>% 
  mutate(Expression = case_when(logFC >= log(1)& adj.P.Val < 0.2 ~ "Up-regulated",
                                logFC <= -log(1) & adj.P.Val < 0.2 ~ "Down-regulated",
                                TRUE ~ "Unchanged"))
top <- 30
top_up <- DE_OS_months_immuno_only %>% 
  filter(Expression == 'Up-regulated') %>% 
  arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
  head(top)
top_down <- DE_OS_months_immuno_only %>% 
  filter(Expression == 'Down-regulated') %>% 
  arrange(desc(-log10(adj.P.Val)), desc(abs(logFC))) %>% 
  head(top)

# idenitfy top genes
top_genes <- rbind(top_up,top_down)

# make volcano plot
volcano_plot <- DE_OS_months_immuno_only %>% 
  ggplot(aes(logFC, -log10(adj.P.Val), color = ifelse(adj.P.Val < 0.2 & logFC > 0, "red",
                                                      ifelse(adj.P.Val < 0.2 & logFC < 0, "blue", "grey")))) +
  geom_point() +
  scale_color_identity() +  # Use identity scale for manual colors
  theme_bw()

volcano_plot <- volcano_plot +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log10(adj.P.Val), label = rownames(top_genes)),
                   size = 2) +
  ggtitle('DE_OS_months, immuno only - adjusted for age, sex, LDH, staging')