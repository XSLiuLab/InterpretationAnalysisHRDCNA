
### correlation matrix for 80 CNA features

library(corrplot)
library(reshape2)
library(RColorBrewer)
library(ggplot2)

alldata <- readRDS("~/HRD/HRDCNA/data/modeldata/alldata.rds")

cor_matrix <- cor(alldata)
cor_df <- melt(cor_matrix)
cor_df_sorted <- cor_df[order(-abs(cor_df$value)), ]
cor_matrix_sorted <- acast(cor_df_sorted, Var1 ~ Var2, value.var = "value")


ggplot(data = cor_df_sorted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#601787", mid = "white", high = "#A54C11",
                       midpoint = 0, limits = c(-1, 1), name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Matrix")



