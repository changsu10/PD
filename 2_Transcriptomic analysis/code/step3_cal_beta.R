library(parallel)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)

### read
trait_df=read.csv('../result/planB/trait_random_slope.csv', check.names = FALSE, row.names = 1)
gene_df=read.csv('../result/planB/gene_random_slope.csv', check.names = FALSE, row.names = 1)

### align sample ids
colnames(trait_df)=paste0('PP-',colnames(trait_df))
gene_df=gene_df[,colnames(trait_df)]

### fit glm
get_glm <- function(gene, trait, gene_data, trait_data) {
    tryCatch({
        formula <- as.formula(paste(trait, "~", gene))
        data <- data.frame(y = t(trait_data[trait, ]), x = t(gene_data[gene, ]))
        model <- glm(formula, data = data)

        beta <- coef(model)[2]  # Extract beta coefficient
        p_value <- summary(model)$coefficients[2, 4]  # Extract p-value
        return(c(beta = as.numeric(beta), p_value = p_value))

  }, error = function(e) {#some gene's beta are all 0 (maybe not converge?)
    message(paste("Error at gene", gene, '& trait',trait,":", e$message))
    return(NULL)  # Return a placeholder or log the error
  })
}

# Reserve one core for the system
no_cores <- detectCores() - 1  
cl <- makeCluster(no_cores)
tasks <- expand.grid(gene = rownames(gene_df), trait = rownames(trait_df))# Create a list of tasks
clusterExport(cl, varlist = c("gene_df", "trait_df", "get_glm", "tasks"))

#Run the tasks in parallel
results <- parLapply(cl, 1:nrow(tasks), function(i) {
  get_glm(tasks$gene[i], tasks$trait[i], gene_df, trait_df)
})

# Stop the cluster
stopCluster(cl)

# Initialize matrices
beta_matrix <- matrix(nrow = nrow(gene_df), ncol = nrow(trait_df), dimnames = list(rownames(gene_df), rownames(trait_df)))
p_value_matrix <- matrix(nrow = nrow(gene_df), ncol = nrow(trait_df), dimnames = list(rownames(gene_df), rownames(trait_df)))

# Fill the matrices
for (i in 1:length(results)) {
  gene <- tasks$gene[i]
  trait <- tasks$trait[i]
  if (length(results[[i]])>0){
    beta_matrix[gene, trait] <- results[[i]]["beta"]
    p_value_matrix[gene, trait] <- results[[i]]["p_value"]
  } else{
    #beta_matrix[gene, trait]=NA
    #p_value_matrix[gene, trait]=NA
  }
}

write.csv(beta_matrix, "../result/planB/beta_result.csv",quote=FALSE)
write.csv(p_value_matrix, "../result/planB/pvalue_result.csv",quote=FALSE)



