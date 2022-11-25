#### Function to fit generalized linear model per splicing event, w and w/ independent cis-sQTL of the sGene ####
glm.cis_models <- function(e){
  
  print(paste0("#################  ", e, "  #################"))
  
  # Event PSI residuals ----
  res.e <- psi_residuals[e,]
  # wide to long
  res.e_long <- data.frame(Sample_ID=names(res.e), residuals=as.numeric(res.e), row.names=NULL) #residuals wide to long
  # Merge residuals with metadata -> res.md ----
  res.md <- merge(res.e_long, metadata, by = 'Sample_ID')
  
  # Recover gene:sVariants (isQTL) ----
  # Retrieve the gene
  g <- events.genes[[e]]
  # Retrieve the sVariants
  variants.tissue_g <- gene_variants.list[[g]]
  print(paste0("The event ", e, " is associated to the gene ", g," and has a total of ",length(variants.tissue_g)," sVariants"))
  
  # Retrieve the genotype of the sVariants
  vg_g <- droplevels(vg.tissue_long[vg.tissue_long$variant_id%in%variants.tissue_g,])
  # Filter out individuals if missing genotypes (NA) ----
  vg_g$Individual_ID <- as.character(vg_g$Individual_ID)
  indv_na <- unique(vg_g[is.na(vg_g$genotype),]$Individual_ID)
  indv_in <- unique(vg_g$Individual_ID)[!unique(vg_g$Individual_ID) %in% indv_na]
  vg_g.filtered <- vg_g[vg_g$Individual_ID %in% indv_in,]
  # genotypes long to wide (to have each snps as a column to do the lm())
  vg_g.filtered.wide <- reshape2::dcast(vg_g.filtered[,c('variant_id','Individual_ID','genotype')],
                                        Individual_ID ~ variant_id, value.var="genotype")
  # Subset res.md (residuals + metadata) to have the same individuals
  res.md.filtered <- droplevels(res.md[res.md$Individual_ID %in% indv_in,])
  res.md.filtered <- res.md.filtered[,-which(colnames(res.md.filtered)=='Sample_ID')]
  
  # Merge res.md.filtered (residuals + metadata) & sVariants (vg_g.filtered.wide) by Individual_ID
  df <- merge(res.md.filtered, vg_g.filtered.wide, by = 'Individual_ID')
  df <- df %>% mutate_if(is.character,as.factor)
  variants_def <- colnames(df)[!colnames(df) %in% c(its,'Individual_ID','residuals')]
  
  # Check variance of individual traits (its)  ----
  its_var <- sapply(its, function(it) var(as.numeric(df[[it]])), simplify=F)
  its_in <- names(its_var)[its_var>0]
  its_in.list <- its_var[names(its_var) %in% its_in]
  
  # Check variance of snps ----
  #snps_var <- sapply(variants_def, function(v) var(as.numeric(df[[v]])), simplify=F)
  snps_table <- sapply(variants_def, function(v) table(as.numeric(df[[v]])), simplify=F)
  snps_var <- sapply(variants_def, function(v) sum(table(as.numeric(df[[v]]))>2)==length(table(as.numeric(df[[v]]))) & 
                       length(table(as.numeric(df[[v]])))>1, simplify=F)  # at least 3 donors of each genotype
  snps_in <- names(snps_var)[which(snps_var==T)]
  snps_in.list <- snps_var[names(snps_var) %in% snps_in]
  snps.out <- NA
  #snps_still_dependent <- NA
  
  if(length(snps_in)>0){
    
    # Build models and fit glm ----
    modelA <- paste(c(its_in[its_in!="Ancestry"],snps_in),collapse = '+')
    modelB <- paste(c(its_in[its_in!="Ancestry"],snps_in,"Ancestry"),collapse = '+')
    lm_formula.y <- paste0('residuals')
    lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
    lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
    multi.fit.modelA <- glm(lm_formula.modelA,
                            family = quasibinomial(link = "logit"),
                            control = list(maxit = 100),
                            data=df)
    multi.fit.modelB <- glm(lm_formula.modelB,
                            family = quasibinomial(link = "logit"),
                            control = list(maxit = 100),
                            data=df)
    
    # Check if there are linearly dependent sVariants ----
    if(length(snps_in)>1){
      # cor matrix ----
      df.snps <- df[,colnames(df)%in%snps_in]
      df.snps %>% mutate_if(is.factor, as.numeric) -> df.snps
      correlationMatrix <- cor(df.snps)
  
      # findCorrelation from {caret} package (filter cutoff=0.9 -> out) ----
      #print('Using findCorrelation from caret package to filter linearly dependent sVariants: cutoff>0.9 -> out ')
      highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.9, names = T)
      snps.out <- ifelse(length(highlyCorrelated)==0, NA, paste(highlyCorrelated,collapse = ":"))
          
      if(!is.na(snps.out)){
        print(paste0(snps.out, " is linearly dependent and thus excluded"))
        # Exclude correlated
        snps.f.out <- unlist(strsplit(snps.out, split = ":"))
        snps_in <- snps_in[!snps_in%in%snps.f.out]    
        # Update models
        modelA <- paste(c(its_in[its_in!="Ancestry"],snps_in),collapse = '+')
        modelB <- paste(c(its_in[its_in!="Ancestry"],snps_in,"Ancestry"),collapse = '+')
        lm_formula.y <- paste0('residuals')
        lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
        lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
        multi.fit.modelA <- glm(lm_formula.modelA,
                                family = quasibinomial(link = "logit"),
                                control = list(maxit = 100),
                                data=df)
        multi.fit.modelB <- glm(lm_formula.modelB,
                                family = quasibinomial(link = "logit"),
                                control = list(maxit = 100),
                                data=df)
        # }else{
        #     print(paste0("No highly correlated sVariants"))
        #     correlated_variants <- FALSE
        #     snps_still_dependent <- NA
        #   }
        }
    }else{
      snps.f.out <- c()
    } 
    # e report ----
    # report data frame
    e.df <- data.frame(Tissue = tissue,
                       event = e,
                       ensembl.id = g,
                       indv_total = length(unique(metadata$Individual_ID)),
                       indv_in = length(indv_in),
                       indv_lost = length(indv_na),
                       EUR_total = table(res.md$Ancestry)[['EUR']],
                       AFR_total = table(res.md$Ancestry)[['AFR']],
                       EUR_in = table(res.md[res.md$Individual_ID%in%indv_in,]$Ancestry)[['EUR']],
                       AFR_in = table(res.md[res.md$Individual_ID%in%indv_in,]$Ancestry)[['AFR']],
                       EUR_lost = table(res.md[res.md$Individual_ID%in%indv_na,]$Ancestry)[['EUR']],
                       AFR_lost = table(res.md[res.md$Individual_ID%in%indv_na,]$Ancestry)[['AFR']],
                       num_snps = length(variants_def),
                       num_snps_in = length(snps_in), #ifelse(snps_in==0, 0,length(snps_in)),
                       num_snps_lost =  length(variants_def)-length(snps_in), #ifelse(snps_in==0, length(variants_def), length(variants_def)-length(snps_in)),
                       num_snps_filtered_out =  ifelse(is.na(snps.out), NA, length(snps.f.out))
                       #num_snps_still_dependent = ifelse(!is.na(snps_still_dependent[1]), length(snps_still_dependent), NA)
    )
    rownames(e.df) <- e
    
    # Return results ----
    #if(is.na(snps_still_dependent) || length(snps_still_dependent)==0){
      
      # hier.part ----
      if(length(snps_in) <= 8){
        hier.part.out <- hier.part.mod(y=df[,"residuals"], x=df[,c(its_in[its_in!="Ancestry"],snps_in,"Ancestry")], 
                                       fam = "quasibinomial", link = "logit", gof = "Rsqu")
        if(length(hier.part.out)==1){# & is.na(hier.part.out)){
          r2.traits <- sapply(its_in, function(i) NA)
          r2.sVariants <- NA
        }else{
          hier.part.results <- hier.part.out$IJ
          r2.traits <- sapply(its_in, function(i) hier.part.results[i,1])
          r2.sVariants <- sum(sapply(rownames(hier.part.results)[!rownames(hier.part.results) %in% its_in], function(i) hier.part.results[i,1]))    
        }
      }else{
        r2.traits <- sapply(its_in, function(i) NA)
        r2.sVariants <- NA
      }
      names(r2.sVariants) <- "sVariants"
      
      # output
      res <- list("report_summary" = e.df,
                  "glm.modelA" = multi.fit.modelA,
                  "glm.modelB" = multi.fit.modelB,
                  "df" = df[,colnames(df) %in% c("Individual_ID", "residuals", its_in, snps_in)],
                  "hier.part" = c(r2.traits, r2.sVariants,
                                  "n.sVariants" = length(snps_in),
                                  "EA" = table(res.md[res.md$Individual_ID%in%indv_in,]$Ancestry)[['EUR']],
                                  "AA" = table(res.md[res.md$Individual_ID%in%indv_in,]$Ancestry)[['AFR']]))
      return(res)
      
    # }else{
    #   print(paste0("We're adding sVariants, but it is not possible to build the model because there are still: ", length(snps_still_dependent), " linearly dependent sVariants (", snps_still_dependent, ")"))
    #   return(NA)
    # }
    
  }else{
    print(paste0(g, " has 0 isQTLs with variance"))
    return(NA)
    break
  }
}

  