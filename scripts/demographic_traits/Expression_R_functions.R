#### Function to fit  linear model per gene, w and w/ independent cis-eQTL of the eGene ####

lm.cis_models <- function(g){
   
    print(paste0("----  ", g, "  ----"))
    
    # Get gene residual expression across samples ----
    if(is.vector(exprs_residuals)){
        res.g <- exprs_residuals
    }else{
        res.g <- exprs_residuals[rownames(exprs_residuals)==g,]    
    }
    res.g_long <- data.frame(Sample_ID=names(res.g), residuals=res.g, row.names=NULL) #residuals wide to long
    
    # Merge residuals with metadata -> res.md
    res.md <- merge(res.g_long, metadata, by = 'Sample_ID')
    res.md.filtered <- res.md[,-which(colnames(res.md)=='Sample_ID')]
    df <- res.md.filtered
    df <- df %>% mutate_if(is.character,as.factor)
        
    paste0('Fitting the model with cis-eQTLs')
    # 1. Select gene ieQTL
    # 2. Select donors with genotyped ieQTL
    # 3. Select traits with variance
    # 4. Select ieQTLs with variance
    # 5. Run models
    
    # 1. Retrieve the eVariants (ieQTLs in gene) ----
    variants.tissue_g <- gene_variants.list[[g]]
    print(paste0("The gene ", g, " has a total of ",length(variants.tissue_g)," eVariants"))
    
    # Retrieve the genotype of the eVariants
    vg_g <- droplevels(vg.tissue_long[vg.tissue_long$variant_id%in%variants.tissue_g,])
    
    # 2. Filter out individuals if any of the snps is NA ----
    vg_g$Individual_ID <- as.character(vg_g$Individual_ID)
    indv_na <- unique(vg_g[is.na(vg_g$genotype),]$Individual_ID) # individuals with missing genotypes
    indv_in <- unique(vg_g$Individual_ID)[!unique(vg_g$Individual_ID)%in%indv_na] 
    vg_g.filtered <- vg_g[vg_g$Individual_ID%in%indv_in,]

    # genotypes long to wide (to have each snps as a column to do the lm())
    vg_g.filtered.wide <- reshape2::dcast(vg_g.filtered[,c('variant_id','Individual_ID','genotype')],
                                          Individual_ID ~ variant_id, value.var="genotype")
    
    # Subset res.md (residuals + metadata) to have the same individuals
    res.md.filtered <- droplevels(res.md[res.md$Individual_ID%in%indv_in,])
    res.md.filtered <- res.md.filtered[,-which(colnames(res.md.filtered)=='Sample_ID')]
    
    # Merge res.md.filtered (residuals + metadata) & variants (vg_g.filtered.wide) by Individual_ID
    df <- merge(res.md.filtered, vg_g.filtered.wide, by = 'Individual_ID')
    df <- df %>% mutate_if(is.character,as.factor)
    # All genotyped individuals are EUR for gene ENSG00000235615.2 !!!!
    if(length(levels(df$Ancestry))==1){
        print(paste0(g, " only has genotyped individuals of one ancestry"))
        return(NA)
        break
    }
    variants_def <- colnames(df)[!colnames(df)%in%c(its,'Individual_ID','residuals')]
    
    
    # 3. check variance of its (its:individual traits) ----
    its_var <- sapply(its, function(it) var(as.numeric(df[[it]])), simplify=F)
    its_in <- names(its_var)[its_var>0]
    its_in.list <- its_var[names(its_var)%in%its_in]
    
    # 4. check variance of snps before being included in the 'basal' model ----
    #snps_var <- sapply(variants_def, function(v) var(as.numeric(df[[v]])), simplify=F)
    #snps_in <- names(snps_var)[snps_var>0]
    snps_table <- sapply(variants_def, function(v) table(as.numeric(df[[v]])), simplify=F)
    snps_var <- sapply(variants_def, function(v) sum(table(as.numeric(df[[v]]))>2)==length(table(as.numeric(df[[v]]))) & 
                           length(table(as.numeric(df[[v]])))>1, simplify=F)  # at least 3 donors of each genotype
    snps_in <- names(snps_var)[which(snps_var==T)]
    snps_in.list <- snps_var[names(snps_var)%in%snps_in]
    snps.out <- NA
    
    # 5. Run models ----
    if(length(snps_in)>0){
        # Build models and fit lm ----
        modelA <- paste(c(its_in[its_in!="Ancestry"],snps_in),collapse = '+')
        modelB <- paste(c(its_in[its_in!="Ancestry"],snps_in,"Ancestry"),collapse = '+')
        lm_formula.y <- paste0('residuals')
        lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
        lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
        multi.fit.modelA <- lm(lm_formula.modelA,
                               data=df)
        multi.fit.modelB <- lm(lm_formula.modelB,
                               data=df)
        #Check if there are linearly dependent sVariants ----
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
                if(length(snps_in)==0){
                    print(paste0(g, " has 0 eVariants with variance"))
                    return(NA)
                    break
                }
                # Update models
                modelA <- paste(c(its_in[its_in!="Ancestry"],snps_in),collapse = '+')
                modelB <- paste(c(its_in[its_in!="Ancestry"],snps_in,"Ancestry"),collapse = '+')
                lm_formula.y <- paste0('residuals')
                lm_formula.modelA <- paste(c(lm_formula.y, modelA),collapse='~')
                lm_formula.modelB <- paste(c(lm_formula.y, modelB),collapse='~')
                multi.fit.modelA <- lm(lm_formula.modelA,
                                       data=df)
                multi.fit.modelB <- lm(lm_formula.modelB,
                                       data=df)
            }
        }else{
            snps.f.out <- c()
        } 
        # 6. g report ----
        # report data frame
        g.df <- data.frame(Tissue = tissue,
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
                           num_snps_in = length(snps_in),
                           num_snps_lost =  length(variants_def)-length(snps_in),
                           num_snps_filtered_out =  ifelse(is.na(snps.out), NA, length(snps.f.out))
        )
        rownames(g.df) <- g
        
        # 7. hier.part ----
        if(length(snps_in) <= 8){
            hier.part.out <- hier.part.mod(y=df[,"residuals"], x=df[,c(its_in[its_in!="Ancestry"],snps_in,"Ancestry")], fam = "gaussian", gof = "Rsqu")
            if(length(hier.part.out)==1){# & is.na(hier.part.out)){
                r2.traits <- sapply(its_in, function(i) NA)
                r2.eVariants <- NA
            }else{
                hier.part.results <- hier.part.out$IJ
                r2.traits <- sapply(its_in, function(i) hier.part.results[i,1])
                r2.eVariants <- sum(sapply(rownames(hier.part.results)[!rownames(hier.part.results) %in% its_in], function(i) hier.part.results[i,1]))    
            }
        }else{
            r2.traits <- sapply(its_in, function(i) NA)
            r2.eVariants <- NA
        }
        names(r2.eVariants) <- "eVariants"
        
        # 8.output ----
        res <- list("report_summary" = g.df,
                "lm.modelA" = multi.fit.modelA,
                "lm.modelB" = multi.fit.modelB,
                "df" = df[,colnames(df) %in% c("Individual_ID", "residuals", its_in, snps_in)],
                "hier.part" = c(r2.traits, r2.eVariants,
                                "n.eVariants" = length(snps_in),
                                "EA" = table(res.md[res.md$Individual_ID%in%indv_in,]$Ancestry)[['EUR']],
                                "AA" = table(res.md[res.md$Individual_ID%in%indv_in,]$Ancestry)[['AFR']])
                )
        
        return(res)
        
    }else{
        print(paste0(g, " has 0 snps with variance"))
        return(NA)
        break
    }
        
}

