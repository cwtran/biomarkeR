getHR <- function(formula_text,data_subset) {
  #get hazard ratio info
  fit_survdiff <- survdiff(as.formula(formula_text), data = data_subset)
  fit_coxph <- coxph(as.formula(formula_text), data = data_subset, ties = "exact")
  pval = pchisq(fit_survdiff$chisq, length(fit_survdiff$n) - 1, lower.tail = FALSE)
  pval.z = summary(fit_coxph)[['coefficients']][,'Pr(>|z|)']
  hr = if (is.null(fit_coxph) || all(is.na(fit_coxph$coefficients))) NA else (summary(fit_coxph)$conf.int)[1:nrow(summary(fit_coxph)$conf.int), 1]
  hr_ci = if (is.null(fit_coxph) || all(is.na(fit_coxph$coefficients))) c(NA, NA) else (summary(fit_coxph)$conf.int)[1:nrow(summary(fit_coxph)$conf.int), 3:4]
  ref_group = fit_coxph$xlevels[[1]][1]
  
  hr_string = if (is.null(dim(hr_ci)) == T) {
    # only 2 groups
    paste0(round(hr,digits=4)," [",round(as.numeric(hr_ci[1]),digits=4),"-",round(as.numeric(hr_ci[2]),digits=4),"]")
  } else { 
    # more than 2 groups
    paste0(round(hr,digits=4)," [",round(as.numeric(hr_ci[,1]),digits=4),"-",round(as.numeric(hr_ci[,2]),digits=4),"]")
  }
  
  hr_string2 = if (is.null(dim(hr_ci)) == T) {
    paste0(round(hr,digits=2)," [",round(as.numeric(hr_ci[1]),digits=2),"-",round(as.numeric(hr_ci[2]),digits=2),"]")
  } else {
    paste0(round(hr,digits=2)," [",round(as.numeric(hr_ci[,1]),digits=2),"-",round(as.numeric(hr_ci[,2]),digits=2),"]")
  }
  
  hrDF=data.frame( hr_value=hr,
                    ci_low=ifelse(is.null(dim(hr_ci)), as.numeric(hr_ci[1]), as.numeric(hr_ci[, 1])),
                    ci_high=ifelse(is.null(dim(hr_ci)), as.numeric(hr_ci[2]), as.numeric(hr_ci[, 2])),
                    pval=pval,
                    pval.z = signif(pval.z, 2),
                  hr_string=hr_string,
                  hr_string2=hr_string2,
                  hr_string3 = paste0(hr_string2, " p=", signif(pval.z, 2)),
                    ref_group = ref_group,
                    stringsAsFactors=F)
  
  # cleanup rownames
  if (nrow(hrDF) == 1) {
    rownames(hrDF) <- fit_coxph$xlevels[[1]][2]
  } else {
    temp_names <- gsub(trimws(strsplit(as.character(formula_text), "~")[[1]][2]), "", rownames(hrDF))
    rownames(hrDF) <- temp_names
  }
  
  
  hrDF <- hrDF %>% tibble::rownames_to_column("group")
  #hrDF$group <- gsub("low", "L", hrDF$group)
  hrDF$group <- gsub("moderate", "mod", hrDF$group)
  #hrDF$group <- gsub("high", "H", hrDF$group)
  
  return(hrDF)
}


# OS / PFS generic function
SurvivalbySig <- function(covar, datasource, addHR, censorcol, survtimecol, survtype, legend_title = NULL, legend_labels = NULL) {
  
  # setting plotting function to allow saving
  # grid.draw.ggsurvplot <- function(x){
  #  survminer:::print.ggsurvplot(x, newpage = FALSE)
  # }
  
  # define variables
  surv.plot.list <- list()
  hr.results.list <- list()
  covar <- covar # covariate/group we want to explore
  ds <- datasource # dataframe
  addHR <- addHR # addHR summary to plot (T or F)
  ccol <- censorcol # which column is the censoring column (OS_CNSR, PFS_CNSR)
  survtime <- survtimecol # which column is the survival time column (OS, PFS)
  survtype <- survtype # OS or PFS
  
  if(!is.null(legend_title)) {
    legend_title <- legend_title
  }
  
  if(!is.null(legend_title)) {
    legend_labels <- legend_labels
  }
  
  if(!is.null(legend_title) & !is.null(legend_labels)) {
      for (d in 1:length(covar)) {
        temp_cv <- covar[d]
        maxOS <- max(ds %>% select(!!sym(survtime)))
        
        surv.plot.list[[d]] <- ggsurvplot(
          fit = surv_fit(as.formula(paste0("Surv(", survtime, ", ", ccol, ") ~ ", temp_cv)), data = ds),
          xlab = "Months", 
          xlim = c(0, maxOS), # dynamically vary by max OS in data
          ylab = paste0(survtype),
          legend.title = paste(legend_title), 
          legend.labs = legend_labels,
          surv.median.line = "v",
          palette = "npg",
          risk.table = T,
          fontsize = 3.5,
          pval = F,
          conf.int = F,
          break.time.by = 3,
          risk.table.y.text.col = T,
          tables.theme = theme_cleantable(),
          tables.height = 0.2,
          tables.y.text = FALSE)
      }
    } else {
      # custom legend title and labels
      for (d in 1:length(covar)) {
        temp_cv <- covar[d]
        maxOS <- max(ds %>% select(!!sym(survtime)))
        
        surv.plot.list[[d]] <- ggsurvplot(
          fit = surv_fit(as.formula(paste0("Surv(", survtime, ", ", ccol, ") ~ ", temp_cv)), data = ds),
          xlab = "Months", 
          xlim = c(0, maxOS), # dynamically vary by max OS in data
          ylab = paste0(survtype),
          # legend.title = paste(legend_title), # use default
          # legend.labs = legend_labels,
          surv.median.line = "v",
          palette = "npg",
          risk.table = T,
          fontsize = 3.5, # risk table font
          pval = F,
          conf.int = F,
          break.time.by = 3,
          risk.table.y.text.col = T,
          tables.theme = theme_cleantable(),
          tables.height = 0.2,
          tables.y.text = FALSE)
    }
  }
  
  # increase legend text size
  surv.plot.list[[d]]$plot <- surv.plot.list[[d]]$plot + theme(legend.text = element_text(size = 14, color = "black"), 
                                                           legend.title = element_text(size = 14, color = "black", face="bold"))
  
  
  if (addHR == T) {
    # getHR
    surv_formula <- paste0("Surv(", survtime, ", ", ccol, ") ~ ", temp_cv)
    hr.results.list[[d]] = getHR(surv_formula, data = ds)
    
  

    # add HR [95% CI] and p value to top right of plot
    for (t in 1:length(hr.results.list[[d]]$hr_string2)) {
      hr_string_title <- "Group          HR (95% CI)"
      if (t == 1) {
      hr_string_final <- paste0(hr_string_title, "\n", hr.results.list[[d]]$ref_group[1], "          ", "Reference")
      hr_string_final <- paste0(hr_string_final, "\n", hr.results.list[[d]]$group[t], "     ", hr.results.list[[d]]$hr_string3[t])
      } else {
      hr_string_final <- paste0(hr_string_final, "\n", hr.results.list[[d]]$group[t], "     ", hr.results.list[[d]]$hr_string3[t])
      }
    }
  
    #df_txt <- data.frame(x=label.custom.x, y=0.92, lab=paste0("HR = ", hr_info$hr_string2,"\n p=", signif(hr_info$pval,digits=2)))
    # df_txt <- data.frame(x=(maxOS*0.7), y=0.85, lab=paste0(hr_string_final,"\n p=", signif(hr.results.list[[d]]$pval,digits=2)))
    df_txt <- data.frame(x=(maxOS*0.7), y=0.85, lab=paste0(hr_string_final))
    
    ngroups = nrow(hr.results.list[[d]])+1
    
    surv.plot.list[[d]]$plot$layers[[5]] <- layer(geom="text", position="identity", stat="identity", 
                              mapping=aes(x=x, y=y, label=lab), data=df_txt,
                              params=list(size=5, color = "black"))
  }
  
  
  # calculate number of rows needed for plot (if more than 1 plot)
  if (length(covar) > 1) {
    adj.row <- ceiling(length(covar) / 4)
    surv.all.plot <- arrange_ggsurvplots(surv.plot.list, print = TRUE, ncol = ifelse(length(covar) < 4, length(covar), 4), nrow = adj.row)
    return(surv.all.plot) 
    
  } else {
    return(surv.plot.list[[1]]) # only 1 plot
  }
    
}
