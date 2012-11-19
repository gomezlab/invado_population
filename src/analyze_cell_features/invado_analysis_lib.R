################################################################################
# invado_analysis_lib.R: the R functions associated with processing the
# invaodopodia data
################################################################################

library(Hmisc);

get_non_nan_data <- function(files,min.lifetime=0) {
    data = list();
    
    for (this_file in files) {
        values = as.matrix(read.csv(this_file,header=F));        
        
        for (cell_num in 1:dim(values)[1]) {
            i_nums = which(!is.na(values[cell_num,]));
            simple_values = na.omit(values[cell_num,]);
            if (length(simple_values) < min.lifetime) {
                next;
            }
            data$values = c(data$values,simple_values)
            data$file = c(data$file, rep(this_file,length(simple_values)));
            data$cell_num = c(data$cell_num, rep(cell_num,length(simple_values)));
            data$i_nums = c(data$i_nums, i_nums);
        }
        # print(paste('Done reading', this_file))
    }
    
    return(as.data.frame(data))
}

get_row_cumprods <- function(files,min.lifetime=0) {
    data = list();

    for (this_file in files) {
        values = as.matrix(read.csv(this_file,header=F));
        values_no_zero = values;
        values_no_zero[is.na(values)] = 0;
        # print(values_no_zero/100 + 1)
        for (cell_num in 1:dim(values)[1]) {
            total_data_points = sum(! is.na(values[cell_num,]))
            if (total_data_points < min.lifetime) {
                next;
            }
            # print(values_no_zero[cell_num,]/100 + 1);
            this_cumprod = cumprod(values_no_zero[cell_num,]/100+1);
            
            data$values = rbind(data$values,this_cumprod);
            data$file = c(data$file, this_file);
            data$cell_num = c(data$cell_num, cell_num);
            data$min_val = c(data$min_val, min(this_cumprod))
        }
    }
    
    return(data)
}

process_gel_diff_percents <- function(files,min.lifetime=0) {
    results = list()

    for (this_file in files) {
        percent_change = as.matrix(read.csv(this_file,header=F))
        for (cell_num in 1:dim(percent_change)[1]) {
            this_cell = percent_change[cell_num,]
            this_cell_na_omit = na.omit(this_cell)
            if (length(this_cell_na_omit) < min.lifetime) {
                next;
            }
            
            this_cumprod = 1-cumprod(this_cell_na_omit/100+1);
            this_cumprod_pad = c(this_cumprod, rep(NA, 51 - length(this_cumprod)));
            
            results$file = c(results$file, this_file)
            results$cell_num = c(results$cell_num, cell_num)
            results$start_time = c(results$start_time, which(!is.na(this_cell))[1])
            results$max_degrade = c(results$max_degrade,max(this_cumprod));
            results$min_degrade = c(results$min_degrade,min(this_cumprod));
                        
            results$seq$degrade_seq[[length(results$seq$degrade_seq) + 1]] = this_cumprod;
            results$seq$degrade_seq_pad[[length(results$seq$degrade_seq_pad) + 1]] = this_cumprod_pad;
        }
    }

    return(results)
}

colMedians <- function(this_mat,na.rm=T) {
    medis = c()
    for (i in 1:dim(this_mat)[2]) {
        medis = c(medis,median(this_mat[,i],na.rm=na.rm))
    }
    return(medis)
}

colConfUpper <- function(this_mat,conf.int=0.95) {
    upper = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            upper = c(upper,NA);
        } else {
            temp = t.test(this_mat[,i],conf.level=conf.int)
            if (is.nan(temp$conf.int[2])) {
                upper = c(upper,0)
            } else {
                upper = c(upper,temp$conf.int[2])
            }
        }
    }
    return(upper)
}

colConfLower <- function(this_mat,conf.int=0.95) {
    lower = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            lower = c(lower,NA);
        } else {
            temp = t.test(this_mat[,i],conf.level=conf.int)
            if (is.nan(temp$conf.int[1])) {
                lower = c(lower,0)
            } else {
                lower = c(lower,temp$conf.int[1])
            }
        }
    }
    return(lower)
}

colSD <- function(this_mat) {
    st_devs = c()
    for (i in 1:dim(this_mat)[2]) {
        if (all(is.na(this_mat[,i]))) {
            st_devs = c(st_devs,NA);
        } else {
            st_devs = c(st_devs,sd(this_mat[,i],na.rm=T))
        }
    }
    return(st_devs)
}

################################################################################
# Plotting
################################################################################

plot_single_exp_mat <- function(time,exp_mat,upper_conf=NA,lower_conf=NA,col=NA,...) {
    if (! is.na(col)) {
        col_light = rgb(t(col2rgb(col))/255,alpha=0.5);
    } else {
        col_light = NA;
    }
    
    for (i in 1:dim(exp_mat)[1]) {
        if (i == 1) {
            plot(time,exp_mat[i,],col=col_light,...)
        } else {
            lines(time,exp_mat[i,],col=col_light,...)
        }
    }
    
    lines(time,colMeans(exp_mat),col=col,lwd=3,...)
    
    if (! is.na(upper_conf[1]) && ! is.na(lower_conf[1])) {
        polygon(c(time,rev(time)), c(lower_conf,rev(upper_conf)),col=col_light,border=NA)
    }
}

plot_multiple_exp_mat <- function(time,exps,exp_confs,cols=NA,exp_names=NA,...) {
    #colors
    if (is.na(cols[1])) {
        cols = rainbow(length(names(exps)));
    }

    cols_light = list()
    for (exp_type in names(exps)) {
        cols_light[[exp_type]] = rgb(t(col2rgb(cols[[exp_type]]))/255,alpha=0.5)
    }

    #exps to include
    if (is.na(exp_names[1])) {
        exp_names = names(exps);
    }
    
    first = 1;
    for (exp_type in exp_names) {
        print(exp_type)
        if (first) {            
            plot(time,colMeans(exps[[exp_type]]),col=cols[[exp_type]],lwd=3,typ='l',...)
            first = 0;
        } else {
            lines(time,colMeans(exps[[exp_type]]),col=cols[[exp_type]],lwd=3,...)
        }
        polygon(c(time,rev(time)), c(exp_confs[[exp_type]]$upper,rev(exp_confs[[exp_type]]$lower)),
            col=cols_light[[exp_type]],border=cols_light[[exp_type]])
    }
}

################################################################################
# Data Loading
################################################################################

read_overall_results_files <- function(overall_dir) {
    files = Sys.glob(file.path(overall_dir,'*'))
    
    data = list()
    for (i in files) {
        if (regexpr('.csv',basename(i)) == -1) {
            next;
        }
        
        prop_name = sub('.csv','',basename(i));
        data[[prop_name]] = as.matrix(read.csv(i,header=F));
    }

    return(data)
}

find_degrade_percent <- function(overall_data,max.gel.minus=0,corrected.max.final=0) {
    gel.minus.degrade = !is.nan(overall_data$gel_minus_surrounding) & 
        overall_data$gel_minus_surrounding < max.gel.minus;
    
    corrected.final.degrade = !is.nan(overall_data$corrected_final_diff) & 
        overall_data$corrected_final_diff < corrected.max.final;
    corrected.final.degrade = corrected.final.degrade[,rep(1,dim(gel.minus.degrade)[[2]])]

    active_degrade = gel.minus.degrade & corrected.final.degrade;
    
    has_degraded = active_degrade;
    for (row_num in 1:dim(has_degraded)[[1]]) {
        degrade_times = which(has_degraded[row_num,]);
        if (length(degrade_times) == 0) {
            next;
        } else {
            has_degraded[row_num,degrade_times[1]:dim(has_degraded)[[2]]] = T
        }
    }
    
    degrade_percent = colSums(has_degraded)/dim(has_degraded)[[1]];
    
    return(degrade_percent)
}

find_average_degrade_percents <- function(data_sets,...) {
    these_percents = c();
    for (i in 1:length(data_sets)) {
        this_degrade_percent = find_degrade_percent(data_sets[[i]],...);
        these_percents = rbind(these_percents,this_degrade_percent)
    }

    return(colMeans(these_percents))
}

find_rmsd <- function(data_sets,...) {
    control_percents = find_average_degrade_percents(data_sets$control,...)
    BB94_percents = find_average_degrade_percents(data_sets$BB94,...)
    DMSO_percents = find_average_degrade_percents(data_sets$DMSO,...)
    FAK_percents = find_average_degrade_percents(data_sets$FAK_test,...)
    
    comp = list();
    comp$BB94_control_rmsd = sqrt(mean((control_percents - BB94_percents)**2));
    comp$BB94_zero_rmsd = sqrt(mean(BB94_percents**2));
    comp$BB94_final = tail(BB94_percents,1);
    comp$BB94_DMSO_final = tail(DMSO_percents,1) - tail(BB94_percents,1);
    comp$BB94_DMSO_rmsd = sqrt(mean((BB94_percents - DMSO_percents)**2));
    
    comp$DMSO_control_final = tail(control_percents,1) - tail(DMSO_percents,1);
    comp$DMSO_control_rmsd = sqrt(mean((control_percents - DMSO_percents)**2));
    
    comp$FAK_DMSO_rmsd = sqrt(mean((FAK_percents - DMSO_percents)**2));
    comp$FAK_DMSO_final = tail(FAK_percents,1) - tail(DMSO_percents,1);

    return(comp)
}
