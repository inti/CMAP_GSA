#!/usr/bin/env Rscript

### DEFINE SOME USEFUL FUNCTIONS

print_OUT<-function(string,LOG){
  print(paste(date(),"      ",string),quote = FALSE)  
}

page_z<-function(Sm,mu,m,rho ){  sqrt(m)*(Sm - mu)/rho  }

#########
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library(pROC,lib.loc="~/scratch/DATA/Rlibraries/"))

# library to use command line options
suppressPackageStartupMessages(library("optparse",lib.loc="~/scratch/DATA/Rlibraries/"))
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false", dest="verbose", help="Print little output"),
    make_option(c("-o","--outfile"), type="character", default="cmap_gsa.out.txt", help="Output file Name", metavar="NULL"),
    make_option(c("-d","--distance"), type="integer", default=20, help="Max SNP-to-gene distance allowed (in kb)", metavar=NULL),
    make_option(c("--gene_file"), type="character", help="File with genes to be analyse", metavar=NULL),
    make_option(c("--clinical_trials"), default= "data/trials_4_cmap_drugs.txt",type="character", help="File with gene conditions", metavar=NULL),
    #make_option(c("--pc_var"), default= 0.05,type="double", help="Threshold of variance explained to include PC in analysis", metavar=NULL),
    make_option(c("--auc"), action="store_true", default= FALSE, help="Calcukate Area Under the Curve", metavar=NULL),	
    make_option(c("-b","--bootstraping"), action="store_true", default= "FALSE", help="Set AUC calculation on bootstraping mode. Default: DeLong method.", metavar=NULL),
    make_option(c("--correct_ci"), action="store_true", default= 'FALSE', help="Set correction for multuple testing on confidence interval of AUC values", metavar=NULL),
    make_option(c("--n.boot"), default= 1000,type="integer", help="Number of bootstraping replicates to calculate AUC values", metavar=NULL),
    make_option(c("--ci_boot"), default= 0.95,type="double", help="Condidence interval value for AUC values", metavar=NULL),
    make_option(c("--plot"), action="store_true",default= FALSE,type="logical", help="Make plots of results. It need library ggplot2 installed", metavar=NULL),
    make_option(c("--gsa"), action="store_true",default= FALSE,type="logical", help="Run enrichemnt analysis using PAGE method", metavar=NULL),
    make_option(c("--n_perm"), default= 0,type="integer", help="Number of permutation to calculate the null distribution and significance for  enrichmet results", metavar=NULL),
    make_option(c("--one_sided"), action="store_true",default= TRUE,type="logical", help="Enrichment test is one sided. Differential regulation regarding is up or down", metavar=NULL),
    make_option(c("--two_sided"), action="store_true",default= FALSE, dest="one_sided",type="logical", help="Enrichment test is two sided. Test both up and down regulation separately", metavar=NULL),    
    make_option(c("--lfdr"), action="store_true",default= FALSE,type="logical", help="Calculate local FDR values for PAGE results", metavar=NULL),
    make_option(c("--hugo_id"), action="store_true",default= TRUE,type="logical", help="User gene ids are Hugo symbols", metavar=NULL),
    make_option(c("--ensembl_id"), action="store_true",dest="hugo_id",default= FALSE,type="logical", help="User gene ids are ensembl gene ids", metavar=NULL),
    make_option(c("--collapse_probes_mean"), action="store_true",default= FALSE,type="logical", help="Collapse probe level information into a single gene-leve value by averaging individual probe values", metavar=NULL),
    make_option(c("--collapse_probes_abs_max"), action="store_true",default= FALSE,type="logical", help="Choose the row with the highest absolute value", metavar=NULL),	
    make_option(c("--collapse_probes_abs_min"), action="store_true",default= FALSE,type="logical", help="Choose the row with the lowest absolute value", metavar=NULL),
    make_option(c("--collapse_probes_eigen"), action="store_true",default= FALSE,type="logical", help="Choose the row with the lowest absolute value", metavar=NULL)
)

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)

# DESCRIPTION
# received as input a file with genes and conditions test is the expression of the genes is affected in the CMAP data.
###########

# Check user has provided all necesary options
if (length(opt$gene_file) == 0){
	print_OUT("Need to provide a files with gene ids and conditions.");
	print_OUT("Finishing execution now.");
	quit();
}
if (length(opt$lfdr) != 0) {
	suppressPackageStartupMessages(library(fdrtool))
}
if (length(opt$plot) != 0) {
        suppressPackageStartupMessages(library(ggplot2))
}


##### NOW PROCEED WITH THE ANALYSIS

# load CMAP data
print_OUT("Loading CMAP data from [ data/CMAP.RData ]. It may take a few minutes.")
load("data/CMAP.RData") 

# read probe annotation 
print_OUT("Loading AFFY probes to gene mapping to parse CMAP data from [ data/probes_annotation.txt ]")
probe_annot<-read.table("data/probes_annotation.txt",header=T,sep="\t")
#print_OUT("Detecting probes that match to more than one gene")
#probe_annot = ddply(probe_annot, .(affy_hg_u133_plus_2), function(d) if (length(unique(d$ensembl_gene_id)) == 1){return(d)},.progress=create_progress_bar(name="text"),.parallel=T)

print_OUT("Filtering out probes matching to multiple genes")
avg_data = avg_data[,unique(probe_annot$affy_hg_u133_plus_2)]

# read drug trials information
print_OUT(paste("Reading conditions for genes from [ ",opt$clinical_trials," ]",sep=""))
trials<-read.table( opt$clinical_trials ,header=T,sep="\t")
trials_phase_4= unique(trials[which(trials$phase == "phase_iv"),])
# filter diseases
min_drugs = 5
trials_phase_4=ddply(trials_phase_4, .(disease), function(data) {
		if (length(unique(data$drug)) >= min_drugs ){
			return(data)
		}
	}
)

# read input file
# two tab separated columns expected: condition_name, gene_id
print_OUT(paste("Reading gene-sets of to be analysed from [ ",opt$gene_file," ]",sep=""))
data=read.table(file=opt$gene_file,sep="\t")
colnames(data)=c("condition","gene_id")
data$gene_id = toupper(data$gene_id)
# add annotation to data.
if (opt$hugo_id == FALSE){
	annotated_data<-merge(data,probe_annot[,c("affy_hg_u133_plus_2","ensembl_gene_id","hgnc_symbol")],by.x="gene_id",by.y="ensembl_gene_id")
} else {
	annotated_data<-merge(data,probe_annot[,c("affy_hg_u133_plus_2","ensembl_gene_id","hgnc_symbol")],by.x="gene_id",by.y="hgnc_symbol")
}

##### Define CMAP DATA

collapse_method = 0;
if ( opt$collapse_probes_mean == TRUE ){
	# summarise the log-fold-change values at the gene level
	collapse_method = "Average"; # default avegare across genes
}
if (opt$collapse_probes_abs_max == TRUE) {
	collapse_method = "absMaxMean"; # max abs value
}
if (opt$collapse_probes_abs_min == TRUE ){
	collapse_method = "absMinMean"; # min absolute value
}
if (opt$collapse_probes_eigen == TRUE ){
        collapse_method = "ME"; # eigen vector
}
if (collapse_method != 0){
	print_OUT(paste("Collapsing probe values at gene level using method [ ",collapse_method," ].",sep=""))
	suppressPackageStartupMessages(library(WGCNA,lib.loc="~/scratch/DATA/Rlibraries/"))
	gene_id = ddply(probe_annot, .(affy_hg_u133_plus_2), summarise, gene_id=unique(ensembl_gene_id),.parallel=T)
	data = t(avg_data)
	data = collapseRows( data,rowGroup = gene_id$gene_id,rowID = rownames(data),method=collapse_method)
	avg_data = t(data$datETcollapsed)
	tmp_cols=data.frame(colnames(avg_data))
	annotated_data2=merge(annotated_data,tmp_cols,by.x="ensembl_gene_id",by.y="colnames.avg_data.")
}


##### PAGE analysis #####
print_OUT("Running PAGE analysis");
print_OUT(paste("   '-> WIll run [ ",opt$n_perm," ] permutations to calculate the statistic under the null",sep=""))
drug_summary_stats = ldply(rownames(avg_data), function(drug){ 
		data = avg_data[drug,];
		data_pow2 = data^2;
		return(c(mean(data,na.rm=T),sd(data,na.rm=T),mean(data_pow2,na.rm=T),sd(data_pow2,na.rm=T)  ));
	}
)
rownames(drug_summary_stats) = rownames(avg_data);
colnames(drug_summary_stats) = c("mean_logFC","sd_logFC","mean_logFC_pow2","sd_logFC_pow2")


drug_de_df=ddply(annotated_data, .(condition), function(d) {
		d=unique(d);
		if (collapse_method != "NULL"){
			cmap_subset=avg_data[,d$ensembl_gene_id];
		} else {
			cmap_subset=avg_data[,d$affy_hg_u133_plus_2];
		}
		if (opt$one_sided == "TRUE") {
			
			observed_page_z = ldply(rownames(cmap_subset), function(drug)  {
					lfc = cmap_subset[drug,];
					Sm = mean(lfc,na.rm=T);
					m = sum(length(lfc) - abs(is.na(lfc)));
					mu = drug_summary_stats[drug,"mean_logFC"];
					rho = drug_summary_stats[drug,"sd_logFC"];
					z = page_z( Sm ,mu,m,rho);
					return(z); 
				}
			,.parallel=T)
		} else { # else is two sided test
			observed_page_z = ldply(rownames(cmap_subset), function(drug)  {  
					lfc = cmap_subset[drug,];
					lfc = lfc^2;
					Sm = mean(lfc,na.rm=T);
					m = sum(length(lfc) - abs(is.na(lfc)));
					mu = drug_summary_stats[drug,"mean_logFC_pow2"];
					rho = drug_summary_stats[drug,"sd_logFC_pow2"];
					z = page_z( Sm ,mu,m,rho);
					return(z); 
				}
			,.parallel=T)
		}
		observed_page_z$drug = rownames(cmap_subset);
		colnames(observed_page_z) = c("page_z","drug");

		if (opt$n_perm > 0) { 
			geneset_size = ncol(cmap_subset); 

			background = ldply(1:opt$n_perm, function(x)   {
					random_index = sample(1:ncol(avg_data), geneset_size);
					random_data = avg_data[,random_index];
					if (opt$one_sided == "TRUE") {
						random_page_z = ldply(rownames(random_data), function(drug)  {
								lfc = random_data[drug,];
								Sm = mean(lfc,na.rm=T);
								m = sum(length(lfc) - abs(is.na(lfc)));
								mu = drug_summary_stats[drug,"mean_logFC"];
								rho = drug_summary_stats[drug,"sd_logFC"];
								z = page_z( Sm ,mu,m,rho);
								return(z);
							}
						,.parallel=T)
					} else {
						random_page_z = ldply(rownames(random_data), function(drug)  {
								lfc = random_data[drug,];
								lfc = lfc^2;
								Sm = mean(lfc,na.rm=T);
								m = sum(length(lfc) - abs(is.na(lfc)));
								mu = drug_summary_stats[drug,"mean_logFC_pow2"];
								rho = drug_summary_stats[drug,"sd_logFC_pow2"];
								z = page_z( Sm ,mu,m,rho);
								return(z);
							}
						,.parallel=T)

					}
					return(t(random_page_z));
				}     
			,.parallel=T)
			background = t(background)
			rownames(background) = rownames(cmap_subset);
			drug_z = ldply(observed_page_z$drug, function(drug)  {
					background_mean = mean(background[drug,],na.rm=T)
					background_sd = sd(background[drug,],na.rm=T)
					z = (observed_page_z[which(observed_page_z$drug == drug),"page_z"] - background_mean)/background_sd; 
					return( c(drug,z,background_mean,background_sd )  );
				} 
			)
			colnames(drug_z) = c("drug","empirical_Z","background_mean","background_sd");
			drug_z$empirical_Z = as.numeric(drug_z$empirical_Z);
			drug_z$background_mean = as.numeric(drug_z$background_mean);
			drug_z$background_sd = as.numeric(drug_z$background_sd);
			drug_z = drug_z[order(-drug_z$empirical_Z),];
			if (opt$one_sided == TRUE){ 
				drug_z$lfdr_empirical_z<-fdrtool( pnorm(drug_z$empirical_Z,lower.tail=F),statistic = "pvalue", plot=F,verbose=F )$lfdr
			} else {
				drug_z$lfdr_empirical_z<-fdrtool( drug_z$empirical_Z,statistic = "normal", plot=F,verbose=F )$lfdr
			}
			observed_page_z = merge(observed_page_z,drug_z,by="drug")
		} 
		## add FDR ###
		if (opt$one_sided == TRUE){ 
			observed_page_z$lfdr_page_z<-fdrtool( pnorm(observed_page_z$page_z,lower.tail=F),statistic = "pvalue", plot=F,verbose=F )$lfdr
		} else {
			observed_page_z$lfdr_page_z<-fdrtool( observed_page_z$page_z,statistic = "normal", plot=F,verbose=F )$lfdr
		}
		return(observed_page_z);
	}
,.progress=create_progress_bar(name="text"),.parallel=F)

if (opt$gsa == "TRUE"){
	gsa_out_file = paste(opt$outfile,".PAGE.txt",sep="")
	print_OUT(paste("   '-> Writing GSA results to [ ", gsa_out_file ," ]",sep=""));
	write.table(drug_de_df,file=gsa_out_file,quote=FALSE,row.names=F,col.names=T,sep="\t");
}



##### AUC analysis #####
# For condition get the gene data across all drugs instances and calculate the eigen values.
# return a list on which each element is a matrix with the PC for each condition.
if (opt$auc == TRUE ) {
	print_OUT("Calculating AUC values for each gene-set");
	#print_OUT(paste("   '-> Selecting for analysis PC explaining more than [ ", opt$pc_var ," ] of variation.",sep=""));
	drug_pca_df=dlply(drug_de_df, .(condition),  function(d) {
			d=unique(d);
                        if (collapse_method != "NULL"){
                                cmap_subset=avg_data[,d$ensembl_gene_id];
                        } else {
                                cmap_subset=avg_data[,d$affy_hg_u133_plus_2];
                        }
			rot=merge(d,trials_phase_4,by="drug");
			return(rot);
                 }
	,.progress=create_progress_bar(name="text"),.parallel=F)

	if (opt$n_perm > 0) {
		print_OUT("   '-> For each condition I will calculate AUC values based on the [ empirical_Z ] column")
		col_val = "empirical_Z";
	} else {
		col_val = "page_z";
	}
	auc_table = llply(names(drug_pca_df), function(condition) {
			rot = drug_pca_df[[condition]];
			target = llply(unique(rot$disease), function(condition)   abs(rot$disease == condition)    );
			names(target) = unique(rot$disease);
			n_test = length(names(target));
			if (opt$correct_ci == TRUE){
				conf.level = 1 - (1 - opt$ci_boot)^n_test; 
			 } else {
				conf.level=opt$ci_boot;
			}
			ldply( names(target), function(condition) {
					if (opt$boot == "FALSE" & sum(response = target[[ condition  ]],na.rm=T) == 1) { return() }                                
					if (opt$boot == "TRUE"){
					   auc_values<-ci.auc(response = target[[ condition  ]], predictor = rot[,col_val],
						boot.stratified=TRUE,boot.n=opt$n.boot,direction="<", conf.level =conf.level,
						method="bootstrap",progress=create_progress_bar(name="none"));
					} else {
					   # Using DeLong method
					   auc_values<-ci.auc(response = target[[ condition  ]], predictor = rot[,col_val],direction="<",conf.level =conf.level,progress=create_progress_bar(name="none"));
					}
					auc_values=t(as.data.frame(unclass(auc_values)));
					auc_values=as.data.frame(auc_values);
					colnames(auc_values)<-c("low_ci","median","up_ci");
					auc_values$disease<-condition;
					auc_values$N_drugs<-sum(target[[condition]],na.rm=T);
					#auc_values$eigen_vector = eigen_vector;
					auc_values$conf.level = conf.level;
					return(auc_values);
				}
			)
		}
	,.progress=create_progress_bar(name="text"))
	names(auc_table)=names(drug_pca_df);

	# PLOT AND WRITTE TO FILES
	if (opt$plot == TRUE){
		print_OUT(paste("Writting output files and plots to [ ",opt$outfile,".condition.* ] with extensions txt and pdf, respectively ",sep=""));
	} else {
		print_OUT(paste("Writting output files to [ ",opt$outfile,".condition.txt ]",sep=""));
	}

	l_ply(names(auc_table), function(condition) {
			data=auc_table[[ condition  ]];
			data$touch_05 = abs( (data$low_ci - 0.5)*( data$up_ci - 0.5 )  > 0 );
			data_filtered = data[which(data$touch_05 > 0),];
			data_filtered$disease<-factor(data_filtered$disease, levels=  unique(factor(data_filtered[order(data_filtered$low_ci),"disease"])));
			data$condition = condition
			#my_manual_col=c("grey","black")
			#if (length(unique( data$eigen_vector)) == 1){
			#	my_manual_col = "black";
			#}
			dodge = position_dodge(width = 0.9, height = NULL)
			if (opt$plot == TRUE){
				#data_filtered[,"disease"]=as.factor(as.character(data_filtered[,"disease"]));
				plot=ggplot(data=data_filtered, aes(x=disease,y=median)) + 
					geom_point(col="black", position=dodge) +
					geom_linerange(aes(ymin=low_ci,ymax=up_ci),col="grey", position=dodge) + 
					#coord_flip() + 
					#scale_colour_manual(values = my_manual_col) + 
					opts(   axis.text.x=theme_text(size=10, angle = 35, vjust=1,hjust=1),
						panel.background=theme_rect(colour=NA),
						panel.grid.minor=theme_blank(),
						panel.grid.major=theme_blank(),
						axis.line=theme_segment(),
						legend.position="none");
				file_plot=paste(opt$outfile,".",condition,".AUC.pdf",sep="");
				ggsave(plot=plot,file=file_plot,width=15,height=7);
			}
			file_table=paste(opt$outfile,".",condition,".AUC.txt",sep="");
			write.table(data,file=file_table,quote=FALSE,row.names=F,col.names=T,sep="\t");
		}
	,.progress=create_progress_bar(name="text"))
}

print_OUT("Analysis finished.")
