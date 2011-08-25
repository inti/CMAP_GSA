#!/usr/bin/env Rscript

### DEFINE SOME USEFUL FUNCTIONS

print_OUT<-function(string,LOG){
  print(paste(date(),"      ",string),quote = FALSE)  
}


#########
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pROC,lib.loc="~/scratch/DATA/Rlibraries/"))

# library to use command line options
suppressPackageStartupMessages(library("optparse",lib.loc="~/scratch/DATA/Rlibraries/"))
option_list <- list(
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quietly"), action="store_false", dest="verbose", help="Print little output"),
    make_option(c("-o","--outfile"), type="character", default="forge.out.txt", help="Output file Name", metavar=NULL),
    make_option(c("-d","--distance"), type="integer", default=20, help="Max SNP-to-gene distance allowed (in kb)", metavar=NULL),
    make_option(c("--gene_file"), type="character", help="File with genes to be analyse", metavar=NULL),
    make_option(c("--clinical_trials"), default= "data/trials_4_cmap_drugs.txt",type="character", help="File with gene conditions", metavar=NULL),
    make_option(c("--pc_var"), default= 0.05,type="double", help="Threshold of variance explained to include PC in analysis", metavar=NULL),
    make_option(c("--auc"), action="store_true", default= FALSE, help="Calcukate Area Under the Curve", metavar=NULL),	
    make_option(c("-b","--bootstraping"), action="store_true", default= "FALSE", help="Set AUC calculation on bootstraping mode. Default: DeLong method.", metavar=NULL),
    make_option(c("--correct_ci"), action="store_true", default= 'FALSE', help="Set correction for multuple testing on confidence interval of AUC values", metavar=NULL),
    make_option(c("--n.boot"), default= 1000,type="integer", help="Number of bootstraping replicates to calculate AUC values", metavar=NULL),
    make_option(c("--ci_boot"), default= 0.95,type="double", help="Condidence interval value for AUC values", metavar=NULL),
    make_option(c("--plot"), action="store_true",default= FALSE,type="logical", help="Make plots of results. It need library ggplot2 installed", metavar=NULL)	
)

# get command line options, if help option encountered print help and exit, 
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

# DESCRIPTION
# received as input a file with genes and conditions test is the expression of the genes is affected in the CMAP data.
###########

# load CMAP data
print_OUT("Loading CMAP data from [ data/CMAP.RData ]")
load("data/CMAP.RData") 

# read probe annotation 
print_OUT("Loading AFFY probes to gene mapping to parse CMAP data from [ data/probes_annotation.txt ]")
probe_annot<-read.table("data/probes_annotation.txt",header=T,sep="\t")

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

# add annotation to data.
annotated_data<-merge(data,probe_annot[,c("affy_hg_u133_plus_2","ensembl_gene_id","hgnc_symbol")],by.x="gene_id",by.y="hgnc_symbol")

# For condition get the gene data across all drugs instances and calculate the eigen values.
# return a list on which each element is a matrix with the PC for each condition.

print_OUT("Calculating PCA for each gene-set");
print_OUT(paste("   '-> Selecting for analysis PC explaining more than [ ", opt$pc_var ," ] of variation.",sep=""));
drug_pca_df=dlply(annotated_data, .(condition), function(d) {
		d=unique(d);
		cmap_subset=avg_data[,d$affy_hg_u133_plus_2];
		drug_pca=prcomp(t(cmap_subset),center=T,scale=T);
		matrix=summary(drug_pca)[["importance"]]
		good_pcs<-which(matrix["Proportion of Variance",] >= opt$pc_var);
		rot = as.data.frame(drug_pca$rotation);
		rot = as.data.frame(rot[,names(good_pcs)]);
		colnames(rot)=names(good_pcs)	
		rot$drug=rownames(avg_data);
		rot=merge(rot,trials_phase_4,by="drug");	
		return(rot);
	} 
,.progress=create_progress_bar(name="text"),.parallel=F)

eigen_vector_auc_table = llply(names(drug_pca_df), function(condition) {
		rot = drug_pca_df[[condition]];
                target = llply(unique(rot$disease), function(condition)   abs(rot$disease == condition)    );
                names(target) = unique(rot$disease);
		pcs = unique(grep("PC",colnames(rot),value=T));
		N_pcs = length(pcs);
		n_test = N_pcs*length(names(target));
		if (opt$correct_ci == TRUE){
			conf.level = 1 - (1 - opt$ci_boot)^n_test; 
		 } else {
			conf.level=opt$ci_boot;
		}	
		table = ldply( pcs,  function(eigen_vector) {
                	#cat("\tWorking on eigen vector [ ",eigen_vector," ]\n",sep="")
                	ldply( names(target), function(condition) {
                               		if (opt$boot == "FALSE" & sum(response = target[[ condition  ]],na.rm=T) == 1) { return() }                                
					if (opt$boot == "TRUE"){
                               		   auc_values<-ci.auc(response = target[[ condition  ]], predictor = rot[,eigen_vector],
                               		        boot.stratified=TRUE,boot.n=opt$n.boot,direction="<", conf.level =conf.level,
                               		        method="bootstrap",progress=create_progress_bar(name="none"));
                               		} else {
                               		   # Using DeLong method
                               		   auc_values<-ci.auc(response = target[[ condition  ]], predictor = rot[,eigen_vector],direction="<",conf.level =conf.level,progress=create_progress_bar(name="none"));
                               		}
                               		auc_values=t(as.data.frame(unclass(auc_values)));
                               		auc_values=as.data.frame(auc_values);
                               		colnames(auc_values)<-c("low_ci","median","up_ci");
                               		auc_values$disease<-condition;
                               		auc_values$N_drugs<-sum(target[[condition]],na.rm=T);
                               		auc_values$eigen_vector = eigen_vector;
					auc_values$conf.level = conf.level;
                               		return(auc_values);
                        	}
			)
		} ,.parallel=T);
		return(table);
	}
,.progress=create_progress_bar(name="text"))
names(eigen_vector_auc_table)=names(drug_pca_df);

# PLOT TO FILES

l_ply(names(eigen_vector_auc_table), function(condition) {
		data=eigen_vector_auc_table[[ condition  ]];
		data$touch_05 = abs( (data$low_ci - 0.5)*( data$up_ci - 0.5 )  > 0 );
		data_filtered = data[which(data$touch_05 > 0),];
		#data_filtered = ddply(data, .(disease), function(d)  if(sum(d$touch_05) > 0) { return(d) });
		#data_filtered$disease<-factor(data_filtered$disease, levels=  unique(factor(data_filtered$disease)));
		data_filtered$disease<-factor(data_filtered$disease, levels=  unique(factor(data_filtered[order(data_filtered$low_ci),"disease"])));
		data$condition = condition
		#my_manual_col=c("grey","black")
		#if (length(unique( data$eigen_vector)) == 1){
		#	my_manual_col = "black";
		#}
		dodge = position_dodge(width = 0.9, height = NULL)
		if (opt$plot == TRUE){
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
			file_plot=paste("test_AUCs.",condition,".pdf",sep="");
			ggsave(plot=plot,file=file_plot,width=15,height=7);
		}
		file_table=paste("test_AUCs.",condition,".txt",sep="");
		write.table(data,file=file_table,quote=FALSE,row.names=F,col.names=T,sep="\t");
	}
,.progress=create_progress_bar(name="text"))



