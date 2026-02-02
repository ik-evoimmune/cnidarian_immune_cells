setwd("/sci/labs/yehum79/itamar273/immune_cells_RNAseq_New_Aanalysis/Arnau_metacell_analysis/scRNAseq_analysis")

library(metacell)
source("Accessory_functions.R")
library(zoo)
library(scales)

##########Calculating metacells
scdb_init("mc_db",force_reinit=T)
mcell_import_multi_scmat_10x(mat_nm="Nvec",dataset_table_fn="10x_import_table.txt",force=F,base_dir=".")

mat=scdb_mat("Nvec")

cz=Matrix::colSums(mat@mat)

#load clicktag classification and remove non-classified cells
ct_info=read.table("clicktag_info.txt",h=T,row.names=1,sep="\t",stringsAsFactors=F)
good_cells=rownames(ct_info[which(ct_info$classification=="cell"),])
bad_cells=setdiff(colnames(mat@mat),good_cells)

mcell_mat_ignore_cells("Nvec_f","Nvec",ig_cells=bad_cells)
mat=scdb_mat("Nvec_f")

#remove low-coverage orphan peak clusters (see https://github.com/sebepedroslab/GeneExt for details)
peak_counts=rowSums(as.matrix(mat@mat[grepl("orphan",rownames(mat@mat)),]))
peaks_kill=names(which(peak_counts < 50))
mcell_mat_ignore_genes("Nvec_f","Nvec_f",ig_genes=peaks_kill)
mat=scdb_mat("Nvec_f")


#calculate gene-level statistics and select genes for clusteirng
mcell_add_gene_stat(gstat_id="gstat", mat_id="Nvec_f", force=T)
mcell_gset_filter_multi("gstat","clust_markers",T_tot=30,T_top3=2,T_szcor=-0.1,T_niche=0.1,force_new=T,bl=c("mCherry_plus_strand","mCherry_minus_strand"))

#calculate metacell solution
mcell_add_cgraph_from_mat_bknn(mat_id="Nvec_f",gset_id="clust_markers",graph_id="graphk100",K=100,dsamp=F)
mcell_coclust_from_graph_resamp(coc_id="coc500_min20",graph_id="graphk100",min_mc_size=8,p_resamp=0.75,n_resamp=500)
mcell_mc_from_coclust_balanced("coc500_min20",mat="Nvec_f",mc_id="mc_k30",K=30,min_mc_size=10,alpha=2)


##########Reordering metacells
mc=scdb_mc("mc_k30")
#compute normalized confusion matrix, max degree? ~100
confu_norm=mc_compute_norm_confu_matrix("mc_k30","graphk100",100)
hc=mc_confusion_clustering(confu_norm)

#visualize tree and decide on a cut height for color assignment
plot(hc)
rect.hclust(hc,h=5.5)
mc_clusts=cutree(hc,h=5.5)

tmp_new_colors=colorRampPalette(c("dodgerblue3", "chocolate3","chartreuse3","brown3", "slateblue3","plum3","darkgoldenrod3", "aquamarine3","coral3"))(length(table(mc_clusts[hc$order])))
#reorder the generated colors
tmp_new_colors=as.character(treemap::treepalette(as.data.frame(tmp_new_colors))$tmp_new_colors)
new_colors=rep(tmp_new_colors, table(factor(mc_clusts[hc$order],levels=unique(mc_clusts[hc$order]))))  ##this allows us to control gradient of colors based on reordering.

mc_reord=mc_reorder(mc,hc$order)
mc_reord@colors=new_colors
names(mc_reord@colors)=colnames(mc_reord@mc_fp)
scdb_add_mc("mc_reord",mc_reord)



##########Gene expression maps
mc=scdb_mc("mc_reord")
mat=scdb_mat("Nvec_f")
scr_plot_cmod_markers(mc,mat,clust_ord=1:ncol(mc@mc_fp),output_file="Nvec_gene_expression.png",gene_annot_file="Nvec_v4_annotation.txt",per_clust_genes=50,height=12000,width=4000,gene_min_fold=2,plot_sc=F,plot_mc_names=TRUE,plot_v_lines=TRUE,add_lateral_genes=T,include_peaks=F,clust_col=mc@colors,pmin=5)




#########Experimental conditions distribution
ct_info=read.table("clicktag_info.txt",h=T,row.names=1,sep="\t",stringsAsFactors=F)
cell_class=as.character(ct_info[colnames(mat@mat),"clicktag_label"])
names(cell_class)=colnames(mat@mat)
niche_order = 1:ncol(mc@mc_fp) 

x=sapply(1:ncol(mc@mc_fp),function(x) table(factor(cell_class[names(which(mc@mc==x))],levels=c("Ctrl","iHCl","tPIC"))))
colnames(x)=as.character(niche_order)
categories=c("Ctrl","iHCl","tPIC")
color=RColorBrewer::brewer.pal(3,"Set1")
names(color)=categories

pdf("Experiment_distribution.pdf",height=8,width=25,useDingbats=F)
barplot(x,las=2,col=color,cex.names=1,cex.axis=2)
legend("topleft",legend=names(color),fill=color,box.lty=0,horiz=F,ncol=round(length(color)/3,0),cex=1.5)
dev.off()


##########Comparison with reference gastrula scRNA-seq atlas
mc_ref=scdb_mc("gastrula_K30")
annot_ref=read.table("mc_db/annotation.gastrula_K30.tsv",h=T,sep="\t",quote="",stringsAsFactors=F,comment.char="")

var_genes=intersect(names(which(apply(mc@mc_fp,1,max)>2)),names(which(apply(mc_ref@mc_fp,1,max)>2)))
x=scr_clustering_comparison(mc_ref@mc_fp,mc@mc_fp,markers=var_genes,out_fn="Reference_vs_new.png",original_ordering_1st=T,original_ordering_2nd=T,h=3000,w=3000)



