mc_compute_norm_confu_matrix=function(mc_id,graph_id,max_deg=NULL){
	cgraph = scdb_cgraph(graph_id)
	if(is.null(max_deg)){
		max_deg = nrow(cgraph@edges)
	}
	confu = mcell_mc_confusion_mat(mc_id, graph_id, max_deg,ignore_mismatch=T)
	r_confu = rowSums(confu)
    c_confu = colSums(confu)
    norm = r_confu %*% t(c_confu)
    confu_n = confu/norm
    confu_nodiag = confu_n
    diag(confu_nodiag) = 0
    confu_n = pmin(confu_n, max(confu_nodiag))
    confu_n = pmin(confu_n, quantile(confu_n, 1 - 3/nrow(confu_n)))
    
	return(confu_n)
}

mc_confusion_clustering=function(confu_n,clust_method="average"){
    epsilon = quantile(confu_n[confu_n != 0], 0.02)
    hc = hclust(as.dist(-log10(epsilon + confu_n)), clust_method)
	return(hc)
}

scr_plot_cmod_markers=function(mc_object,mat_object,black_list=c(),output_file,clust_ord = c(),height=8000, width=3000,sub_list_mc=NULL,plot_sc=F, smoothen=5,
                                                clust_col=NULL,per_clust_genes=20,gene_min_fold=3,gene_annot_file,gene_list=NULL,plot_mc_names=FALSE,plot_v_lines=FALSE,
                                                transverality_N=ncol(mc_object@mc_fp),transv_excluded_niches=NULL,order_genes=T,pmin_v=5,add_lateral_genes=FALSE,include_peaks=FALSE)
  {    
    
        annot=read.table(gene_annot_file,header=T,sep="\t",fill=TRUE,quote="",row.names=1)
    #mat_niche = g_fp[scr_markers,];colnames(mat_niche)=seq(1:ncol(mat_niche))
        if(is.null(sub_list_mc)){
        niche_geomean_n= mc_object@mc_fp
        }else{
             niche_geomean_n= mc_object@mc_fp[,sub_list_mc]
             clust_ord=sub_list_mc
        }

        if(is.null(gene_list)){
        	genes=unique(as.vector(unlist(apply(niche_geomean_n, 2, function(x) names(head(sort(-x[x>gene_min_fold]),n=per_clust_genes))))))
        	transversal_genes=names(which(apply(niche_geomean_n[,setdiff(as.character(colnames(niche_geomean_n)),transv_excluded_niches)], 1, function(x) sort(x,decreasing=T)[transverality_N]>1.4)))
        	genes=setdiff(genes, transversal_genes)
        }else{
        	genes=gene_list
        }
		
		if(add_lateral_genes){
			lateral_genes=names(which(apply(niche_geomean_n, 1, function(x) sort(x,decreasing=T)[round(ncol(niche_geomean_n)/10)]>1.5)))  ##genes overexpressed in 1/10th of the metacells
			genes=union(genes,lateral_genes)
		}
		if(include_peaks==FALSE){
			genes=genes[!grepl("peak_",genes)]
		}
		genes=intersect(genes,rownames(niche_geomean_n))
		genes=setdiff(genes, black_list)
		
		message("Will use ",length(genes)," genes")

        mat_niche=niche_geomean_n[genes,]
                  
    if(length(clust_ord)==0) {
			message("recomputing cell ord")

			hc1 = hclust(dist(cor(mat_niche,method="pearson")), "ward.D2")
			clust_ord = as.character(hc1$order)
			write.table(clust_ord,file="tmp_cell_clusts_ordered_by_scr_markers_plot.txt",quote=FALSE,col.names = FALSE,row.names=FALSE)
			png(paste0(output_file,"_TREE.png"),h=500,w=1000)
			plot(hc1,xlab="",xaxt='n',hang=-1,ylab="",main="",cex=1)
			dev.off()  
			scr_tmp_niche_order<<-as.character(hc1$order)
		}

    #hc2 = hclust(dist(cor(t(mat_niche), m="spearman")), "ward.D2")
    #hc2 = hclust(dist(cor(t(lus))), "ward.D2")
    #hc2$order = order(apply(mat_niche[,as.character(clust_ord)],1,function(x) which.max(rollmean(x,5))))
    if(order_genes){
                gene_ord=order(apply(mat_niche[,as.character(clust_ord)],1,function(x) which.max(rollmean(x,1))))
        }else{
                gene_ord=1:length(genes)
        }
        #hc_gmods = cutree(hc2, 50)
        #mean_e = as.matrix(apply(t(lus_cl_fp), 1, function(x) tapply(x, hc_gmods, mean)))
        #gene_ord = order(apply(mean_e[hc_gmods,clust_ord],1,function(x) which.max(rollmean(x,3))))
    write.table(genes[gene_ord],file="tmp_markers_ordered.txt",quote=FALSE,col.names = FALSE,row.names=FALSE)
    png(paste0(output_file,".png"), h=height, w=width)
    par(mar=c(0,0,0,0))
    par(fig=c(0.25,0.75,0.1,0.9))
    #mat_niche = g_fp[scr_markers,];colnames(mat_niche)=seq(1:ncol(mat_niche))
    #mat_niche_to_plot = pmin(pmax(mat_niche[hc2$order,], 0), 30)
    
    shades=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
    #image(t(pmax(log2(mat_niche_to_plot[, clust_ord]),0)), col=shades,xaxt="n",yaxt="n")
    image(t(pmin(log2(niche_geomean_n[genes[gene_ord],as.character(clust_ord)]+1),pmin_v)), col=shades,xaxt="n",yaxt="n")
        mtext(annot[genes[gene_ord],2], side=4, 
          at=seq(0,1,length.out=length(genes[gene_ord])), las=1, line=1)
    mtext(paste(annot[genes[gene_ord],1],genes[gene_ord],sep="||"), side=2, 
          at=seq(0,1,length.out=length(genes[gene_ord])), las=1, line=1)
    
    if(plot_v_lines){
    	n_mcs_per_color=table(clust_col)[unique(clust_col)]
    	x=0
    	for(i in 1:length(n_mcs_per_color)) {
			abline(v=(n_mcs_per_color[i]+x)/(sum(n_mcs_per_color)-1)-1/(2*sum(n_mcs_per_color)), lwd=2,col=alpha("black",0.4))      
			x=x+n_mcs_per_color[i]
    	}
    }
    
    if(plot_mc_names==TRUE){
    	mtext(clust_ord,side=1,at=seq(0,1,length.out=sum(ncol(mat_niche))), las=2, line=1,cex=2)
    	mtext(clust_ord,side=3,at=seq(0,1,length.out=sum(ncol(mat_niche))), las=2, line=1,cex=2)
    }
      
	if(!is.null(clust_col)){
		par(fig=c(0.25,0.75,0.05,0.08),new=TRUE)
		if(is.null(names(clust_col))){names(clust_col)=as.character(clust_ord)}
        image(as.matrix(1:length(clust_ord)),col=clust_col[as.character(clust_ord)], axes = F,xaxt='n',yaxt='n')
        }
    dev.off()

        print(length(genes))

        ##################PLOT SINGLE-CELL PROFILE########################
    if(plot_sc){
	cell_order=c()  
    for (niche in clust_ord){
			cells=names(mc_object@mc[which(mc_object@mc==niche)])    
          cell_order=c(cell_order,cells)
    }
        cluster_cell_count=as.matrix(table(mc_object@mc))
    n_cells_cluster=cluster_cell_count[clust_ord,1]

        umis=as.matrix(mat_object@mat[,names(mc_object@mc)])
        mat = umis[genes, cell_order]
        totu = colSums(umis[, cell_order])
        mat = t(t(mat)/totu)*800

        lus_1 = log2(1+7*mat[genes[gene_ord], cell_order])
        lus = apply(lus_1 - apply(lus_1, 1, median),2, function(x) pmax(x,0))
        lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x,smoothen, fill=0)))

        #lus_smoo = t(apply(lus[genes[gene_ord],cell_order], 1, function(x) rollmean(x,5, fill=0)))

    png(paste0(output_file,"_sc_cells.png"), h=height, w=width*2)
    par(mar=c(0,0,0,0))
    par(fig=c(0.2,0.90,0.1,0.9))
    shades=colorRampPalette(c("white","white","orange","red","purple","black"))(1000)
    image(t(pmin(lus_smoo,4)), col=shades,xaxt="n",yaxt="n")
        #print(quantile(lus_smoo,seq(0,1,by=0.01)))
    x=0
    for(i in 1:length(n_cells_cluster)) {
      abline(v=(n_cells_cluster[i]+x)/(sum(n_cells_cluster)-1)-1/(2*sum(n_cells_cluster)), lwd=2)      
      #mtext(clust_ord[i], side=1, at=((n_cells_cluster[i]/1.5+x)/(sum(n_cells_cluster)-1)), adj=1, las=2, line=1,cex=4)
      mtext(clust_ord[i], side=3, at=((n_cells_cluster[i]/1.5+x)/(sum(n_cells_cluster)-1)), las=2, line=1,cex=4)     
      x=x+n_cells_cluster[i]
    }


    if(!is.null(clust_col)){          
		par(fig=c(0.2,0.90,0.05,0.08),new=TRUE)
		if(is.null(names(clust_col))){names(clust_col)=as.character(clust_ord)}
        image(as.matrix(1:length(mc_object@mc)),col=clust_col[as.character(sort(mc_object@mc))], axes = F,xaxt='n',yaxt='n')
	}
	dev.off()
	}
	graphics.off()
}


scr_clustering_comparison=function(matrix1,matrix2,markers=NULL,out_fn,master_gene=NULL,w=600,h=3000,original_ordering_1st=F,original_ordering_2nd=F,pmin=1,cor_method="pearson",log_t=TRUE){
		library(RColorBrewer)

		if(!is.null(markers)){
			markers=intersect(markers,rownames(matrix1))
			markers=intersect(markers,rownames(matrix2))
			mat1=matrix1[markers,]
			mat2=matrix2[markers,]
		}
		if(is.null(markers)){
			mat1=matrix1
			mat2=matrix2
		}
		
		#colnames(scr_markers_km$centers)
		if(log_t==TRUE){
			cl_cor=cor(log2(mat1),log2(mat2),method=cor_method)
		}else{
			cl_cor=cor(mat1,mat2,method=cor_method)
		}
		
		
		rownames(cl_cor)=as.character(colnames(mat1))
		colnames(cl_cor)=as.character(colnames(mat2))
		
		
		
		if(is.null(master_gene) & !original_ordering_1st){
			hc=hclust(dist(cor(t(cl_cor))),method="ward.D2")
			tmp_cor=cl_cor[hc$order,]; rownames(tmp_cor)=seq(1,nrow(tmp_cor))
			order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
		} else if(original_ordering_1st & !original_ordering_2nd){
			hc=list()
			hc$order=colnames(mat1)
			tmp_cor=cl_cor[,]; rownames(tmp_cor)=seq(1,nrow(tmp_cor))
			order_cols=names(sort(apply(tmp_cor[,],2,function(x) which.max(x))))
		}else if(original_ordering_1st & original_ordering_2nd){
			hc=list()
			hc$order=colnames(mat1)
			order_cols=colnames(matrix2)
		} else if (!is.null(master_gene)){
			order_cols=names(sort(matrix2[master_gene,],decreasing=T))
			hc=hclust(dist(cor(t(cl_cor[,order_cols]))),method="ward.D2")
		}
		
		cor_color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
		
		png(out_fn,h=h,w=w)
		par(mar=c(0,0,0,0))
		par(fig=c(0.2,0.95,0.1,0.95))		
		image(t(pmin(pmax(cl_cor[hc$order,order_cols],0),pmin)),col=cor_color,axes=F,xaxt='n',yaxt='n')
		big_k = ncol(cl_cor)
		for(i in 1:big_k) {	abline(v=(i-0.5)/(big_k-1), lwd=0.5)}
		big_k = nrow(cl_cor)
		for(i in 1:big_k) {	abline(h=(i-0.5)/(big_k-1), lwd=0.5)}
		
		par(fig=c(0.2,0.95,0.09,0.1),new=TRUE)
		mtext(colnames(cl_cor[hc$order,order_cols]),side=1,cex=1.5,las=2,at=seq(0,1,length.out=ncol(cl_cor)))
		
		par(fig=c(0.2,0.95,0.95,0.96),new=TRUE)
		mtext(colnames(cl_cor[hc$order,order_cols]),side=3,las=2,cex=1.5,at=seq(0,1,length.out=ncol(cl_cor)))
		
		par(fig=c(0.19,0.2,0.1,0.95),new=TRUE)
		mtext(rownames(cl_cor[hc$order,order_cols]),side=2,las=2,,cex=1.5,at=seq(0,1,length.out=nrow(cl_cor)))
		
		par(fig=c(0.04,0.16,0.07,0.09),new=TRUE)
		image(x=seq(0,max(cl_cor),max(cl_cor)/(length(cor_color)-1)), 
		y=c(0,1), col=cor_color, yaxt="n", z=matrix(nrow=length(cor_color),ncol=1,data=c(1:length(cor_color))),
		cex.axis=0.9)
		mtext("Pearson corr",side=3)
		dev.off()
		#pheatmap(cl_cor[hc$order,order_cols],cluster_rows=F,cluster_cols=F)
		
		return(cl_cor[hc$order,order_cols])

}