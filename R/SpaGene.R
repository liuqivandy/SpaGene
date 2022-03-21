#' Identify spatially variable genes
#'
#' @description Identify spatial variable genes based on spatial connectness of spots with high expression compared to random permutation

#' @param expr gene expression matrix, the row is the gene and the column is the spot/cell
#' @param location location matrix, the row number of location should match the column number of expr
#' @param normalize whether to normalize the data (default: TRUE)
#' @param topn the number of spots/cells considered high expression (default: 20 percent of the total spots/cells)
#' @param knn the number of nearest neighbours to search (default: 8)
#' @param perm the number of random permutations (default: 500)
#' @param minN the minimum number of spots/cells with gene expression. Genes expressed equal to or less than minN spots/cells are excluded (default:0)
#' @param sizefactor the size factor for normalization (default:10000)
#' @return a list containing results of each gene (spagene_res) and normalized gene expression matrix (normexp)

#' @export

SpaGene <- function(expr,location,normalize=T,topn=floor(0.2*dim(location)[1]),knn=8,perm=500,minN=0,sizefactor=10000) {
  set.seed(1)
  expr<-expr[Matrix::rowSums(expr>0)>minN,]

  ncell<-dim(location)[1]
  ngene<-dim(expr)[1]

  if (dim(expr)[2]!=ncell) {stop("the ncol of expr should match the nrow of location")}

  if (is.null(rownames(expr))){rownames(expr)<- paste0("gene",1:ngene)}

   nnmatrix<-RANN::nn2(location,k=knn)$nn.idx

   rand_result<-unlist(lapply(1:perm,function(x){ind<-sample(1:ncell,topn);return(Caldegree(ind,nnmatrix,knn))}))

   mean_rand<-mean(rand_result)
   sd_rand<-sd(rand_result)

   spagene_res<-data.frame(score=rep(NA,ngene),row.names=rownames(expr),stringsAsFactors = FALSE)


  if (is(expr,"sparseMatrix")){
    exprt<-Matrix::t(expr)
    colind<-exprt@i+1
    dp<-exprt@p
    expval<-exprt@x

    if (normalize==TRUE) {
       lib_size<-Matrix::rowSums(exprt)
       expval<-log(expval/lib_size[colind]*sizefactor+1)
       exprt@x<-expval
       expr<-Matrix::t(exprt)
     }



    for (geneind in 1:ngene) {

      geneexp<-rep(0,ncell)
      subind<-(dp[geneind]+1):dp[geneind+1]
      geneexp[colind[subind]]<-expval[subind]
      highind<-order(geneexp,sample(ncell,ncell),decreasing=T)[1:topn]
      spagene_res$score[geneind]<-Caldegree(highind,nnmatrix,knn)


    }
  } else{

         if (normalize==TRUE) {
             expr<-log(t(t(expr)/(colSums(expr))*sizefactor)+1)
          }

        for (geneind in 1:ngene) {
              geneexp<-expr[geneind,]

              highind<-order(geneexp,sample(ncell,ncell),decreasing=T)[1:topn]
              spagene_res$score[geneind]<-Caldegree(highind,nnmatrix,knn)
         }
    }




  spagene_res$zval<-(spagene_res$score-mean_rand)/sd_rand


  spagene_res$pval<-pnorm(spagene_res$score,mean=mean_rand,sd=sd_rand)
  spagene_res$adjp<-p.adjust(spagene_res$pval,method="BH")
  return(list(normexp=expr,spagene_res=spagene_res))
}

#' Find spatial patterns
#' @description Find spatial patterns
#' @param spagene_res result from SpaGene, a list containing normexp and spagene_res
#' @param cutoff the adjp cutoff to select spatially variable genes (default: 0.01)
#' @param nPattern the number of patterns (default:8)
#' @return a list containing the pattern (pattern), gene similarity with the pattern (genepattern), and the pattern weight (patternw)
#' @export

FindPattern<-function(spagene_res,cutoff=0.01,nPattern=8){

  genes<-rownames(spagene_res$spagene_res)[spagene_res$spagene_res$adjp<cutoff]
  data<-spagene_res$normexp
  data_g<-as.matrix(data[rownames(data)%in%genes,])
  set.seed(16)
  model<-RcppML::nmf(data_g,nPattern,verbose = FALSE)
  cellload<-model$h
  patternw<-model$d
  genew<-model$w
  rownames(genew)<-rownames(data_g)
  rownames(cellload)<-paste0("Pattern",1:nPattern)
  genepattern<-cor(t(data_g),t(cellload),method="spearman")
  return(list(pattern=cellload,genepattern=genepattern,patternw=patternw))
}


#' Identify spatially colocalized ligand-receptor pairs
#'
#' @description Identify spatially colocalized ligand-receptor pairs

#' @param expr gene expression matrix, the row is the gene and the column is the spot/cell
#' @param location location matrix, the row number of location should match the column number of expr
#' @param normalize whether to normalize the data (default: TRUE)
#' @param topn the number of spots/cells considered high expression (default: 20 percent of the total spots/cells)
#' @param knn the number of nearest neighbours to search (default: 8)
#' @param perm the number of random permutations (default: 500)
#' @param minN the minimum number of spots/cells with gene expression. Genes expressed equal to or less than minN spots/cells are excluded (default:0)
#' @param sizefactor the size factor for normalization (default:10000)
#' @param LRpair ligand-receptor pair
#' @return a data frame containing the result of each ligand-receptor pair

#' @export



SpaGene_LR<-function(expr,location,normalize=T, topn=floor(0.2*dim(location)[1]),knn=8,perm=500,minN=0,sizefactor=10000,LRpair=LRpair) {

  set.seed(1)
  expr<-expr[Matrix::rowSums(expr>0)>minN,]

  ncell<-dim(location)[1]


  if (dim(expr)[2]!=ncell) {stop("the ncol of expr should match the nrow of location")}

  if (is.null(rownames(expr))){rownames(expr)<- paste0("gene",1:ngene)}

  nnmatrix<-RANN::nn2(location,k=knn)$nn.idx

  rand_result<-unlist(lapply(1:perm,function(x){return(Caldegree_pair(sample(1:ncell,topn),sample(1:ncell,topn),nnmatrix,knn))}))

  mean_rand<-mean(rand_result)
  sd_rand<-sd(rand_result)

  if (normalize==TRUE) {
    expr<-Matrix::t(Matrix::t(expr)/(Matrix::colSums(expr))*sizefactor)
  }

  npair<-dim(LRpair)[1]
  lr_result<-data.frame(score=rep(NA,npair),comm=rep(NA,npair),row.names=rownames(LRpair),stringsAsFactors = FALSE)

  for (pairid in 1:dim(LRpair)[1]) {
    ligand<-LRpair[pairid,1]
    receptor<-LRpair[pairid,2]
    if (sum(rownames(expr) %in% c(ligand,receptor))==2) {
      ligandind<-order(expr[rownames(expr)==ligand,],sample(ncell,ncell),decreasing=T)[1:topn]
      receptorind<-order(expr[rownames(expr)==receptor,],sample(ncell,ncell),decreasing=T)[1:topn]
      lr_result$score[pairid]<-Caldegree_pair(ligandind,receptorind,nnmatrix,knn)
      lr_result$comm[pairid]<-length(intersect(ligandind,receptorind))
    }
  }



  lr_result<-lr_result[!is.na(lr_result$score),]

  lr_result$zval<-(lr_result$score-mean_rand)/sd_rand


  lr_result$pval<-pnorm(lr_result$score,mean=mean_rand,sd=sd_rand)
  lr_result$adjp<-p.adjust(lr_result$pval,method="BH")

  return(lr_result)

}


#' Plot patterns
#'
#' @description plot patterns from spatially variable genes

#' @param pattern pattern result from FindPattern
#' @param location location matrix
#' @param max.cutoff the maximum value cutoff (default:0.9)
#' @param pt.size the point size (default:2)
#' @param alpha.min the alpha value of the minimum value (default:0.1)
#' @return a list of ggplot
#' @export

PlotPattern<-function(pattern,location,max.cutoff=0.9,pt.size=2,alpha.min=0.1) {

  if(!requireNamespace("RColorBrewer", quietly = TRUE)){install.packages("RColorBrewer")}

   colnames(location)<-c("x","y")
   npattern<-dim(pattern$pattern)[1]
   plist<-list()


  for (i in 1:npattern) {

    feature=pattern$pattern[i,]
    max.use<-quantile(feature,max.cutoff)
    feature[feature>max.use]<-max.use
    alpha=(feature-min(feature))/(max(feature)-min(feature))*(1-alpha.min)+alpha.min
    tmp<-as.data.frame(cbind(location,exp=feature,alpha=alpha))

    p1<-ggplot(tmp,aes(x=x,y=y,col=exp,alpha=alpha))+geom_point(size=pt.size)+scale_y_reverse()+scale_color_gradientn(colours=rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu")))+xlab("")+ylab("")+theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+guides(color = "none", alpha = "none")+ggtitle(paste0("Pattern",i))
    plist[[i]]<-p1

  }
  patchwork::wrap_plots(plist)
}

#' plot one specific ligand-receptor pair
#'
#' @description plot one specific ligand-receptor pair to find the colocalized region

#' @param expr gene expression matrix, the row is the gene and the column is the spot/cell
#' @param location location matrix, the row number of location should match the column number of expr
#' @param normalize whether to normalize the data (default: TRUE)
#' @param topn the number of spots/cells considered high expression (default: 20 percent of the total spots/cells)
#' @param knn the number of nearest neighbours to search (default: 8)
#' @param LRpair the ligand-receptor pair for plot
#' @param pt.size the point size (default:2)
#' @param alpha.min the alpha for the minimum value (default:0.1)
#' @return a data frame containing the result of each ligand-receptor pair
#' @export

plotLR<-function(expr,location,normalize=T,topn=floor(0.2*dim(location)[1]),knn=8,LRpair=c("Ptn","Ptprz1"),pt.size=2,alpha.min=0.1){
  if (sum(rownames(expr) %in% LRpair)!=2) { stop("ligand or receptor are not expressed")}
  nnmatrix<-RANN::nn2(location,k=knn)$nn.idx
  countsum<-Matrix::colSums(expr)

  ncell<-dim(expr)[2]
  if (normalize==TRUE) {
    expr<-Matrix::t(log(Matrix::t(expr)/countsum*median(countsum)+1))
  }


  ligand<-expr[LRpair[1],]
  receptor<-expr[LRpair[2],]
  LRexp<-rbind(ligand,receptor)
  neighexp<-apply(nnmatrix,1,function(x){apply(LRexp[,x[2:knn]],1,max)})

  LRexp<-t(scale(t(LRexp)))
  neighexp<-t(scale(t(neighexp)))
  LRexp[LRexp<0]<-0
  neighexp[neighexp<0]<-0
  LRadd<-pmax(LRexp[1,]*neighexp[2,],LRexp[2,]*neighexp[1,])

  if (sum(ligand>0)>topn) {n1<-order(ligand,sample(ncell,ncell),decreasing=T)[1:topn]} else{n1<-which(ligand>0)}
  if (sum(receptor>0)>topn) {n2<-order(receptor,sample(ncell,ncell),decreasing=T)[1:topn]} else{n2<-which(receptor>0)}
  expcol<-rep(0,ncell)
  expcol[n1]<-1
  expcol[n2]<-2
  expcol[intersect(n1,n2)]<-3
  tmp<-data.frame(x=location[,1],y=location[,2],Exp=as.factor(expcol))
  tmpLRadd<-data.frame(x=location[,1],y=location[,2],LR=LRadd)

  alpha=(LRadd-min(LRadd))/(max(LRadd)-min(LRadd))*(1-alpha.min)+alpha.min

  p1<-ggplot(tmp,aes(x=x,y=y,col=Exp))+geom_point(size=pt.size)+scale_color_manual(values=c("gray","red","green","blue"),labels=c("Both low","Ligand high","Recptor High","Both High"))+ggtitle(paste0(LRpair,collapse="_"))+xlab("")+ylab("")+theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
  p2<-ggplot(tmpLRadd,aes(x=x,y=y,col=LR))+geom_point(size=pt.size,alpha=alpha)+scale_color_gradient2(midpoint=quantile(LRadd,probs=0.5),low="gray",high="red",mid="gray")+xlab("")+ylab("")+theme(axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+labs(color = "LR")
  p1+p2&scale_y_reverse()
}


Caldegree<-function(nodelist,nnmatrix,knnnum) {
  nodenum<-length(nodelist)

  nnmatrix_sub<-nnmatrix[nodelist,-1]

  deg<-rep(0,nodenum)
  matchres<-matrix(match(nnmatrix_sub,nodelist),ncol=knnnum-1,byrow=F)
  ind<-!is.na(matchres)
  deg<-deg+rowSums(ind)

  matchres<-matchres[ind]
  for (i in 1:length(matchres)) deg[matchres[i]]<-deg[matchres[i]]+1

  dis<-sum(cumsum(tabulate(deg+1,nbins=2*knnnum+1)/nodenum))
  return(dis)
}

Caldegree_pair<-function(nodelist1,nodelist2,nnmatrix,knnnum) {
  nodenum<-length(nodelist1)

  nnmatrix_sub<-nnmatrix[nodelist1,-1]

  deg1<-rep(0,nodenum)
  matchres<-matrix(match(nnmatrix_sub,nodelist2),ncol=knnnum-1,byrow=F)
  ind<-!is.na(matchres)
  deg1<-deg1+rowSums(ind)
  deg2<-rep(0,length(nodelist2))
  matchres<-matchres[ind]
  for (i in 1:length(matchres)) deg2[matchres[i]]<-deg2[matchres[i]]+1

  deg<-c(deg1,deg2)

  dis<-sum(cumsum(tabulate(deg+1,nbins=2*knnnum+1)/nodenum))
  return(dis)
}



