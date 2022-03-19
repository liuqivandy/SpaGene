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
#' @param sizefactor the size factor for normalization
#' @return a list containing results of each gene (spagene_res) and normalized gene expression matrix (normexp)

#' @export

SpaGene <- function(expr,location,normalize=T,topn=floor(0.2*dim(location)[1]),knn=8,perm=500,minN=0,sizefactor=10000) {

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


  spagene_res$score<-spagene_res$score[!is.na(spagene_res$score)]

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

  model<-RcppML::nmf(data_g,nPattern)
  cellload<-model$h
  patternw<-model$d
  genew<-model$w
  rownames(genew)<-rownames(data_g)
  rownames(cellload)<-paste0("Pattern",1:nPattern)
  genepattern<-cor(t(data_g),t(cellload),method="spearman")
  return(list(pattern=cellload,genepattern=genepattern,patternw=patternw))
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




