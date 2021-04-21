
#' main function of this package, to calculate cell proportion
#'
#' @param bulk.eset ExpressionSet for bulk data, raw counts as the input
#' @param method Three deconvolution methods are provided: MuSiC, NNLS and bseqsc (based on CIBERSORT), MuSiC is the defual method
#' @return est.prop, A cell proportion matrix, each row represent a sample, each col is a cell type
#' @import tidyverse, reshape2, xbioc, MuSiC,bseqsc
#' @export
#'
BronCell.prop<-function(bulk.eset,method="MuSiC"){

  if(!exists("sc.basis")){
    data(sc.basis)
  }

  if(method=="MuSiC"){
    Est.prop<-music_prop_bron(bulk.eset = bulk.eset, sc.basis = sc.basis, clusters = 'cellType',
                                samples = 'sampleID', verbose = F)

    est_prop<-as.data.frame(Est.prop$Est.prop.weighted)

  } else if (method=="nnls"){
    Est.prop<-music_prop_bron(bulk.eset = bulk.eset, sc.basis = sc.basis, clusters = 'cellType',
                                samples = 'sampleID', verbose = F)

    est_prop<-as.data.frame(Est.prop$Est.prop.allgene)
  } else if (method=="bseq"){
    fit <- bseqsc_proportions(bulk.eset, sc.basis$M.theta, verbose = TRUE)
    est_prop<-as.data.frame(t(coef(fit)))
  }

  return(est_prop)
}


#' plot cell proporions for each bronchial cell types
#' @param est_prop cell proportion matrix, estimated by funtion BronCell.prop
#' @import tidyverse, reshape2
#' @export
BronCell.plot<-function(est_prop){
  est_prop$id<-rownames(est_prop)
  df1<-melt(est_prop,id.vars = "id")
  colnames(df1)<-c("ID","Cell_type","Estimated_proportion")
  plot<-ggplot(df1, aes(x = Cell_type, y = Estimated_proportion,fill=Cell_type)) +
    geom_boxplot()+geom_jitter(size=0.5)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x=element_text(angle=45, hjust=1))

  return(plot)
}


#' modified function of "music_prop" from MuSiC package using the signature matrix as the input
#' original code is on: https://github.com/xuranw/MuSiC
#'
#' This function is to calculate the MuSiC deconvolution proportions
#'
#' @param bulk.eset ExpressionSet for bulk data
#' @param sc.basis Gene signature matrix generated from bronchial biopsy scRNAseq data, using function music_basis from MuSiC package
#' @param clusters character, the phenoData of single cell dataset used as clusters;
#' @param samples character,the phenoData of single cell dataset used as samples;
#' @param select.ct vector of cell types, default as NULL. If NULL, then use all cell types provided by single cell dataset;
#' @param cell_size data.frame of cell sizes. 1st column contains the names of cell types, 2nd column has the cell sizes per cell type. Default as NULL. If NULL, then estimate cell size from data;
#' @param ct.cov logical. If TRUE, use the covariance across cell types;
#' @param verbose logical, default as TRUE.
#' @param iter.max numeric, maximum iteration number
#' @param nu regulation parameter, take care of weight when taking recipical
#' @param eps Thredshold of convergence
#' @param centered logic, substract avg of Y and D
#' @param normalize logic, divide Y and D by their standard deviation
#' @return a list with elements:
#'    * Estimates of MuSiC
#'    * Estimates of NNLS
#'    * Weight of MuSiC
#'    * r.squared of MuSiC
#'    * Variance of MuSiC estimates
#' @export
music_prop_bron = function(bulk.eset, sc.basis, clusters, samples, select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE,
                           iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE, ... ){

  bulk.gene = rownames(bulk.eset)[rowMeans(exprs(bulk.eset)) != 0]
  bulk.eset = bulk.eset[bulk.gene, , drop = FALSE]
  markers<-rownames(sc.basis$Disgn.mtx)
  cm.gene = intersect( rownames(sc.basis$Disgn.mtx), bulk.gene )

  if(length(cm.gene)< 0.2*length(markers))
    stop("Too few common genes!")

  if(verbose){message(paste('Used', length(cm.gene), 'common genes...'))}

  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx)); m.bulk = match(cm.gene, bulk.gene)
  D1 = sc.basis$Disgn.mtx[m.sc, ];
  M.S = colMeans(sc.basis$S, na.rm = T);

  if(!is.null(cell_size)){
    if(!is.data.frame(cell_size)){
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }else if(sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))){
      stop("Cell type names in cell_size must match clusters")
    }else if (any(is.na(as.numeric(cell_size[, 2])))){
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }

  Yjg = relative.ab(exprs(bulk.eset)[m.bulk, ]); N.bulk = ncol(bulk.eset);
  if(ct.cov){
    Sigma.ct = sc.basis$Sigma.ct[, m.sc];

    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL

    for(i in 1:N.bulk){
      if(sum(Yjg[, i] == 0) > 0){
        D1.temp = D1[Yjg[, i]!=0, ];
        Yjg.temp = Yjg[Yjg[, i]!=0, i];
        Sigma.ct.temp = Sigma.ct[, Yjg[,i]!=0];
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
      }else{
        D1.temp = D1;
        Yjg.temp = Yjg[, i];
        Sigma.ct.temp = Sigma.ct;
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
      }

      lm.D1.weighted = music.iter.ct(Yjg.temp, D1.temp, M.S, Sigma.ct.temp, iter.max = iter.max,
                                     nu = nu, eps = eps, centered = centered, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg)); weight.gene.temp[Yjg[,i]!=0] = lm.D1.weighted$weight.gene;
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }else{
    Sigma = sc.basis$Sigma[m.sc, ];

    valid.ct = (colSums(is.na(Sigma)) == 0)&(colSums(is.na(D1)) == 0)&(!is.na(M.S))

    if(sum(valid.ct)<=1){
      stop("Not enough valid cell type!")
    }

    if(verbose){message(paste('Used', sum(valid.ct), 'cell types in deconvolution...' ))}

    D1 = D1[, valid.ct]; M.S = M.S[valid.ct]; Sigma = Sigma[, valid.ct];

    Est.prop.allgene = NULL
    Est.prop.weighted = NULL
    Weight.gene = NULL
    r.squared.full = NULL
    Var.prop = NULL
    for(i in 1:N.bulk){
      if(sum(Yjg[, i] == 0) > 0){
        D1.temp = D1[Yjg[, i]!=0, ];
        Yjg.temp = Yjg[Yjg[, i]!=0, i];
        Sigma.temp = Sigma[Yjg[,i]!=0, ];
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
      }else{
        D1.temp = D1;
        Yjg.temp = Yjg[, i];
        Sigma.temp = Sigma;
        if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
      }

      lm.D1.weighted = music.iter(Yjg.temp, D1.temp, M.S, Sigma.temp, iter.max = iter.max,
                                  nu = nu, eps = eps, centered = centered, normalize = normalize)
      Est.prop.allgene = rbind(Est.prop.allgene, lm.D1.weighted$p.nnls)
      Est.prop.weighted = rbind(Est.prop.weighted, lm.D1.weighted$p.weight)
      weight.gene.temp = rep(NA, nrow(Yjg)); weight.gene.temp[Yjg[,i]!=0] = lm.D1.weighted$weight.gene;
      Weight.gene = cbind(Weight.gene, weight.gene.temp)
      r.squared.full = c(r.squared.full, lm.D1.weighted$R.squared)
      Var.prop = rbind(Var.prop, lm.D1.weighted$var.p)
    }
  }
  colnames(Est.prop.weighted) = colnames(D1)
  rownames(Est.prop.weighted) = colnames(Yjg)
  colnames(Est.prop.allgene) = colnames(D1)
  rownames(Est.prop.allgene) = colnames(Yjg)
  names(r.squared.full) = colnames(Yjg)
  colnames(Weight.gene) = colnames(Yjg)
  rownames(Weight.gene) = cm.gene
  colnames(Var.prop) = colnames(D1)
  rownames(Var.prop) = colnames(Yjg)

  return(list(Est.prop.weighted = Est.prop.weighted, Est.prop.allgene = Est.prop.allgene,
              Weight.gene = Weight.gene, r.squared.full = r.squared.full, Var.prop = Var.prop))
}


