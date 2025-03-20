RABC_path <- "C:/Users/yy193/Desktop/GEN/表达谱分析结果"
setwd(RABC_path)
RABC_name <- list.files(RABC_path , pattern = ".txt")
read <- function(x){
  read.table(x , header = T , comment.char = "" , quote = "" , strip.white = T , sep='\t')
}
listRABC <- lapply(RABC_name , read)
names(listRABC) <- basename(RABC_name)

i <- 1
col_name <- c("symbol" , "logFC" , "AveExpr" , "t" , "P.value" , "adj.P.value" , "B")
name_length <- length(RABC_name)
while (i <= name_length) {
  RABC <- data.frame(listRABC[i])
  names(RABC) <- col_name
  output_path1 <- "C:/Users/yy193/Desktop/GEN/表达谱分析差异基因"
  dir.create(output_path1, showWarnings = FALSE)
  setwd("C:/Users/yy193/Desktop/GEN/表达谱分析差异基因")
  abslogFC <- abs(RABC$logFC)
  RABC $ abslogFC <- abslogFC
  order <- RABC[order(RABC$abslogFC , decreasing = T) , ]
  quan <- order[1:round(nrow(order)*0.2) , ]
  result <- quan[which(quan$P.value < 0.05),]
  name <- RABC_name[i]
  write.table(result , file = name  , row.names = F , quote = F , sep = '\t')
  i = i + 1
}

DRABC_path <- "C:/Users/yy193/Desktop/GEN/表达谱分析差异基因"
setwd(DRABC_path)
DRABC_name <- list.files(DRABC_path , pattern = ".txt")
read1 <- function(x){
  read.table(x , header = T ,sep = '\t')
}
DRABC_list <- lapply(DRABC_name , read1)
names(DRABC_list) <- basename(DRABC_name)

GEN_path <- "C:/Users/yy193/Desktop/GEN/功能基因集"
setwd(GEN_path)
GEN_name <- list.files(GEN_path , pattern = ".txt",full.names = T)
GEN_list <- list()
i <- 1
while (i <= length(GEN_name)) {
  data <- readLines(GEN_name[i])
  split_list <- lapply(data, function(x) unlist(strsplit(x,split="\t")))
  max_length <- max(sapply(split_list, length))
  filled_list <- lapply(split_list, function(x) {
    if (length(x) < max_length) {
      c(x, rep("NA", max_length - length(x)))
    } else {
      x
    }
  })
  result_matrix <- do.call(rbind, filled_list)
  GEN_list[[i]] <- result_matrix
  i = i + 1
}
names(GEN_list) <- basename(GEN_name)

#output_path <- "C:/Users/yy193/Desktop/GEN/输出文件"
#dir.create(output_path, showWarnings = FALSE)
#i <- 1
#while (i <= length(GEN_list)) {
#  output_file <- file.path(output_path, names(GEN_list)[i])
#  write.table(GEN_list[[i]], file = output_file, row.names = F, col.names = F, sep = "\t",quote = F)
#  i=i+1
#}

all_homeworks <- list()
i <- 1
while (i <= length(GEN_list)) {
  Geneset <- GEN_list[[i]]
  gen_set_name <- sub("\\.txt$", "",names(GEN_list)[i])
  a <- 1
  while (a <= length(listRABC)) {
    DRABC <- DRABC_list[[a]]#差异基因
    drabc_set_name <- sub("\\.txt$", "",names(DRABC_list)[a])
    Gen <- t(Geneset)
    index <- nrow(Gen)
    Gen_1 <- Gen[3:index,]#功能基因
    RABC_GEN <- DRABC$symbol#差异基因
    RABC1 <- listRABC[[a]]
    RABC1_GEN <- RABC1$symbol#背景基因 
    col <- ncol(Gen_1)
    homework <- data.frame(row.names = c("DataID","GeneID","FuncGene_DEGs","FuncGene_nonDEGs","nonFuncgene_DEGs","nonFuncgene_nonDEGs","pvalue"))
    b <- 1
    while (b <= col) {
      Gen_101 <- Gen_1[,b] #功能基因
      I <- intersect(RABC_GEN,Gen_101) #功能基因与差异基因交集
      I1 <- intersect(RABC1_GEN,Gen_101) #功能基因与背景基因交集
      G <- length(I1) #功能基因个数
      M <- length(RABC1_GEN) #背景基因个数
      D <- length(RABC_GEN) #差异基因个数
      DEG <- length(I) #功能基因与差异基因交集个数
      p = 1-phyper(DEG,D,M-D,G)
      DataID <- drabc_set_name
      GeneID <- Gen[1,b]
      r <- c(DataID,GeneID,DEG,G-DEG,D-DEG,M-G-D+DEG,p)
      homework[,b] <- r
      b = b + 1
    }
    
    all_homeworks[[paste0(gen_set_name, "_", drabc_set_name)]] <- homework
    a = a + 1
  }
  i = i + 1  
}

output_path <- "C:/Users/yy193/Desktop/GEN/homework"
dir.create(output_path, showWarnings = FALSE) 

for (name in names(all_homeworks)) {
  file_path <- file.path(output_path, paste0(name, ".txt")) 
  write.table(t(all_homeworks[[name]]), file = file_path, row.names = F, col.names = TRUE, sep = "\t", quote = FALSE)
}

