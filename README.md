## 环境配置  
在安装R包之前需要先安装依赖包：openxlsx、seqinr、plyr、randomForestSRC、glmnet、plsRglm、gbm、caret、mboost、e1071、BART、MASS、snowfall、xgboost、ComplexHeatmap、RColorBrewer、pROC  
## OmniLearn  
*OmniLearn* 是一个功能强大且用户友好的 R 包，专为构建基于 *113 种机器学习算法* 的诊断模型而设计。该包集成了特征选择、模型训练、性能评估和可视化工具，为处理复杂、高维数据集的研究人员提供了一站式解决方案。  
## 使用方法  
#### 1、查看支持的方法
```R
## 查看支持的方法
print(methods_use[1:10])
```
![alt text]([image.png](https://github.com/Feng-Rommel/OmniLearn/blob/main/fig/image.png))

#### 2、数据处理  
*训练数据*  
```R
# 训练集表达谱是行为基因（感兴趣的基因集），列为样本的表达矩阵（基因名与测试集保持相同类型，如同为SYMBOL或ENSEMBL等）
Train_expr <- read.table("./test_data/Training_expr.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# 行为样本，列包含至少一个需要预测的二分类变量(仅支持[0，1]格式)
Train_class <- read.table("./test_data/Training_class.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
# 统一表达谱的样本和类别信息的样本
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = F]
```
![alt text](image-1.png)
![alt text](image-2.png)

*验证队列*  
```R
## 处理与训练数据相同  
Test_expr <- read.table("./test_data/Testing_expr.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.table("./test_data/Testing_class.txt", header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam]; Test_class <- Test_class[comsam,,drop = F]
# 提取相同基因
comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
Test_expr <- t(Test_expr[comgene,]) # 输入模型的表达谱行为样本，列为基因
```
*验证队列的样本信息多了Cohort列，因为可能不止一个验证集，所以使用Cohort来区分不同的验证集*  
![alt text](image-3.png)
![alt text](image-4.png)

*数据标准化（可选）*
```R
Train_set = scaleData(data = Train_expr, centerFlags = T, scaleFlags = T) 
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort)) # 注意测试集标准化顺序与此一致
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = T, scaleFlags = T)
```

#### 3、模型训练
*我们提供了一个灵活的模型训练方法，如果想直接完整运行113种算法可以直接使用下面的代码*  
```R
## 首先使用筛选算法进行基因筛选
classVar = "outcome" # 设置所要预测的变量名（仅支持[0,1]二元变量格式）
min.selected.var = 3 # 设置模型最少纳入的变量数
methods = methods_use
## Pre-training 将各方法所用到的变量筛选过程汇总，以减少计算量
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+") # 检视所有方法，分析各方法是否需要进行变量预筛选(pre-training)
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1]) # 删除各方法用于构建分类模型的算法，保留用于变量筛选的算法
preTrain.method = unique(unlist(preTrain.method)) # 汇总所有变量筛选算法，去除重复计算
preTrain.var <- list() # 用于保存各算法筛选的变量
set.seed(seed = 777) # 设置建模种子，使得结果可重复
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, # 变量筛选所需要的机器学习方法
                                 Train_set = Train_set, # 训练集有潜在预测价值的变量
                                 Train_label = Train_class, # 训练集分类标签
                                 mode = "Variable",       # 运行模式，Variable(筛选变量)和Model(获取模型)
                                 classVar = classVar) # 用于训练的分类变量，必须出现在Train_class中
}
preTrain.var[["simple"]] <- colnames(Train_set)# 记录未经筛选的变量集（以便后续代码撰写），可视为使用simple方法（无筛选功能）的变量筛选结果
```

```R
## 进行模型训练
#### Model training
min.selected.var = 3 # 设置模型最少纳入的变量数
model <- list() # 用于保存各模型的所有信息
set.seed(seed = 777) # 设置建模种子，使得结果可重复
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method # 本轮算法名称
  method <- strsplit(method, " \\+ ")[[1]] # 各步骤算法名称
  
  if (length(method) == 1) method <- c("simple", method) # 如果本方法没有预筛选变量，则认为本方法使用simple方法进行了变量筛选
  
  Variable = preTrain.var[[method[1]]] # 根据方法名称的第一个值，调用先前变量筛选的结果
  ## 如果输入的特征数小于预设值则不训练
  if(length(Variable)<min.selected.var){
    next
  }
  Train_label = Train_class            # 所以此处需要修改变量名称，以免函数错误调用对象
  model[[method_name]] <- RunML(method = method[2],        # 根据方法名称第二个值，调用构建的函数分类模型
                                Train_set = Train_set[,Variable],     # 训练集有潜在预测价值的变量
                                Train_label = Train_label, # 训练集分类标签
                                mode = "Model",            # 运行模式，Variable(筛选变量)和Model(获取模型)
                                classVar = classVar)       # 用于训练的分类变量，必须出现在Train_class中
  
  # 如果最终模型纳入的变量数小于预先设定的下限，则认为该算法输出的结果是无意义的
  if(length(ExtractVar(model[[method_name]])) < 2) {
    model[[method_name]] <- NULL
  }
}
saveRDS(model, "./model.rds")
```

## 4、模型应用
```R
# 根据给定表达量计算样本风险评分
# 预测概率
methodsValid <- names(model)
RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, "./RS_mat.txt",sep = "\t", row.names = T, col.names = NA, quote = F) # 输出风险评分文件

# 根据给定表达量预测分类
Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set))
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
#Class_mat <- cbind.data.frame(Test_class, Class_mat[rownames(Class_mat),]) # 若要合并测试集本身的样本信息文件可运行此行
write.table(Class_mat, "./Class_mat.txt", # 测试集经过算法预测出的二分类结果
            sep = "\t", row.names = T, col.names = NA, quote = F)
```

## 5、提取特征
```R
# 提取所筛选的变量（列表格式）
fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

# 提取所筛选的变量（数据框格式）
fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, "./fea_df.txt", # 两列，包含算法以及算法所筛选出的变量
            sep = "\t", row.names = F, col.names = T, quote = F)
```

## 6、模型评估
```R
# 对各模型计算C-index
AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(fit = model[[method]],     # 分类预测模型
                                Test_set = Test_set,      # 测试集预测变量，应当包含训练集中所有的变量，否则会报错
                                Test_label = Test_class,   # 训练集分类数据，应当包含训练集中所有的变量，否则会报错
                                Train_set = Train_set,    # 若需要同时评估训练集，则给出训练集表达谱，否则置NULL
                                Train_label = Train_class, # 若需要同时评估训练集，则给出训练集分类数据，否则置NULL
                                Train_name = "TCGA",       # 若需要同时评估训练集，可给出训练集的标签，否则按“Training”处理
                                cohortVar = "Cohort",      # 重要：用于指定队列的变量，该列必须存在且指定[默认为“Cohort”]，否则会报错
                                classVar = classVar)       # 用于评估的二元分类变量，必须出现在Test_class中
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, "./AUC_mat.txt",
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------

AUC_mat <- read.table("./AUC_mat.txt",sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           # 计算每种算法在所有队列中平均AUC
avg_AUC <- sort(avg_AUC, decreasing = T)     # 对各算法AUC由高到低排序
AUC_mat <- AUC_mat[names(avg_AUC), ]      # 对AUC矩阵排序
fea_sel <- fea_list[[rownames(AUC_mat)[1]]] # 最优模型（测试集AUC均值最大）所筛选的特征
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3)) # 保留三位小数

if(ncol(AUC_mat) < 3) { # 如果用于绘图的队列小于3个
  CohortCol <- c("red","blue") # 则给出两个颜色即可（可自行替换颜色）
} else { # 否则通过brewer.pal赋予超过3个队列的颜色
  CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired") # 设置队列颜色
}
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, # 主矩阵
                    avg_AUC, # 侧边柱状图
                    CohortCol, "steelblue", # 列标签颜色，右侧柱状图颜色
                    cellwidth = cellwidth, cellheight = cellheight, # 热图每个色块的尺寸
                    cluster_columns = F, cluster_rows = F) # 是否对行列进行聚类

pdf(file.path(fig.path, "AUC.pdf"), width = cellwidth * ncol(AUC_mat) + 3, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())
```
![alt text](image-5.png)



