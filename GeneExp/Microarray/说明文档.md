# microarray-analysis
From .CEL to expression matrix  
This document is a pre-proposs pipeline of Affymetrix microarray.

以 GSE35977为例的处理流程

####################################  
一 数据下载

		1. GEO数据库下载
			1. 简单版可以直接使用拼接方式 注:sample information可用matrix中注释行转置  
			   分组信息可http://www.ncbi.nlm.nih.gov/geo/browse/?view=samples&series=#GSE号码#&mode=csv下载 
		2. r包GEOquery
        		getGEO(GEO = NULL, filename = NULL, destdir = tempdir(), GSElimits=NULL,GSEMatrix=TRUE,AnnotGPL=FALSE)

二 数据处理

-------------------------------------------------AFFY 预处理---------------------------------------------------------

	* 芯片质量检测
		1. 若提供若提供了image file 应对其进行观察，看是否有因实验操作不当而产生的气泡 (最直接的质量控制)
		2. 下载所有原始文件(.CEL)并解压到统一同一文件夹中
		3. R 操作进行质量分析
		     a. 芯片扫描图像（灰度 ）b.灰度值箱线图 c. histogram曲线 d. MA-plot分析 e. RNA降解分析
		4. 删除有问题的芯片(芯片筛选)

	* 数据预处理[详情见]:
		1. 数据预处理一般有三个步骤: 

			* - 背景处理（background adjustment）
			* - 归一化处理（normalization，或称为“标准化处理”）
			* - 汇总（summarization）。

		2. 芯片基本情况初探(可了解pm(pefect match probe)以及mm(mismatch probe)的个数, 探针所映射的基因名)
		3. 背景处理 

			* rma  
				* RMA方法的原理比较复杂，  
				可以参看文献：  
				R. A. Irizarry, B. Hobbs, F. Collin, et al. 
				Exploration, normalization, and summaries of high density oligonucleotide array probe   
				level data. Biostatistics, 4:249–64, 2003b. 11, 18, 27, 232, 241, 432, 443。
				* RMA方法仅使用PM探针数据，背景调整后MM的信号值不变。
				* 得到表达量的信息

			* Mas5
				* MAS方法将芯片分为k（默认值为16）个网格区域，用每个区域使用信号强度最低的2%探针去计算背景值和噪声
				* MAS方法应用后PM和MM的信号强度都被重新计算。
				* 得到detection p-value

		4. 归一化处理(根据具体情况进行选择)

			* 线性缩放方法
				* 线性缩放方法以第一块芯片为参考，它的数值没有被处理，而其他芯片都被缩放了。  
				  对同一块芯片，不同探针的缩放倍数是一个常数。PM和MM的缩放方法完全一样。  
				  这是Affy公司在其软件（4.0和5.0版本）中使用的方法。这种方法先选择一个芯片作为参考，  
				  将其他芯片和参考芯片逐一做线性回归分析，用回归直线（没有截距）对其他芯片的信号值做缩放。  
				  Affy公司的软件做回归分析前去除了2%最强和最弱信号。  
				  
			* 非线性缩放方法
				* 非线性拟合时不是取整张芯片而仅取部分（一列）作为基线。它的数值没有被处理，而其他数据都被缩放了。

			* 分位数（quantile）方法--分位数标准化（Quantile Normalization）
				* 一般芯片的杂交实验很容易产生误差，所以经常一个样本要做 3~6 次的重复实验。  
				  平行实验间的数据差异可以通过 Quantile Normalization去处掉。
			  	  总平行实验的前提条件是假设 n次实验的数据具有相同的分布，其算法主要分为三步  
					（1）对每张芯片的数据点排序。  
					（2）求出同一位置的几次重复实验数据的均值，并用该均值代替该位置的基因的表达量。  
					（3）将每个基因还原到本身的位置上。
					
			* 其他 （如循环局部加权回归法（Cyclic loess）和 Contrasts方法）
				* 芯片内的数据标准化，主要是去除每张芯片的系统误差，这种误差主要是由荧光染色差异，  
				  点样机器(arrayer print-tip)，或者杂交试验所产生的，通过标准化，使每个基因的表达点都具有独立性。

		5. 汇总 最后一步汇总是使用合适的统计方法通过probeset（包含多个探针）的杂交信号计算出计算表达量。
			* 需要注意的是computeExprSet函数除需要指定统计方法外还需要指定PM校正的方式   
			  常用的汇总方法是medianpolish, liwong和mas. liwong方法仅使用PM做背景校正（pmcorrect.method="pmonly")

	* 数据过滤
		1. 探针过滤
			* 去掉含SNP的探针[在哪里下载? 官网下载]
			* detection- P<0.06的数据量≥80%[]
		2. 样本过滤
			* 利用PCA对样本进行主成分分析, 删除outliar
		3. 去掉AFFX开头的做质控的探针
			* genefilter 包中 nsFilter函数，注意参数设置
		4. 去掉探针匹配效果不佳的探针
			* 在芯片的annotation file 中找到探针匹配效果不佳探针
		5. 去掉gene-assignement为空的数据
		6. impute missing value
			* x[is.na(x)]=0.0001 就是把missing value替换成0.0001;
			* 用impute把missing value按照数据分布补回来 ，数值大小与缺失值周围数值大小有关 Impute.knn()
			* 其他方法待补充

	* 批次校正
		1. 读取探针及样本过滤后文件,检查是否有missing value
			* missing value 处理方法：imputation
		2. 制作batch文件【利用excel】
			* 复制sample名
			* 转置粘贴到新的工作表【使得原本一行的数据放到一列】
			* 观察sample名，利用数据-分列-分隔符选择“_”进行分列
			* 删除其他信息，保存为.txt文件
		3. 利用sva包COMBAT函数进行批次校正
	
	* 其他协变量校正[待补充]
		1. 利用线性回归lm()对其他协变量进行校正

------------需要思考的问题  
1.质量分析的每一项作用是什么?  
2.质量分析后什么样的结果才是符合预期的?什么样的结果证明芯片质量不合格?  
3.rma和mas5的区别是什么?  
4.normalization各种方法的应用条件是什么?  
5.在处理表达数据时先过滤探针还是先过滤样本?为什么?  
