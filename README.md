# DCS第二组-Presentation
![MasterHead](pic/banner(2).png)
[![Typing SVG](https://readme-typing-svg.herokuapp.com?font=Fira+Code&pause=1000&width=435&lines=%F0%9F%98%8DDCS%E7%AC%AC%E4%BA%8C%E7%BB%84%E5%85%A8%E4%BD%93%E6%88%90%E5%91%98%E7%A5%9D%E6%82%A8%E8%82%A0%E9%81%93%E5%81%A5%E5%BA%B7%EF%BC%8C%E6%AF%8F%E5%A4%A9%E5%BC%80%E5%BF%83~)](https://git.io/typing-svg)


组长：李乔凡

组员：于子喧、周驰、达圣文

# Contents
- [🌈参考论文网址](#参考论文网址)
	- [达圣文部分](#达圣文部分)
	- [李乔凡部分](#李乔凡部分)
	- [周驰部分](#周驰部分)
	- [于子喧部分](#于子喧部分)

- [🌈宏基因组物种比对流程-详解&code](#流程详解)

- [🌈常用网站&公众号推荐](#常用网站-微信公众号推荐)
	- [常用网站](#网站)
	- [公众号](#公众号)

# 🌈参考论文网址
## 达圣文部分
1.人类基因组计划（HGP）：https://www.nature.com/nature/volumes/409/issues/6822

2.人类微生物组计划（HMP）：https://www.hmpdacc.org/

3.微生物组：https://doi.org/10.1016/j.cell.2018.02.044

4.人体肠道宏基因组计划（MetaHIT）：http://www.metahit.eu

## 李乔凡部分
1.人类微生物组研究指南：https://doi.org/10.1097/CM9.0000000000000871

2.孟德尔随机化分析验证因果关系：https://doi.org/10.1038/s41588-021-00968-y

3.肠道微生物组的长期遗传稳定性和个体特异性：https://doi.org/10.1016/j.cell.2021.03.024

4.肠道微生物组中的细菌 SNP 与宿主的 BMI 相关：https://doi.org/10.1038/s41591-023-02599-8

5.宏基因组基础流程：https://mp.weixin.qq.com/s/xHe1FHLm3n0Vkxz0nNbXvQ

6.微生物-疾病关联研究范式（MWAS）：
https://doi.org/10.1038/nature11450
https://doi.org/10.1038/nrmicro.2016.83

7.物种-代谢关联分析：
https://doi.org/10.1038/s41591-019-0458-7
https://doi.org/10.1038/s41586-022-04567-7

8.SV&SNP~代谢关联分析策略：
https://10.1016/j.cell.2021.03.024
https://10.103a8/s41591-023-02599-8

## 周驰部分
1.宏基因组的应用.临床：
https://genome.cshlp.org/content/29/5/831

https://www.pnas.org/doi/full/10.1073/pnas.1809700115

2.宏基因组的应用.环境：
https://www.nature.com/articles/s41586-024-07891-2

## 于子喧部分
1.鼻腔微生物关联：https://www.nature.com/articles/s42003-024-05822-5

2.口腔微生物关联：https://www.nature.com/articles/s41421-021-00356-0

3.皮肤微生物关联：https://onlinelibrary.wiley.com/doi/10.1002/advs.202300050

4.阴道微生物关联：https://www.nature.com/articles/s41467-024-52102-1

5.肠道微生物关联：https://www.nature.com/articles/s41588-021-00968-y

6.肠道微生物与代谢/表型相关：https://www.nature.com/articles/s41467-017-00900-1

7.描述性/关联性→因果性/机制性：https://www.nature.com/articles/s41588-021-00968-y

8.深入到基因层面研究：https://www.nature.com/articles/s41591-023-02599-8

# 🌈流程详解
样本👉提取DNA👉宏基因组测序👉数据预处理👉质控、比对（基于reads）👉组装/拼接分析👉物种组成、功能分析。

强烈建议参考我们组开发的metapi，包含了宏基因组的所有常规流程的软件，无需一个一个安装：https://github.com/ohmeta/metapi

下文流程基于我们组开发的metaprof流程：https://github.com/weiting-liang/metaprof

软件支持：

fastp：https://github.com/OpenGene/fastp

bowtie2：: https://github.com/BenLangmead/bowtie2

MetaPhlAn4：https://github.com/biobakery/MetaPhlAn
## environment
  - biopython>=1.76
  - bowtie2>=2.3.5.1
  - fastp>=0.20.1
  - metaphlan>=4.1.0
  - numpy>=1.18.4
  - pandas>=1.0.3
  - pigz>=2.3.4
  - samtools>=1.9
  - seqkit>=0.12.1
  - snakemake>=5.14.0
  - kraken2>=2.1.1
  - bracken>=2.5

## install

```
git clone https://github.com/weiting-liang/metaprof.git

cd metaprof
conda env create -n metaprof -f ./rules/env.yaml
conda activate metaprof
```


#database prepare  
human reference
https://github.com/marbl/CHM13
```
mkdir /path/database
cd database && mkdir humanhost && cd humanhost 
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip chm13v2.0.fa.gz -d
mkdir bowtie2_index && cd bowtie2_index
bowtie2-build ../chm13v2.0.fa chm13v2
```

metaphlan3
```
cd /path/database
mkdir metaphlan && cd metaphlan
metaphlan --install --index mpa_vOct22_CHOCOPhlAnSGB_202212 --bowtie2db metaphlan_database
```

kraken2
```
cd /path/database
mkdir kraken2 && cd kraken_pub
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20210517.tar.gz
tar -xzvf k2_standard_20210517.tar.gz -C ./k2_standard_20210517
```

## run
- change the `samples.txt` to adapt to your data: separator should be "\t": id^Ifq1^Ifq2$
- custom the `config.yaml` database's path and parameters

```
#dry run
snakemake --snakefile rules/profile.smk -n
#test
snakemake --snakefile rules/profile.smk --core 16 2> smk.log &
#cluster: custom the cluster.yaml
nohup sh snakemake.sh &
```

## output

2.results/  
  filter_summary.txt  
  metaphlan4.profile.merge.txt 
  bracken.merged.abundance.profile.*.tsv  

#assay:  
1.assay  
  01.trimming/  
  02.rmhost/  
  03.profile/  
  benchmarks/   #check the cpu's time and max_vms to optimize the cluster's parameters  
  cluster_logs/   
  logs/         #find programs' errors  


## references
https://github.com/ohmeta/metapi/blob/master/metapi/Snakefile

https://github.com/liu930724/meta_profile

## pipeline WDL
见文件[README1](metaphlan4_Pu_1.3.0.wdl)

```
如有需要使用，请联系我，我会向你提供docker
```
# 🌈常用网站-微信公众号推荐
## 网站
组学原始数据归档（GSA）：http://gsa.big.ac.cn

Qiita：https://qiita.ucsd.edu

MGnify：https://www.ebi.ac.uk/metagenomics

gcMeta：https://gcmeta.wdcm.org

R Markdown：https://rmarkdown.rstudio.com

R Graph Gallery：https://www.r-graph.gallery.com

SangerBox绘图网站：http://sangerbox.com/home.html

HiPlot绘图网站：https://hiplot.cn/basic

BIC绘图网站：https://www.bic.ac.cn/BIC/#/


## 公众号

生信益站    |    宏基因组    |    iNature    |    生信技能树    |    生信通
<p float="left">
  <img src="pic/生信益站.png" width="150" />
  <img src="pic/宏基因组.png" width="150" /> 
  <img src="pic/iNature.png" width="150" />
  <img src="pic/生信技能树.png" width="150" />
  <img src="pic/生信通.png" width="150" />
</p>

## 🌈致谢

感谢我们敬爱的刘小敏老师、张涛老师、肖亮老师、杨焕明老师。

感谢周驰、达圣文、于子喧同学的共同努力。

感谢梁卫婷师姐提供的metaprof流程。

特别感谢亲爱的于子喧同学，My love。

2024.10.23
