整条pipeline包括前半部分直接在服务器或集群上用CFEA数据库的WGBS/RRBS流程跑，最后得到的bed文件在服务器上导入RnBeads（R包）跑完默认R包的全部流程，最最后会得到一个启动子，基因，tilling，cpg岛的差异情况。（后面其实还可以设置阈值和补充其他必要分析得到最终的候选biomarker，这里暂时不用）

# Note
* CFEA pipeline脚本里用的参考基因组是hg19，我们用的是hg38。
* 每次关闭R后rnbeads的参数才会重置为默认。
* 可能是环境问题hub服务器上偶尔报错，优先使用cnode服务器，或者集群（没有账号可以和我说一下）。
* 使用cnode和hub时注意在提交命令前htop命令看一下服务器负载情况，调整脚本的线程数等。
        
注意检查PATH和软件版本
# DATA & SCRIPTS
## 1.数据表格：https://docs.qq.com/sheet/DSER0bmF1UkxSc1hF?tab=BB08J2,详细原始文件位置，大小，类别等信息都在里面
> 我已经做完了第一个数据集GSE149438里胃癌GC和正常人Normal的测序数据，因为第一遍主要是配置环境和安装软件R包，调整和理解脚本和参数，认识软件大致输入输出和原理（大概），所以我在这个repo里加了个csv文件，包括5个正常人和5个胃癌，这里建议大家先尝试重现下小部分数据，避免浪费计算资源和时间（我前期用全部数据跑了很长时间。）

## 2.CFEApipeline：https://github.com/lemonsky123/CFEA-pipeline
> 这个pipeline整体脚本比较完善，但是有一些输入输出文件目录的框架需要注意一下，详见具体脚本描述，否则很出现很多报错。
> 还有一点就是最开始他让用conda配置python2.7的环境，实际上其脚本都是用更高等级的python写的，很多语法python2.7并不支持，下面列出了我在服务器上为这个流程测试成功的相关软件的版本信息。有的是系统自带的可以不用额外安装，有的不一定，总之坑很多。
|name|version|note|
|perl|v5.26.2|系统自带不需要额外安装|
|python|3.6.10|需要在anaconda里配置对应版本的Python环境|
|fastqc|v0.11.9| |
|multiqc|version 1.9|在cnode上是正常的，但在hub上因为语系变量的问题而无法运行，需要暂时export统一语系变量，长期更改需要加入.bashrc文件(即重新登录仍然生效)|
|bismark|v0.22.3| |
|bowtie2|gccversion 7.5.0| |
|samtools|1.7|本身在下载时库选择错误，需要在下载后添加软连接，否则后面的几个步骤报错|
|macs2|2.2.7.1| |	
|picard|NA| |	
|deeptools|3.4.3| |

## 3.RnBeads，强大且完备的R包，可以出来sequencing和array的数据，第二版本去年刚发表，用的人好像很少？参见官网：https://rnbeads.org/index.html
> R包安装和依附包安装是重点，我安装了近一周时间，多个版本的R上测试过，建议使用conda配置一个新的R环境，优先使用能自动并行的微软的MRO版本。（其实cnode server上配置的就是比较完备的MRO版本，还已经配置了rstudio-server，可以通过浏览器端口访问，实在自己安装不了也可以使用系统自带版本，比如我习惯在浏览器访问，用的就是系统的MRO）

# CFEA 
前面预处理可以不分GC和NC文件夹
注意：
* Bismark的index目录为index的上层目录
* Bowtie2的路径不要加上bowtie2
##01mapping bam/sam
Option1(advised): using available script
Help:
python mapping_WGBS.py -i [path fastq] -p [numer of processes] -m 
        [mapping result dir] -l [log dir] -b [path to put bam files] 
        -n [project name for multiQC] -t [bowtie2 path] -d [bismark index]
        -c [path to put trimmed fastq files]
Cmd:
python ../CFEA-pipeline/RRBS_WGBS/mapping_WGBS.py -i rawdata/  -m output/mapping/  -l output/log/mapping/ -b output/bam/  -n GSE149438  -d ref/genome/homo/hg38/  -c output/clean2 -p 14 -t ~/anaconda3/envs/CFEA/bin/
* 消耗时间较长，46对共400G的raw reads在hub上按-p 5(5个一起跑)整整跑了近两天，应该是脚本里面没有确定bismark内部的多线程处理参数。最后仅mapping的记录用于fastqc和multiqc（应该是甲基化测序数据处理的特色，直接trim后仍然是碱基不平衡的状态，可能不适合qc）。最后整体单一位点比对率均>60%
* bug: 遇到了samtools的附属包不存在或者不是最新版本的错误，导致samtools无法使用，以上脚本产生的bam为空, 后来发现是bismark与CFEA的版本不同导致新版本附带了一个旧版的samtools，该samtools在CFEA环境中优先识别，但无法打开，所以无法产生bam文件。https://github.com/bioconda/bioconda-recipes/issues/12100，已经通过软连接解决两个失败的重新跑一遍后正常
* 最好能大致阅读一下bismark的说明书
* 使用的基因组是lulab的贡献基因组hg38，位置在lulab intranet上有

Option2: MANUAL
for ncname in `ls -1 output/clean/NC/*.fq | cut -d/ -f4 | cut -d. -f1 | sort | uniq`; do bismark --parallel 6 -o output/sam/NC/ --temp_dir tmp/NC/   ref/genome/homo/hg38/  -1 output/clean/NC/$ncname.sra_1_val_1.fq -2 output/clean/NC/$ncname.sra_2_val_2.fq ; done

##02extract methylation coverage site

python ../CFEA-pipeline/RRBS_WGBS/extract_methylation_coverage.py -i output/bam/gc/ -c output/extract/gc/ -l output/log/extract/gc/ -p 5

产生.bismark.cov.gz符合RnBeads的测序数据对输入格式的要求，应该可以直接输入，后面RnBeads会自行merge
.bedGraph.gz仅保留了染色体，区段，甲基化比例三列，其中区段从1-based变成0-based
这里的结果其实已经可以直接用于rnbeads了，后面的##03merge minus/plus strand cpg和###04run merge其实###run merge不用接着跑了，得到的bigwig等其实只有数据库可视化会用

# RnBeads
> 具体安装和使用参见官网https://rnbeads.org/index.html

##导入RnBeads
rnb.run.analysis(dir.reports="/BioII/lulab_b/baopengfei/RnBeads/proj1/results/report",data.source=c("/BioII/lulab_b/baopengfei/RnBeads/proj1/data","/BioII/lulab_b/baopengfei/RnBeads/proj1/data/sample.csv","barcode"),data.type="bed.dir")
注意sample.csv包括barcode（文件名或绝对路径）和source_name（SRAmeta table默认的分组依据），且与bed文件在相同路径，内容如下：

barcode;Run;Age;alternate_ID;Assay Type;source_name
SRR11615795.bismark.cov;SRR11615795;47;Normal_40;Bisulfite-Seq;Normal
SRR11615796.bismark.cov;SRR11615796;57;GC_6;Bisulfite-Seq;GC
...

CMD:

rm(list = ls())

library(RnBeads)

设置rnb.options(assembly = "hg38",filtering.sex.chromosomes.removal = TRUE, differential.enrichment.lola = TRUE, differential.enrichment.go = TRUE, import.table.separator = ";", import.bed.style = "bismarkCov")

设置并行处理cpu核数（和MRO自动并行的关系比较复杂，我没有设置，之前尝试设置了好像也没有用？）：parallel.setup(6)

设置大数据可以使用硬盘补充RAM: rnb.options(disk.dump.big.matrices=TRUE, disk.dump.bigff=TRUE)，好像需要提前安装simpleCache包，rnbeads没有附带。  

rnb.run.analysis(dir.reports="/BioII/lulab_b/baopengfei/RnBeads/proj1/results/report_cnode_R35",data.source=c("/BioII/lulab_b/baopengfei/RnBeads/proj1/data","/BioII/lulab_b/baopengfei/RnBeads/proj1/data/sample.csv","barcode"),data.type="bed.dir")
