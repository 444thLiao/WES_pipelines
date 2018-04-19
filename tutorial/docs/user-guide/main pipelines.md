# main pipelines description

 1. Download 数据
 2. 编辑合适的setting.py
 2. 完成germline、somatic流程
 3. 通过扩增试剂盒的**bed文件**进行统计流程中的**bam文件**的深度信息(.info)。
 4. 从深度信息(.info)统计得到评估文件(assessment file)，作为初步评估该次测序好坏的结果。
 5. 从深度信息(.info)计算得到用以绘制 coverage-depths的散点图的数据。
 6. Download `csv、vcf`结果文件。
 7. 本地进行自定义的筛选流程
 8. 将数量大大减少的输出结果(筛选后.csv)上传到服务器，利用服务器的内存+未下载的深度信息(.info)，使用脚本将其深度信息加入`csv文件`生成`with_info.csv`
 9. 重新将`with_info.csv`下载到本地
 10. 并将其转化成`bed文件`
 11. 根据`bed`文件从`vcf文件`中提取位点，并且需要合并多个vcf文件。
 12. 合并后使用`pcgr`的流程将其整理成网页报告
 13. 额外的分析要求(haplotype分析、拷贝数热图展现)
 14. 整理成常规的分析结构

## germline SNV Calling

基本的命令较为简单，只要使用以下命令即可开始germline calling的流程。

`sample_Id` 可以是本次样本的Id，最简的ID名称即可，如`S0707_09A_CHG030266-hunhewenku-XK-32T-AGGCAGAA_L007_R2.fastq.gz` 应该写成`XK-32T`，至于R1、R2会自己进行检测。<Br>
如果本次的样本过多，可以使用`auto_detect`，自动进行检测输入地址下的所有测序文件。更多的参数可见[main.py API](../../api#mainpy)
```
python Whole_pipelines/main.py -A germline_gatk4 -arg {sample_Id} -set setting.py
```

## Somatic SNV Calling

与germline SNV Calling的操作相似，所以直接使用以下命令开始somatic calling的流程。**值得注意的是，由于somatic和germline的流程有部分重合，请勿同时跑同样本的两个流程**

更多的参数可见[main.py API](../../api#mainpy)
```
python Whole_pipelines/main.py -A somatci_gatk4 -arg {sample_Id} -set setting.py
```
