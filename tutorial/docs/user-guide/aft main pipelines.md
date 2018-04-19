# After main pipelines

基本上说，其实是[main pipelines description](main pipelines#main-pipelines-description)的第3步开始的往后的流程。

由于流程的复杂，所以我同样的写了一个较为简单的主要入口文件，并且写成了一个可交互的脚本，所以会涉及到许多需要输入/确认的地方。

```shell
python aft_pipelines_analysis/aft_pipelines_main.py 1,2,3,4,5,6,7
```

唯一需要输入的参数即`1,2,3,4,5,6,7`，由于流程的复杂，所以用最简单的方法进行了模块化的操作，以便在分析到一半是需要重新单独执行某个模块。


# 入口文件结构

1. 配置参数，包括了文件目录结构以及存储地址
2. 上述各模块的函数化部分
3. `main`主要执行步骤与判断


# 具体模块的步骤
 1. 根据项目名称，到服务器上下载所有的csv和vcf，到指定的文件夹。
 2. 根据本地项目文件内的`setting文件`，执行针对somatic SNV的筛选pipelines
 3. 上传筛选后的文件，并**提醒需要**在服务器上执行**以下命令**，然后将`with_info.csv`下载到本地<Br>
  ```
  python Whole_pipelines/special_fun/add_per_info_into_csv_v2.py -i filtered.csv -o with_info.csv -tb XK-*T/XK-*T.recal_reads.bam -nb XK-*W/XK-*W.recal_reads.bam
  ```
 4. 通过`Whole_pipelines/aft_pipelines_analysis/csv2bed.py`批量的将刚刚下载的`with_info.csv`转成`bed`文件
 5. 利用`bed`文件，并结合`Whole_pipelines/aft_pipelines_analysis/extracted_pos_from_vcf.py`脚本，将多个vcf文件合并并提取出特定的位点，从而生成新的一个vcf文件，该文件集合了Tumor样本和Normal样本的这些位点的基本信息。
 6. 通过这个vcf,直接进行`pcgr`所指定的流程，最后可生成一个文件夹，内有基因组的网页报告。
