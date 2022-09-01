# cellranger

## cellranger安装及参考基因组下载

1. 下载

    #参考基因组下载（human）

        curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    cellranger下载 https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest 的curl网页链接
    
        curl -o cellranger-6.1.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.1.2.tar.gz?Expires=1649559567&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjEuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NDk1NTk1Njd9fX1dfQ__&Signature=Zwk8l6ickwwGCQCAcsHziaoVtnS3MI0yjZMiaqV8UiyL6rW4vTjpyMHFLl04KDWpKpuMi~6D6RfHrJlZKTli---KBfC05b5u8mVcA28uDEZuSyHiOnEpI-cDdwtJ5SuGyjuNbvsUBxTCyJ~mMkrNxTivGC5XCpx6dj312qL4d4RImwEWmwMl0Nm7L8OcRGVLlujgODH51bwB03LCq1VoYIE-ECu7IhVCHlkXRzG9jLpnP98b4xf5x50WitToM4BcsW3kuQBk6w8AXjNLk3zLJwLbqKPxjAhWGZofYIKV4p-3C6G0wuJm-6XRYE0q25KrTYuDVm6Hhr~Yk46XM7Z1MA__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

2. 解压
3. 添加到环境（略）
4. 验证安装是否正确

    cellranger testrun --id=tiny

如果运行成功会显示 `Pipestance completed successfully!`

## cellranger定量环节
在目录下编写一个cellranger运行脚本，`run-cellranger.sh`，内容如下：

    db=/data1/tanhw/testsinglecell/demo/raw/refdata-gex-GRCh38-2020-A  
    ls $db  
    fq_dir=/data1/tanhw/testsinglecell/demo/raw
    cellranger count --id=$1 \
    --localcores=32 \
    --transcriptome=$db \
    --fastqs=$fq_dir \
    --sample=$1 \
    --nosecondary \
    --expect-cells=5000 

**注意修改**
  
>db后面为你所下载的参考基因组目录。    
fq_dir后面为你得原始fastq地址    
--localcores为指定最大使用线程数  
--nosecondary为不进行聚类分群分析，因为后续使用seurat分析  
--expect-cells为指定最大细胞数，根据项目决定

### 批量并行run-cellranger.sh进行比对定量

    cat SraAccList.txt |while read id;do (nohup  bash run-cellranger.sh $id 1>log-$id.txt 2>&1 & );done


### 批量顺序运行run-cellranger.sh进行比对定量

写一个脚本命名为muti.sh

    cat SRR_Acc_List.txt |while read id;do ( bash run-cellranger.sh $id  );done

后台执行脚本

程序是否跑完可以看log，执行完的结果在output，`web_summary.html`是质量分析，`filtered_feature_bc_matrix`中是下游分析的材料


