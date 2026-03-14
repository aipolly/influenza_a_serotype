# influenza_a_serotype



We identified significant slowdowns when processing large sample sets. To address this, we have implemented the following optimizations to improve computational efficiency: 

1. **Pre-built `minimap2` Index:** We modified the workflow to build the `minimap2` index separately. This eliminates the overhead of repeatedly generating the index for every sample, significantly reducing runtime. 
2. **Removal of `seqkit stat`:** The `seqkit stat` step has been removed from the pipeline to streamline the process and reduce unnecessary I/O operations. 
3. **Optimized Read Classification (Single-Pass FASTQ Reading):** We refactored the reads classification method to ensure that FASTQ files are read only **once**. 
3. replace pysam with samtools in `minimap2_sr` function.

If quality-controlled reads and classified read outputs are not required for downstream analysis, we recommend enabling the performance mode parameters below:

​	`--qual` False

​	`--get-reads` False



This repository is for personal study and research purposes only. For full implementation details and original design, please refer to the original project.

original project url:  <https://github.com/mtisza1/influenza_a_serotype>



#### use conda/mamba

```bash
# create env
git clone https://github.com/aipolly/influenza_a_serotype.git
git checkout ref_analysis
mamba create -n iav_serotype -f environment/iav_serotype.yaml

# download database
cd influenza_a_serotype
wget https://zenodo.org/records/17354032/files/Influenza_Lite_DB.tar.gz
tar -zxf Influenza_Lite_DB.tar.gz

# run test
cd test_data
tar -zxf SRR28752446.tar.gz


mamba run -n iav_serotype iav_serotype \
 -r iav_test_data_v020/SRR28752446-1.10k.fastq iav_test_data_v020/SRR28752446-2.10k.fastq \
 -s sample \
 -o sample_out \
 --db ../liteDB_v1.1/ \
 -t 4 \
 --read_format short_paired
 
 mamba run -n iav_serotype  serotype_ref_analysis \
 -b sample_out/sample/sample_influenza_A.sorted.bam \
 -s sample \
 -o sample_out \
 --db ../liteDB_v1.1 \
 --cov-thresh 1

```



#### using pip

comfirm you already have samtools seqkit and fastp ...

```bash

pip install git+https://github.com/aipolly/influenza_a_serotype.git

```

