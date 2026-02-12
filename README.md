# influenza_a_serotype



We identified significant slowdowns when processing large sample sets. To address this, we have implemented the following optimizations to improve computational efficiency: 

1. **Pre-built `minimap2` Index:** We modified the workflow to build the `minimap2` index separately. This eliminates the overhead of repeatedly generating the index for every sample, significantly reducing runtime. 
2. **Removal of `seqkit stat`:** The `seqkit stat` step has been removed from the pipeline to streamline the process and reduce unnecessary I/O operations. 
3. **Optimized Read Classification (Single-Pass FASTQ Reading):** We refactored the reads classification method to ensure that FASTQ files are read only **once**. 
3. replace pysam with samtools

If quality-controlled reads and classified read outputs are not required for downstream analysis, we recommend enabling the performance mode parameters below:

​	`--qual` False

​	`--get-reads` False



This repository is for personal study and research purposes only. For full implementation details and original design, please refer to the original project.

original project url:  <https://github.com/mtisza1/influenza_a_serotype>, 
