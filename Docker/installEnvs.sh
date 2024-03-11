conda create \
--file envs/bcftools.yml \
--file envs/bedtools.yml \
--file envs/biopet.yml \
--file envs/bwa.yml \
--file envs/cutadapt.yml \
--file envs/erko.yml \
--file envs/py36.yml \
--file envs/samtools.yml \
--file envs/vardict.yml \
--file envs/varscan.yml \
--file envs/vep.yml \
&& conda clean --all -f -y