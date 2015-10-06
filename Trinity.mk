#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#	./Trinity.mk all CPU=16 SAMP=10


DIR := ${CURDIR}
CPU=16
IWORM_CPU=10
TRIM=2
RUN=trinity
READ1=left.fastq
READ2=right.fastq
BCODES=barcodes.fa
min_kmer_cov=1
KMER_SIZE=25
jelly_hash_size=10G
grid_node_max_memory=1G
MIN_LEN=200

TRINITY ?= $(shell which 'Trinity.mk')
TRINDIR := $(dir $(firstword $(TRINITY)))
PATH:=$(MAKEDIR):$(PATH)
JELLYFISH_DIR := $(TRINDIR)/trinity-plugins/jellyfish/bin/
TRIMMOMATIC_DIR := $(TRINDIR)/trinity-plugins/Trimmomatic/


all: mkdirs $(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa $(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100 \
	$(DIR)/$(RUN)_out_dir/$(RUN)_bwa_index.sa $(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam \
	$(DIR)/$(RUN)_out_dir/chrysalis/iworm_scaffolds.txt $(DIR)/$(RUN)_out_dir/chrysalis/GraphFromIwormFasta.out


mkdirs:
	mkdir -p $(DIR)/$(RUN)_out_dir
	mkdir -p $(DIR)/$(RUN)_out_dir/chrysalis

$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa: $(READ1) $(READ2)
	seqtk mergepe $(READ1) $(READ2) \
	| skewer -m pe -l $(KMER_SIZE) --quiet -Q $(TRIM) -t $(CPU) -x $(TRIMMOMATIC_DIR)/adapters/TruSeq3-PE.fa - -1 \
	| tee $(DIR)/$(RUN)_out_dir/both.fq \
	| $(JELLYFISH_DIR)/jellyfish count -t $(CPU) -m $(KMER_SIZE) -s $(jelly_hash_size) -o /dev/stdout /dev/stdin 2> /dev/null \
	| $(JELLYFISH_DIR)/jellyfish dump -L $(min_kmer_cov) /dev/stdin -o $(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa

$(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100:$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa
	cd $(DIR)/$(RUN)_out_dir/ && \
	$(TRINDIR)/Inchworm/bin/inchworm --kmers jellyfish.kmers.fa --run_inchworm -K $(KMER_SIZE) -L $(KMER_SIZE) --monitor 1 \
	--DS --keep_tmp_files --num_threads $(IWORM_CPU) --PARALLEL_IWORM  > $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa.tmp 2>/dev/null
	mv $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa.tmp $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa
	$(TRINDIR)/util/misc/fasta_filter_by_min_length.pl $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa 100 \
	> $(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100

$(DIR)/$(RUN)_out_dir/$(RUN)_bwa_index.sa :$(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100
	cd $(DIR)/$(RUN)_out_dir/chrysalis/ && \
	bwa index -p $(RUN)_bwa_index $(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100

$(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam:$(DIR)/$(RUN)_out_dir/$(RUN)_bwa_index.sa
	cd $(DIR)/$(RUN)_out_dir/chrysalis/ && \
	bwa mem -v 1 -p -t $(CPU) $(RUN)_bwa_index $(DIR)/$(RUN)_out_dir/both.fq \
	| awk -F '\t' 'BEGIN{OFS="\t";} /^@/ {print;next;} {$$1=sprintf("%s/%d",$$1,and($$2,64)==0?2:1); print;}' \
	| samtools view  -T . -bu - \
	| samtools sort -l 0 -O bam -T tmp -@ $(CPU) -m $(grid_node_max_memory) -o $(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam -

$(DIR)/$(RUN)_out_dir/chrysalis/iworm_scaffolds.txt:$(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam
	$(TRINDIR)/util/support_scripts/scaffold_iworm_contigs.pl \
	$(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam \
	$(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa > $(DIR)/$(RUN)_out_dir/chrysalis/iworm_scaffolds.txt

$(DIR)/$(RUN)_out_dir/chrysalis/GraphFromIwormFasta.out:$(DIR)/$(RUN)_out_dir/chrysalis/iworm_scaffolds.txt
	$(TRINDIR)/Chrysalis/GraphFromFasta -i $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa -r $(DIR)/$(RUN)_out_dir/both.fq \
	-min_contig_length $(MIN_LEN) -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 \
	-t $(CPU) -k 24 -kk 48  -scaffolding $(DIR)/$(RUN)_out_dir/chrysalis/iworm_scaffolds.txt  > $(DIR)/$(RUN)_out_dir/chrysalis/GraphFromIwormFasta.out