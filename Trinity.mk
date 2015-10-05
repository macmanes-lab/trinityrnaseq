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


TRINITY ?= $(shell which 'Trinity.mk')
TRINDIR := $(dir $(firstword $(TRINITY)))
PATH:=$(MAKEDIR):$(PATH)
JELLYFISH_DIR := $(TRINDIR)/trinity-plugins/jellyfish/bin/
TRIMMOMATIC_DIR := $(TRINDIR)/trinity-plugins/Trimmomatic/


all: mkdirs $(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa iworm


mkdirs:
	mkdir $(DIR)/$(RUN)_out_dir
	mkdir $(DIR)/$(RUN)_out_dir/chrysalis

$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa: $(READ1) $(READ2)
	seqtk mergepe $(READ1) $(READ2) \
	| skewer -m pe -l 25 --quiet -Q 5 -t 12 -x $(TRIMMOMATIC_DIR)/adapters/TruSeq3-PE.fa - -1 \
	| tee both.fq \
	| $(JELLYFISH_DIR)/jellyfish count -t $(CPU) -m $(KMER_SIZE) -s $(jelly_hash_size) -o /dev/stdout /dev/stdin 2> /dev/null \
	| $(JELLYFISH_DIR)/jellyfish dump -L $(min_kmer_cov) /dev/stdin -o $(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa

iworm:$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa
	cd $(DIR)/$(RUN)_out_dir/ && \
	$(TRINDIR)/Inchworm/bin/inchworm --kmers jellyfish.kmers.fa --run_inchworm -K $(KMER_SIZE) -L 25 --monitor 1 \
	--DS --keep_tmp_files --num_threads $(IWORM_CPU) --PARALLEL_IWORM  > $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa.tmp 2>/dev/null