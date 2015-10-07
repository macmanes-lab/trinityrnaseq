#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#	./Trinity.mk all CPU=16 SAMP=10


DIR := ${CURDIR}
CPU=10
IWORM_CPU=10
TRIM=2
RUN=trinity
READ1=
READ2=
min_kmer_cov=1
KMER_SIZE=25
jelly_hash_size=10G
grid_node_max_memory=1G
MIN_LEN=200
max_mem_reads=50000000
TRINITY ?= $(shell which 'Trinity.mk')
TRINDIR := $(dir $(firstword $(TRINITY)))
PATH := $(MAKEDIR):$(PATH)
JELLYFISH_DIR := $(TRINDIR)/trinity-plugins/jellyfish/bin/
TRIMMOMATIC_DIR := $(TRINDIR)/trinity-plugins/Trimmomatic/
bflyHeapSpaceInit = 1G
bflyHeapSpaceMax = 4G
bflyGCThreads=2
max_frag_len = 1000
num_reads=100000

.SILENT:inchworm2 graph2 bundle2 read2comp2 sort2 FastaToDeBruijn partition qgraph bfly step2

step2: inchworm2 graph2 bundle2 read2comp2 sort2 FastaToDeBruijn partition qgraph bfly concatenate

inchworm2: $(RUN)/inchworm.K25.L25.DS.fa
graph2:$(RUN)/chrysalis/GraphFromIwormFasta.out
bundle2:$(RUN)/bundled_iworm_contigs.fasta
read2comp2:$(RUN)/readsToComponents.out
sort2:$(RUN)/chrysalis/readsToComponents.out.sort
FastaToDeBruijn:$(RUN)/chrysalis/bundled_iworm_contigs.fasta.deBruijn
partition:$(RUN)/chrysalis/component_base_listing.txt
qgraph:$(RUN)/chrysalis/Component_bins/Cbin0/c0.graph.out



$(RUN)/inchworm.K25.L25.DS.fa:
	mkdir -p $(RUN) && cd $(RUN) && \
	$(TRINDIR)/Inchworm/bin/inchworm --reads $(READ1) \
	--run_inchworm -K $(KMER_SIZE) -L $(KMER_SIZE) --monitor 1 \
	--DS --num_threads 1 --PARALLEL_IWORM  > $(RUN)/inchworm.K25.L25.DS.fa 2>/dev/null


$(RUN)/chrysalis/GraphFromIwormFasta.out:
	mkdir -p $(RUN)/chrysalis/ && \
	$(TRINDIR)/Chrysalis/GraphFromFasta -i $(RUN)/inchworm.K25.L25.DS.fa \
	-r $(READ1) -min_contig_length $(MIN_LEN) -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 \
	-t $(CPU) -k 24 -kk 48   \
	> $(RUN)/chrysalis/GraphFromIwormFasta.out


$(RUN)/bundled_iworm_contigs.fasta:
	$(TRINDIR)/Chrysalis/CreateIwormFastaBundle \
	-i $(RUN)/chrysalis/GraphFromIwormFasta.out \
	-o $(RUN)/chrysalis/bundled_iworm_contigs.fasta -min $(MIN_LEN) 2>/dev/null


$(RUN)/readsToComponents.out:
	$(TRINDIR)/Chrysalis/ReadsToTranscripts -i $(READ1) \
	-f $(RUN)/chrysalis/bundled_iworm_contigs.fasta \
	-o $(RUN)/chrysalis/readsToComponents.out \
	-t $(CPU) -max_mem_reads $(max_mem_reads)  2>/dev/null


$(RUN)/chrysalis/readsToComponents.out.sort:
	/usr/bin/sort --parallel=6 -T . -S 10G -k 1,1n \
	$(RUN)/chrysalis/readsToComponents.out \
	> $(RUN)/chrysalis/readsToComponents.out.sort 2>/dev/null

$(RUN)/chrysalis/bundled_iworm_contigs.fasta.deBruijn:
	$(TRINDIR)/Inchworm/bin//FastaToDeBruijn \
		--fasta $(RUN)/chrysalis/bundled_iworm_contigs.fasta \
		-K 24 --graph_per_record --threads $(CPU) \
		> $(RUN)/chrysalis/bundled_iworm_contigs.fasta.deBruijn

$(RUN)/chrysalis/component_base_listing.txt:
	$(TRINDIR)/util/support_scripts/partition_chrysalis_graphs_n_reads.pl \
		--deBruijns $(RUN)/chrysalis/bundled_iworm_contigs.fasta.deBruijn \
		--componentReads $(RUN)/chrysalis/readsToComponents.out.sort -N 1000 -L 200

$(RUN)/chrysalis/Component_bins/Cbin0/c0.graph.out:
	/share/trinityrnaseq/Chrysalis/QuantifyGraph -g $(RUN)/chrysalis/Component_bins/Cbin0/c0.graph.tmp  \
		-i $(RUN)/chrysalis/Component_bins/Cbin0//c0.reads.tmp \
		-o $(RUN)/chrysalis/Component_bins/Cbin0/c0.graph.out -max_reads 200000  -k 24

bfly:
	java -Xmx$(bflyHeapSpaceMax) -Xms$(bflyHeapSpaceInit) \
	-XX:ParallelGCThreads=$(bflyGCThreads) \
	-jar $(TRINDIR)/Butterfly/Butterfly.jar \
	-N $(num_reads) -L $(MIN_LEN) -F $(max_frag_len) \
	-C $(RUN)/chrysalis/Component_bins/Cbin0/c0.graph

concatenate:
	find read_partitions/  -name '*inity.fasta'  | $(TRINDIR)/util/support_scripts/partitioned_trinity_aggregator.pl TRINITY_DN > $(DIR)/$(RUN)_out_dir/Trinity.fasta


