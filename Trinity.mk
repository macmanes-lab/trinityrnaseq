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


all: mkdirs jellyfish inchworm index bwa iworm_scaffolds graph bundle read2comp sort list recursive
step2: inchworm2 graph2 bundle2 read2comp2 sort2 FastaToDeBruijn list concatenate

jellyfish:$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa
inchworm: $(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100
index:$(DIR)/$(RUN)_out_dir/chrysalis/$(RUN)_bwa_index.sa
bwa:$(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam
iworm_scaffolds:$(DIR)/$(RUN)_out_dir/chrysalis/iworm_scaffolds.txt
graph:$(DIR)/$(RUN)_out_dir/chrysalis/GraphFromIwormFasta.out
bundle:$(DIR)/$(RUN)_out_dir/chrysalis/bundled_iworm_contigs.fasta
read2comp:$(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out
sort:$(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort
list:$(DIR)/$(RUN)_out_dir/partitioned_reads.files.list
recursive:$(DIR)/$(RUN)_out_dir/recursive_trinity.cmds
concatenate:$(DIR)/$(RUN)_out_dir/Trinity.fasta
inchworm2: $(RUN)/inchworm.K25.L25.DS.fa
graph2:$(RUN)/chrysalis/GraphFromIwormFasta.out
bundle2:$(RUN)/bundled_iworm_contigs.fasta
read2comp2:$(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out
sort2:$(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out.sort
FastaToDeBruijn:$(DIR)/$(RUN)_out_dir/read_partitions/bundled_iworm_contigs.fasta.deBruijn


mkdirs:
	mkdir -p $(DIR)/$(RUN)_out_dir
	mkdir -p $(DIR)/$(RUN)_out_dir/chrysalis/
	mkdir -p $(DIR)/$(RUN)_out_dir/read_partitions


$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa: $(READ1) $(READ2)
	seqtk mergepe $(READ1) $(READ2) \
	| skewer -m pe -l $(KMER_SIZE) --quiet -Q $(TRIM) -t $(CPU) -x $(TRIMMOMATIC_DIR)/adapters/TruSeq3-PE.fa - -1 \
	| tee $(DIR)/$(RUN)_out_dir/both.fq \
	| $(JELLYFISH_DIR)/jellyfish count -t $(CPU) -m $(KMER_SIZE) -s $(jelly_hash_size) -o /dev/stdout /dev/stdin 2> /dev/null \
	| $(JELLYFISH_DIR)/jellyfish dump -L $(min_kmer_cov) /dev/stdin -o $(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa
	seqtk seq -A $(DIR)/$(RUN)_out_dir/both.fq > $(DIR)/$(RUN)_out_dir/both.fa &

$(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100:$(DIR)/$(RUN)_out_dir/jellyfish.kmers.fa
	cd $(DIR)/$(RUN)_out_dir/ && \
	$(TRINDIR)/Inchworm/bin/inchworm --kmers jellyfish.kmers.fa --run_inchworm -K $(KMER_SIZE) -L $(KMER_SIZE) --monitor 1 \
	--DS --keep_tmp_files --num_threads $(IWORM_CPU) --PARALLEL_IWORM  > $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa.tmp 2>/dev/null
	mv $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa.tmp $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa
	$(TRINDIR)/util/misc/fasta_filter_by_min_length.pl $(DIR)/$(RUN)_out_dir/inchworm.K25.L25.DS.fa 100 \
	> $(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100

$(DIR)/$(RUN)_out_dir/chrysalis/$(RUN)_bwa_index.sa:$(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100
	cd $(DIR)/$(RUN)_out_dir/chrysalis/ && \
	bwa index -p $(RUN)_bwa_index $(DIR)/$(RUN)_out_dir/chrysalis/inchworm.K25.L25.DS.fa.min100

$(DIR)/$(RUN)_out_dir/chrysalis/iworm.bowtie.nameSorted.bam:$(DIR)/$(RUN)_out_dir/chrysalis/$(RUN)_bwa_index.sa
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

$(DIR)/$(RUN)_out_dir/chrysalis/bundled_iworm_contigs.fasta:$(DIR)/$(RUN)_out_dir/chrysalis/GraphFromIwormFasta.out
	$(TRINDIR)/Chrysalis/CreateIwormFastaBundle -i $(DIR)/$(RUN)_out_dir/chrysalis/GraphFromIwormFasta.out -o $(DIR)/$(RUN)_out_dir/chrysalis/bundled_iworm_contigs.fasta -min $(MIN_LEN) 2>/dev/null

$(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out:$(DIR)/$(RUN)_out_dir/chrysalis/bundled_iworm_contigs.fasta
	$(TRINDIR)/Chrysalis/ReadsToTranscripts -i $(DIR)/$(RUN)_out_dir/both.fa -f $(DIR)/$(RUN)_out_dir/chrysalis/bundled_iworm_contigs.fasta -o $(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out \
	-t $(CPU) -max_mem_reads $(max_mem_reads)  2>/dev/null

$(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort:$(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out
	/usr/bin/sort --parallel=6 -T . -S 10G -k 1,1n $(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out > $(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort 2>/dev/null

$(DIR)/$(RUN)_out_dir/partitioned_reads.files.list:$(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort
	cat $(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort | sort -k1,1 -u | awk '{print "$(DIR)/$(RUN)_out_dir/read_partitions/c"$$1".trinity.reads.fa"}' > $(DIR)/$(RUN)_out_dir/partitioned_reads.files.list
	for i in $$(awk '{print $$1}' $(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort | sort -u); do grep -w ^$$i $(DIR)/$(RUN)_out_dir/chrysalis/readsToComponents.out.sort | awk '{print $$2 "\n" $$4}' \
	> $(DIR)/$(RUN)_out_dir/read_partitions/c$$i.trinity.reads.fa ; done

$(DIR)/$(RUN)_out_dir/recursive_trinity.cmds:$(DIR)/$(RUN)_out_dir/partitioned_reads.files.list
	$(TRINDIR)/util/support_scripts/write_partitioned_trinity_cmds.pl --reads_list_file $(DIR)/$(RUN)_out_dir/partitioned_reads.files.list --CPU 1 --max_memory 2G  --full_cleanup --seqType fq \
	--trinity_complete > $(DIR)/$(RUN)_out_dir/recursive_trinity.cmds

$(DIR)/$(RUN)_out_dir/Trinity.fasta:$(DIR)/$(RUN)_out_dir/recursive_trinity.cmds
	@echo '\n\n'
	@echo --------------------------------------------------------------------------------
	@echo ------------ Trinity Phase 2: Assembling Clusters of Reads ---------------------
	@echo --------------------------------------------------------------------------------
	@echo '\n\n'
	$(TRINDIR)/trinity-plugins/parafly/bin/ParaFly -c $(DIR)/$(RUN)_out_dir/recursive_trinity.cmds -CPU $(CPU) -v
	find read_partitions/  -name '*inity.fasta'  | $(TRINDIR)/util/support_scripts/partitioned_trinity_aggregator.pl TRINITY_DN > $(DIR)/$(RUN)_out_dir/Trinity.fasta


#/share/trinityrnaseq/Inchworm/bin//inchworm --reads single.fa --run_inchworm -K 25 -L 25 --monitor 1  --DS  --num_threads 1  --PARALLEL_IWORM  > /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/inchworm.K25.L25.DS.fa.tmp
#DONE!!!
$(RUN)/inchworm.K25.L25.DS.fa:
	mkdir -p $(RUN) && cd $(RUN) && \
	$(TRINDIR)/Inchworm/bin/inchworm --reads $(READ1) \
	--run_inchworm -K $(KMER_SIZE) -L $(KMER_SIZE) --monitor 1 \
	--DS --num_threads 1 --PARALLEL_IWORM  > $(RUN)/inchworm.K25.L25.DS.fa 2>/dev/null

#/share/trinityrnaseq/Chrysalis/GraphFromFasta -i /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/inchworm.K25.L25.DS.fa -r single.fa -min_contig_length 200 -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 -t 1 -k 24 -kk 48  > /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/GraphFromIwormFasta.out

#graph2
#DONE!
$(RUN)/chrysalis/GraphFromIwormFasta.out:
	$(TRINDIR)/Chrysalis/GraphFromFasta -i $(RUN)/inchworm.K25.L25.DS.fa \
	-r $(READ1) -min_contig_length $(MIN_LEN) -min_glue 2 -glue_factor 0.05 -min_iso_ratio 0.05 \
	-t $(CPU) -k 24 -kk 48   \
	> $(RUN)/chrysalis/GraphFromIwormFasta.out

#/share/trinityrnaseq/Chrysalis/CreateIwormFastaBundle -i /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/GraphFromIwormFasta.out -o /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/bundled_iworm_contigs.fasta -min 200

#bundle2
$(RUN)/bundled_iworm_contigs.fasta:
	$(TRINDIR)/Chrysalis/CreateIwormFastaBundle -i $(RUN)/GraphFromIwormFasta.out -o $(RUN)/bundled_iworm_contigs.fasta -min $(MIN_LEN) 2>/dev/null

#/share/trinityrnaseq/Chrysalis/ReadsToTranscripts -i single.fa -f /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/bundled_iworm_contigs.fasta -o /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/readsToComponents.out -t 1 -max_mem_reads 50000000

$(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out:
	$(TRINDIR)/Chrysalis/ReadsToTranscripts -i $(DIR)/$(RUN)_out_dir/both.fa -f $(DIR)/$(RUN)_out_dir/chrysalis/bundled_iworm_contigs.fasta -o $(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out \
	-t $(CPU) -max_mem_reads $(max_mem_reads)  2>/dev/null

#/usr/bin/sort --parallel=1 -T . -S 1G -k 1,1n /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/readsToComponents.out > /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/readsToComponents.out.sort

$(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out.sort:
	/usr/bin/sort --parallel=6 -T . -S 10G -k 1,1n $(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out > $(DIR)/$(RUN)_out_dir/read_partitions/readsToComponents.out.sort 2>/dev/null

#/share/trinityrnaseq/Inchworm/bin//FastaToDeBruijn --fasta /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/bundled_iworm_contigs.fasta -K 24 --graph_per_record --threads 1 > /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/bundled_iworm_contigs.fasta.deBruijn

#FastaToDeBruijn
$(DIR)/$(RUN)_out_dir/read_partitions/bundled_iworm_contigs.fasta.deBruijn:
	$(TRINDIR)/Inchworm/bin//FastaToDeBruijn --fasta $(DIR)/$(RUN)_out_dir/read_partitions/bundled_iworm_contigs.fasta \
	-K 24 --graph_per_record --threads $(CPU) > $(DIR)/$(RUN)_out_dir/read_partitions/bundled_iworm_contigs.fasta.deBruijn

#/share/trinityrnaseq/util/support_scripts/partition_chrysalis_graphs_n_reads.pl --deBruijns /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/bundled_iworm_contigs.fasta.deBruijn --componentReads /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/readsToComponents.out.sort -N 1000 -L 200

#/share/trinityrnaseq/trinity-plugins/parafly/bin/ParaFly -c /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/quantifyGraph_commands -CPU 1 -failed_cmds failed_quantify_graph_commands.50578.txt -shuffle

#/share/trinityrnaseq/trinity-plugins/parafly/bin/ParaFly -c /home/macmanes/trinityrnaseq/trinity_out_dir/read_partitions/Fb_0/CBin_0/c32.trinity.reads.fa.out/chrysalis/butterfly_commands -shuffle -CPU 1 -failed_cmds failed_butterfly_commands.50578.txt






















