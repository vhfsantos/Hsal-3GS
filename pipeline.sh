#!/bin/bash
# Pipeline para montagem do genoma de H. salinarum NRC-1 usando reads longas
# ---------------------------------------------------------------------------
# Estudos da instabilidade genômica em Halobacterium salinarum NRC-1 via sequenciamento de long-reads.,
# Vinícius H F Santos; Tie Koide

helpFunction()
{
   echo ""
   echo "
╔════════════════════════════════════════╗
║   Pipeline para montagem do genoma de  ║
║ H. salinarum NRC-1 usando reads longas ║
╚════════════════════════════════════════╝
         Vinícius H F Santos"
   echo ""
   echo ""
   echo ""
   echo "Uso: $0 -f5 <fast5_dir> -t <n_threads>"
   echo -e "\t-f5 Caminho para o diretório dos arquivos fast5"
   echo -e "\t-t Número de threads"
   exit 1
}

while getopts "f5:t:" opt
do
   case "$opt" in
      f5 ) fast5_dir="$OPTARG" ;;
      t ) threads="$OPTARG" ;;
      ? ) helpFunction ;;
   esac
done

if [ -z "$fast5_dir" ]
then
   echo "Erro: Caminho para o diretório dos arquivos fast5 não definido. Saindo...";
   helpFunction
fi

if [ -z "$threads" ]
then
   echo "Erro: Numero de threads não definido. Saindo...";
   helpFunction
fi

# 0) Especificando os diretórios, número de threads e barcodes

barcodes=`echo barcode{01..06}`
albacore_dir="albacore_output"
deepbinner_dir="deepbinner"
porechop_dir="porechop"
canu_dir="canu"
nanopolishmakerange="/opt/nanopolish/scripts/nanopolish_makerange.py"
nanopolish_dir="nanopolish_output"
threadsparallel=4
threadsnanopolish=4

# ***

if [[ ! -d $nanopolish_dir ]] ; then mkdir $nanopolish_dir ; fi
if [[ ! -d $canu_dir ]] ; then mkdir $canu_dir ; fi
if [[ ! -d $porechop_dir ]] ; then mkdir $porechop_dir ; fi
if [[ ! -d $deepbinner_dir ]] ; then mkdir $deepbinner_dir ; fi
if [[ ! -d $albacore_dir ]] ; then mkdir $albacore_dir ; fi

# 1) Basecalling com Albacore
read_fast5_basecaller.py -r \
        -i $fast5_dir \
        -t $threads \
        -f FLO-MIN106 \
        -k SQK-LSK108 \
        -o fastq \
        -q 0 \
        -n 0 \
        -s $albacore_dir \
        --disable_filtering

# 2) Demultiplex com DeepBinner

#### 2.1) Concatenando os arquivos *fastq gerados pelo albacore

find ${albacore_dir}/workspace -name "*.fastq" -exec cat {} + > ${deepbinner_dir}/all_basecalled_reads.fastq

#### 2.2) Classificando as reads (demultiplex)
deepbinner classify -s /opt/Deepbinner/models/EXP-NBD103_read_starts \
                    -e /opt/Deepbinner/models/EXP-NBD103_read_ends \
		    --intra_op_parallelism_threads $threads \
		    --omp_num_threads $threads \
		    --verbose \
		    $fast5_dir > ${deepbinner_dir}/classifications 2> ${deepbinner_dir}/deepbinner_classify.err

#### 2.3) Binning
deepbinner bin --classes ${deepbinner_dir}/classifications \
	       --reads ${deepbinner_dir}/all_basecalled_reads.fastq \
	       --out_dir ${deepbinner_dir} > ${deepbinner_dir}/deepbinner_bin.log 2>&1 

# 3) Trimando adaptadores com Porechop
for i in $barcodes ; do
	porechop -t $threads \
		 --format fastq.gz \
		 -i ${deepbinner_dir}/${i}.fastq.gz \
		 -o ${porechop_dir}/${i}.fastq.gz > ${porechop_dir}/${i}.log 2>&1	
done

# 4) Montando o genoma com Canu
for i in $barcodes ; do
	canu -p $i \
		-d ${canu_dir}/$i \
		genomeSize=2.57m \
		-nanopore-raw \
		${porechop_dir}/$i".fastq.gz" > ${canu_dir}/${i}_canu_assembly.log 2>&1
done

# 5) Gerando consensus com Nanopolish
for i in $barcodes; do
	# 5.1) Indexando reads
	nanopolish index -d $fast5_dir $fastq_dir/${i}.fastq.gz \
			 -s sequencing_summary.txt \
			 > $nanopolish_dir/${i}_nanopolishIndex.log \
			 2> $nanopolish_dir/${i}_nanopolishIndex.err

	# 5.2) Alinhando e ordenando
        minimap2 -ax map-ont \
		 -t $threads \
                 ${canu_dir}/${i}/${i}.unitigs.fasta \
                 ${porchop_dir}/${i}.fastq.gz | \
        samtools sort -@ $threads -o ${nanopolish_dir}/${i}.bam > $nanopolish_dir/${i}_minimap2.log 2> $nanopolish_dir/${i}_minimap2.err

	# 5.3) Indexando arquivos `bam`
	samtools index $nanopolish_dir/${i}.bam

	# 5.4) Gerando `vcf` consensus
	python $nanopolishmakerange ${canu_dir}/${i}/${i}.unitigs.fasta |\
	parallel --results $nanopolishdir/${i}_nanopolish_results -P $threadsparallel \
        nanopolish variants --consensus -o $nanopolish_dir/${i}_{1}_nanopolish.vcf \
                                        -t $threadsnanopolish \
					-q dcm,dam \
                                        -r ${porechop_dir}/${i}.fastq.gz \
                                        -b $nanopolish_dir/${i}.bam \
                                        -g ${canu_dir}/${i}/${i}.unitigs.fasta \
                                        -w {1}

	# 5.5) Convertendo `vcf` para `fasta`
        nanopolish vcf2fasta -g ${canu_dir}/${i}/${i}.unitigs.fasta \
                             $nanopolish_dir/${i}_*_nanopolish.vcf > $nanopolish_dir/${i}_polished.fasta 2> $nanopolish_dir/${i}_vcf2fasta.err
done

