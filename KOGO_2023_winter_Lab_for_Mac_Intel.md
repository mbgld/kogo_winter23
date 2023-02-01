# ***ATTENTION***

아래 코드는 <br>
2023년 02월 02일 한국유전체학회 동계워크샵에서 실시하는<br>
Shotgun Metagenomics 분석 실습 수업을 위해 준비한 것입니다. <br>
---
강의자료 제작: <span style="color:blue">성균관대학교 삼성융합의과학원/강북삼성병원 **김한나**</span> (문의: 147942@hanmail.net) 
 
---

----

파란 글씨는 클릭이 가능하며, **오른쪽 클릭으로 새창으로 열기**로 여시면, 해당 웹페이지로 연결됩니다. Github 에서 코드를 Copy&Paste 하여 진행하며, 차후 다시보려면 html 파일을 사파리나 [크롬브라우저](https://www.google.com/intl/ko/chrome/)를 이용하기 바랍니다.


## Download 실습파일: KOGO_winter
- 실습에 사용할 PC의 OS 확인
    - Window (64bit)
    - #### Mac: Intel chip
    - Mac: M1 chip
    - Linux

- 실습 파일 (이메일로 사전 안내된 구글드라이브 링크로 다운로드)
- 압축파일 통째로 **바탕화면** 에 다운로드 후, 바탕화면에서 압축을 푸세요.
    	
**본 실습파일은 말 그대로 실습용으로서, 실제 연구나 논문작성 등에 사용할 수 없으며, 차후 연구자의 실제 연구 데이터로 본 메뉴얼을 참고하여 분석 시 활용할 수 있으나, 프로그램이나 DB의 업데이트에 따라 일부 명령어가 달라질 수 있습니다.**

---

## Contents ##


* [Raw data processing: QC/Fitering](#raw-data-processing-qcfitering)
    * [1. Raw data (fastq) quality 확인](#1-raw-data-fastq-quality-확인)
    * [2. Quality Filtering](#2-quality-filtering)
    * [3. Removal of unwanted reads (Human 시퀀스 제거)](#3-removal-of-unwanted-reads-human-시퀀스-제거))
   
* [Taxonomic & Functional Profiling](#taxonomic-functional-profiling)
    * [1. Taxonomic Profiling](#1-taxonomic-profiling)
    * [2. Functiona Profiling](#2-functional-profiling)
        
<br/>


# Raw data processing: QC/Fitering
### 1. Raw data (fastq) quality 확인
#### - Tool: [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 이용
---
- 작업폴더 (KOGO_workshop)로 이동
```sh
cd ~/Desktop/KOGO_winter
```
```sh
cd ./Inputs
```
```sh
ls
```
- 총 4명(M1, M2, F1, F2)의 paired-end fastq 파일 4세트
  - 각 1만 reads 시퀀스 데이터
  ```
  ├── Inputs
      ├── F1_10Krc_1.fastq.gz
      ├── F1_10Krc_2.fastq.gz
      ├── F2_10Krc_1.fastq.gz
      ├── F2_10Krc_2.fastq.gz
      ├── M1_10Krc_1.fastq.gz
      ├── M1_10Krc_2.fastq.gz
      ├── M2_10Krc_1.fastq.gz
      ├── M2_10Krc_2.fastq.gz  
  ```

- Fastq 파일 하나를 구경해봅시다
```sh
gzip -dc M2_10Krc_2.fastq.gz | less -S     # q 를 입력해서 빠져나옵니다.
```
- M2_10Krc_2.fastq.gz 1개 파일만 돌려보기
- ``cd ../`` 를 사용하여 한단계 위 폴더(KOGO_winter 폴더)위치로 이동

```sh
./Mac_tools/FastQC.app/Contents/MacOS/fastqc \
./Inputs/M2_10Krc_2.fastq.gz \
-o ./my_results/
```
- FastQC 결과파일
  - html 파일 1개에 모든 QC 결과가 저장됨
  - Basic Statistics: 시퀀스에 대한 기본적인 통계
  - Per base sequence quality: read 의 각 base position 이 가지는 평균 quality score
    - Warning: 어떤 base의 low quartile 이 <10 일때, 또는 median 이 <25 일때
    - Failure: low quartile <5 또는 median <20
  - Per title sequence quality
  - Per sequence quality scores: 각 reads 가 가지는 평균 quality score 의 분포도
    - Warning: mean quality <27
    - Failure: mean quality <20
  - Per base sequence contents: ATGC 구성차이가 >10% 일때 Warning, 20% 일때 Failure
    - 일부 구간에서 Failure 라면, library contamination (overrepresented sequence)등을 의심
    - 전 구간에서 일관된 bias: original library 시퀀스 문제이거나 시퀀싱동안 systematic 문제가 있는 것으로 의심
  - Per sequence GC content
  - Per base N content: N 으로 calling된 base 의 %
  - Sequence Length Distribution
  - Sequence Duplication Levels
  - Overrepresented sequences: 예상되는 것보다 많이 관측되는 특정 sequence
  - Adapter Content: 시퀀스에 adapter 시퀀스가 포함되어 있는지 확인
  
- FastQC 프로그램은 Java 기반으로 만들어져있으므로, 아이콘을 더블클릭해도 실행 가능하나, 여러파일을 동시에 실행하기 위해서는 터미널 이용이 효율적
- 8개의 input files 의 Fastqc 결과
  - Outputs 폴더 > 01_fastqc 폴더에 결과 저장
  - Finder로 확인하세요.
  - 실습에 사용한 fastq 파일은, 짧은시간의 실습진행을 위해, 작은 사이즈로 임의로 만들어진 파일이므로, Qaultiy 에 대한 판단은 불가능하다는 것을 참고하여 보시길 바랍니다.

```sh
├── Outputs
    ├── 01_fastqc
        ├── F1_10Krc_1_fastqc.html
        ├── F1_10Krc_1_fastqc.zip
        ├── F1_10Krc_2_fastqc.html
        ├── F1_10Krc_2_fastqc.zip
        ├── F2_10Krc_1_fastqc.html
        ├── F2_10Krc_1_fastqc.zip
        ├── F2_10Krc_2_fastqc.html
        ├── F2_10Krc_2_fastqc.zip
        ├── M1_10Krc_1_fastqc.html
        ├── M1_10Krc_1_fastqc.zip
        ├── M1_10Krc_2_fastqc.html
        ├── M1_10Krc_2_fastqc.zip
        ├── M2_10Krc_1_fastqc.html
        ├── M2_10Krc_1_fastqc.zip
        ├── M2_10Krc_2_fastqc.html
        ├── M2_10Krc_2_fastqc.zip
```
<br/>
  
### 2. Quality Filtering

#### - Tool: [Trimmomatic]("http://www.usadellab.org/cms/?page=trimmomatic") 이용
---

- 우리는 1. 단계에서 raw data 상태를 확인만 해본 것이지, 아직 어떤 QC 도 진행하지 않았습니다.
- 이제 Adapter 시퀀스 제거 및 low quality 시퀀스를 제거해보겠습니다.
- 작업폴더 (KOGO_workshop)로 이동
```sh
cd ~/Desktop/KOGO_winter
```

- Trimmomatic 을 이용하여 QC 진행

```bash
java -jar ./Mac_tools/Trimmomatic-0.39/trimmomatic-0.39.jar \
   PE \
   ./Inputs/M2_10Krc_1.fastq.gz \
   ./Inputs/M2_10Krc_2.fastq.gz \
   ./my_results/02_M2_TM_R1_paired.fastq \
   ./my_results/02_M2_TM_R1_unpaired.fastq \
   ./my_results/02_M2_TM_R2_paired.fastq \
   ./my_results/02_M2_TM_R2_unpaired.fastq \
   ILLUMINACLIP:./Mac_tools/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 \
   LEADING:3 \
   TRAILING:3 \
   SLIDINGWINDOW:4:20 \
   MINLEN:105
```
- output 파일 순서 유의
- ILLUMINACLIP: adapter 제거 옵션
- LEADING: 5'쪽 triming (Phred score <3 or N)
- TRAILING: 3'쪽 triming (Phred score <3 or N)
- SLIDINGWINDOW: 5'쪽부터 4 bases (:4)씩 스캐닝하면서, window 내 평균 quality <20 (:20)이면 clips 하는 trimming 방법
- MINLEN: read lenghs (150bp 경우) trimming 전후 70% 미만 길이(105bp)이면 제거


- Trimmomatic 에 대한 자세한 옵션 설명은 [Trimmomatic manual](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) 참고
- Trimmomatic 결과 파일

```
├── My_results
    ├── 02_M2_TM_R1_paired.fastq
    ├── 02_M2_TM_R1_unpaired.fastq
    ├── 02_M2_TM_R2_paired.fastq
    ├── 02_M2_TM_R2_unpaired.fastq
```

- 1차 QC가 끝난 파일을 FastQC를 이용해서 다시 살펴보기
```sh
./Mac_tools/FastQC.app/Contents/MacOS/fastqc \
./my_results/02_M2_TM_R2_paired.fastq \
-o ./my_results/
```
- FastQC 결과파일
```
├── My_results
    ├── 02_M2_TM_R2_paired_fastqc.html
    ├── 02_M2_TM_R2_paired_fastqc.zip
```
- raw data 의 FastQC 결과 파일과 비교해보세요. 무엇이 달라졌나요?

- 8개의 input files 의 Trimmomatic 결과
  - Outputs 폴더 > 02_trimmomatic 폴더에 결과 저장
  - Finder로 확인하세요.
```sh
├── Outputs
    ├── 01_fastqc
        ├── F1_TM_R1_paired.fastq
        ├── F1_TM_R1_unpaired.fastq
        ├── F1_TM_R2_paired.fastq
        ├── F1_TM_R2_unpaired.fastq
        ├── F2_TM_R1_paired.fastq
        ├── F2_TM_R1_unpaired.fastq
        ├── F2_TM_R2_paired.fastq
        ├── F2_TM_R2_unpaired.fastq
        ├── M1_TM_R1_paired.fastq
        ├── M1_TM_R1_unpaired.fastq
        ├── M1_TM_R2_paired.fastq
        ├── M1_TM_R2_unpaired.fastq
        ├── M2_TM_R1_paired.fastq
        ├── M2_TM_R1_unpaired.fastq
        ├── M2_TM_R2_paired.fastq
        ├── M2_TM_R2_unpaired.fastq
```
<br/>

### 3. Removal of unwanted reads (Human 시퀀스 제거)
####  - Tool: [Bowtie2]("https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml") 와 [samtools]("https://samtools.sourceforge.net") 이용
---

1). QC 가 끝난 cleaned fastq 파일을 human reference genome 에 mapping
```sh
./Mac_tools/bowtie2-2.5.0-macos-x86_64/bowtie2 \
 --very-sensitive-local \
 -p 2 \
 --seed 99 \
 -x ./DB/GRCh38_noalt_as/GRCh38_noalt_as \
 -1 ./my_results/02_M2_TM_R1_paired.fastq \
 -2 ./my_results/02_M2_TM_R2_paired.fastq \
 -S ./my_results/03_M2_mapped_and_unmapped.sam
```
- 결과파일: 03_M2_mapped_and_unmapped.sam
  - human reference genome 에 mapped 된 리드와, unmapped 리드들이 모두 이 하나의 파일에 들어있음
  
2). SAM 파일을 BAM 파일로 변환
```sh
./Mac_tools/samtools-1.16.1/samtools view \
 -@ 2 \
 -bS ./my_results/03_M2_mapped_and_unmapped.sam \
 > ./my_results/04_M2_mapped_and_unmapped.bam
```        

3). Human genome 에 mapped 되지 않은 (Metagenome으로 예상되는) unmapped reads 만 추출
- samtools flag 사용: [SAM flags]("https://broadinstitute.github.io/picard/explain-flags.html") 설명
  - -f12: Extract (-f) only alignments with both reads unmapped: "read unmapped"/"mate unmapped"
  - -F256: Do not extract (-F) alignments which are: "not primary alignment"
- R1 read 와 R2 read 모두 동시에 unmap 인 reads 만 추출
```sh
./Mac_tools/samtools-1.16.1/samtools view \
 -@ 2 \
 -b -f 12 -F 256 \
 ./my_results/04_M2_mapped_and_unmapped.bam \
 > ./my_results/05_M2.unmapped.bam
```

4). Sort the reads
```sh
./Mac_tools/samtools-1.16.1/samtools sort \
 -@ 2 \
 -n ./my_results/05_M2.unmapped.bam \
 -o ./my_results/06_M2.unmap.sorted.bam
```

5). Sorted reads 를 fastq 포멧으로 변환
```sh
./Mac_tools/samtools-1.16.1/samtools bam2fq \
 -@ 2 \
 ./my_results/06_M2.unmap.sorted.bam \
 > ./my_results/07_M2_metagenome.fastq
``` 

6). 최종 metagenome.fasq 파일 확인
```sh
less -S ./my_results/07_M2_metagenome.fastq
```

- PhiX 시퀀스 제거 역시 위 1)-5) 과정을 한번 더 똑같이 하면 됩니다.
- 8개의 input files (QCed fastq files) 의 Bowtie2 & Samtools 수행 결과
  - Outputs 폴더 > 03_filter_out_human 폴더에 결과 저장
  - Finder로 확인하세요.
  - low quality 시퀀스가 제거되고, host contaminant 도 제거된 메타지놈 시퀀스 파일입니다.
```sh
├── Outputs
    ├── 03_filter_out_human
       ├── F1.unmap.sorted.bam
       ├── F1_mapped_and_unmapped.bam
       ├── F1_mapped_and_unmapped.sam
       ├── F1_metagenome.fastq
       ├── F1.unmapped.bam
       ├── F2_mapped_and_unmapped.bam
       ├── F2_mapped_and_unmapped.sam
       ├── F2_metagenome.fastq
       ├── F2.unmap.sorted.bam
       ├── F2.unmapped.bam
       ├── M1_mapped_and_unmapped.bam
       ├── M1_mapped_and_unmapped.sam
       ├── M1_metagenome.fastq
       ├── M1.unmap.sorted.bam
       ├── M1.unmapped.bam
       ├── M2_mapped_and_unmapped.bam
       ├── M2_mapped_and_unmapped.sam
       ├── M2_metagenome.fastq
       ├── M2.unmap.sorted.bam
       ├── M2.unmapped.bam
```
<br/>

---
(옵션) [Kneaddata]("https://huttenhower.sph.harvard.edu/kneaddata/") 프로그램을 사용하여, Trimommatic 과 bowtie2 를 한꺼번에 돌릴 수도 있음
#### - ``본 실습에서는 진행하지 않겠습니다.`` 필요하신 분들만 워크샵 이후에 참고하세요!
###### - Kneaddata 설치
```sh
#pip install kneaddata   #실제 사용할때는 명령어 앞 #지우고 사용

```
###### - Kneaddata 실행
```sh
#kneaddata \
#-i ./Inputs/M1_10Krc_1.fastq.gz -i ./Inputs/M1_100Krc_2.fastq.gz \
#-db ./DB/GRCh38_noalt_as/GRCh38_noalt_as \
#-db ./DB/PhiX_bt/genome \
#-o ./Outputs/kneaddata_M1 \
#--fastqc ./Mac/FastQC.app/Contents/MacOS/fastqc \
#--run-trim-repetitive \
#--trimmomatic ./Trimmomatic-0.39 \
#--trimmomatic-options="SLIDINGWINDOW:4:20 MINLEN:105" \
#--sequencer-source TruSeq3 \
#--bypass-trf \
#--store-temp-output    #실제 사용할때는 명령어 앞 #지우고 사용
```
---

<br>


# Taxonomic & Functional Profiling

#### - Tool: [HUMAnN3]("https://huttenhower.sph.harvard.edu/humann/") 이용
- HUMAnN은 Functional profiling 을 하기위한 단계 중, Marker gene mapping 단계를 거치기때문에, MetaPhlAn (본 실습에서 사용하는 버전은 [MetaPhlAn4]("https://huttenhower.sph.harvard.edu/metaphlan/")) 결과가 중간 과정에서 생성됩니다.
- 한번의 run 으로 taxonomic & functional profiling 을 동시에 얻을 수 있다는 장점이 있습니다.
---

- input 파일에는, fastq, fastq.gz, fasta, fasta.gz 포멧 모두 가능
- 본인의 컴퓨터 사양(CPU)에 따라 가능한 core 수를 ``--threads`` 적용하면, 많은 core를 사용할 수록 당연히 분석시간이 짧아집니다.
- Translated search (DIAMOND)에 사용되는 DB 연구(자,환경)에 따라 여러가지를 선택할 수 있으며, 본 실습에서는 EC-filtered UniRef50 (0.3GB)의 작은 사이즈의 DB를 사용함
  - full UniRef90 (20.7GB. 실제 분석 시 추천)
  - EC-filtered UniRef90 (0.9GB)
  - full UniRef50 (6.9GB)
  - ``humann_databases --download uniref $DB명 $INSTALL_LOCATION`` 으로 다운로드 후,  ``protein-database $INSTALL_LOCATION`` 옵션을 추가하여 사용
  - ``--bypass-transtated-search`` 옵션을 사용하여 translated search를 skip 할 수도 있음
  

- HUMAnN 실행
```sh
humann \
--input ./Outputs/03_filter_out_human/M1_metagenome.fastq \
--output ./my_results/08_10K_M1_humann \
--threads 2 \
--protein-database ./DB/Humann/uniref
```
- 4명 샘플을 동시에 HUMAnN 실행 방법 ``(오늘 돌리면 집에 못갑니다 ㅠㅠ 실습시간에는 skip 하세요)``
  - 샘플 1개씩이 아닌, 폴더내 모든 샘플을 한꺼번에 돌리고 싶을때에는, filtered_fastq 파일들을 모두 한 폴더에 넣고, 해당 폴더명을 input 으로 지정할 수 있음
  - filtered_fastq 파일들을 QCed_fastq 폴더로 몰아 넣기
```sh
#cd ./Outputs/03_filter_out_human/   #실제 사용할때는 명령어 앞 #지우고 사용

#mkdir QCed_fastq

#mv *.fastq ./QCed_fastq

#humann \
#--input ./Outputs/03_filter_out_human/QCed_fastq \
#--output ./my_results/04_humann_results_all
```


### - HUMAaN Outputs
- my_results 폴더 > 08_10K_M1_humann 폴더를 Finder 에서 확인하세요.
- Outputs 폴더 > 04_taxaNfunc 폴더
  - 총 4명 샘플에 대한 HUMaN 결과파일
  - 이 결과는 10만 리드짜리 input files 을 이용하여, full UniRef90 DB 로 translated search된 결과입니다.
- Finder로 확인하세요.
```sh
├── Outputs
    ├── 04_taxNfunc
        ├── 100K_M1_humann   
            ├── *_genefamilies.tsv   # gene families (Uniref90) 결과 (단위: RPK, reads per kilobase)
            ├── *_pathabundance.tsv   # pathways (MetaCyc) 결과 
            ├── *_pathcoverage.tsv
            ├── *_humann_temp
                ├── *_bowtie2_alinged.sam  # bwotie2 로 align 한 모든 alignments
                ├── *_bowtie2_aligned.tsv  # sam 파일의 작은 버전
                ├── *_bowtie2_index.*.bt2  # chochophlan DB 로부터 만들어진 index 파일 (6개)
                ├── *_bowtie2_unaligned.fa  # unmapped reads -> translated search (DIAMOND) 에 input 파일로 사용
                ├── *_custom_chocophlan_database.ffn 
                ├── *_diamond_aligned.tsv  # DIAMOND (translated alignment step) 실행 후 aglignment 결과
                ├── *_diamond_aligned.fa   # DIAMOND 실행 후 unaligned reads
                ├── *_metaphlan_bowtie2.txt  # metaphlan 실행 후 bowtie2 결과 -> Metaphlan 을 다른 옵션으로 다시 돌리고 싶을 때 input 으로 사용 가능 
                ├── *_metaphlan_bugs_list.tsv # metaphlan 실행후 Taxonomy profiling 결과파일!* (단위: relative abundance, 0-100 range)
                ├── *.log  # HUMAnN 실행 로그파일: 각 단계의 분석 소요시간, unaligned reads 의 %, alignment statistics 등 확인 가능
        ├── 100K_M2_humann  
        ├── 100K_F1_humann        
        ├── 100K_F2_humann    
```
<br/>

## 1. Taxonomic Profiling

#### - MetaPhlAn 결과 확인
- 작업폴더 (KOGO_workshop)에 있는 지 확인
```sh
cd ~/Desktop/KOGO_winter
```

```sh
less -S ./Outputs/04_taxaNfunc/100K_F1_humann/F1_metagenome_humann_temp/F1_metagenome_metaphlan_bugs_list.tsv   # 확인했으면 q 를 누르고 빠져나오세요
```


#### - MetaPhlAn profile 합치기
```sh
merge_metaphlan_tables.py \
./Outputs/05_taxonomy/*_metagenome_metaphlan_bugs_list.tsv \
> ./my_results/09_merged_taxa_table.txt

```
```sh
less -S ./my_results/09_merged_taxa_table.txt
```
#### - Taxonomic profiles 로 Heatmap 그리기
- Species abundance 만 추출하기
```sh
grep -E "s__|clade" ./Outputs/05_taxonomy/merged_taxa_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/clade_name/male_female/g" \
| sed "s/_metagenome_metaphlan_bugs_list//g" \
> ./my_results/10_merged_taxa_table_species.txt
```
```sh
less -S ./my_results/10_merged_taxa_table_species.txt
```

- hclust2.py 실행
- species level 만 그려보기 (top 10)
```sh
python ./Mac_tools/hclust2.py \
-i ./Outputs/05_taxonomy/merged_taxa_table_species.txt \
-o ./my_results/11_metaphlan4_heatmap_species.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--log_scale \
--ftop 10 \
--flabel_size 10 --slabel_size 10 \
--max_flabel_len 100 --max_slabel_len 100 \
--dpi 300
```
- 모든 taxa level에 대해서 그려보기 (top 10)

```sh
python ./Mac_tools/hclust2.py \
-i ./Outputs/05_taxonomy/merged_taxa_table.txt \
-o ./my_results/12_metaphlan4_heatmap_alltaxa.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--log_scale \
--ftop 10 \
--flabel_size 10 --slabel_size 10 \
--max_flabel_len 100 --max_slabel_len 100 \
--dpi 300
```
- plotting 에 대한 기타 상세옵션은 python hclust2.py --help 에서 확인가능
- [MaAsLin]("https://huttenhower.sph.harvard.edu/maaslin2/maaslin") 같은 통계분석 툴을 이용하여 보다 정확한 통계분석 가능

* (옵션) SGB profiles 을 Genome Taxonomy database ([GTDB]("https://gtdb.ecogenomic.org")) taxonomy 로 바꾸기
```sh
sgb_to_gtdb_profile.py -i $metaphlan_output -o $GTDB_output
```


 <br/>
 
## 2. Functional Profiling
```sh
├── Outputs
    ├── 04_taxaNfunc
        ├── 100K_M1_humann   
            ├── *_genefamilies.tsv   # gene families (Uniref90) 결과 (단위: RPK, reads per kilobase)
            ├── *_pathabundance.tsv   # pathways (MetaCyc) 결과 (단위: 
            ├── *_pathcoverage.tsv
```
- 각 output 의 보다 상세한 설명은 [HUMAnN 튜토리얼]("https://github.com/biobakery/humann") 참고
- Gene famlies: Groups of evolutionarily-related protein-coding sequences that often perform similar functions.
- Pathway abundance: Abundances of the pathway's component reactions with each reaction's abundance computed as the sum over abundances of genes catalyzing the reaction.
  - Pathway abundance is proportional to the number of complete "copies" of the pathway in the community. Thus, for a simple linear pathway RXN1 + RXN2 + RXN3 + RXN4, if RXN1 is 10 times as abundant as RXNs 2-4, the pathway abundance will be driven by the abundances of RXNs 2-4.
  - Unlike gene abundance, a pathway's community-level abundance is not necessarily the sum of its stratified abundance values. For example, continuing with the simple linear pathway example introduced above, if the abundances of RXNs 1-4 are [5, 5, 10, 10] in Species_A and [10, 10, 5, 5] in Species_B, HUMAnN 3.0 would report that Species_A and Species_B each contribute 5 complete copies of the pathway. However, at the community level, the reaction totals are [15, 15, 15, 15], and thus HUMAnN 3.0 would report 15 complete copies.
- Pathway coverage: Alternative description of the presence (1) and absence (0) of pathways in a community, independent of their quantitative abundance. A confidence score to each reaction detected in the community (1= confidently detected pathway, 0=less confidently detected). 

<br/>

#### - Normalization of output files
- 작업폴더 (KOGO_workshop)로 현재위치 이동
```sh
cd ~/Desktop/KOGO_winter
```
- genefamilies 파일 결과 확인해보기
```sh
less -S ./Outputs/04_taxaNfunc/100K_M1_humann/M1_metagenome_genefamilies.tsv
```

```sh
humann_renorm_table \
--input ./Outputs/04_taxaNfunc/100K_M1_humann/M1_metagenome_genefamilies.tsv \
--output ./my_results/13_M1_genefamilies_relab.tsv \
--units relab
```
- renomalized genefamilies 파일 결과 확인해보기
```sh
less -S ./my_results/13_M1_genefamilies_relab.tsv
```

<br/>


#### - 여러명의 Individual 샘플 결과를 하나의 파일로 합치기
1. Gene families
```sh
humann_join_tables \
--input ./Outputs/06_humann_renorm/genefam_relab \
--output ./my_results/14_humann_genefamilies_relab_merged.tsv \
--file_name genefamilies_relab
```
```sh
less -S ./my_results/14_humann_genefamilies_relab_merged.tsv
```

2. Pathways coverage
```sh
humann_join_tables \
--input ./Outputs/06_humann_renorm/pathcov \
--output ./my_results/15_humann_pathcoverage_merged.tsv \
--file_name pathcoverage
```
```sh
less -S ./my_results/15_humann_pathcoverage_merged.tsv
```

3. Pathway abundance
```sh
humann_join_tables \
--input ./Outputs/06_humann_renorm/pathabun_relab \
--output ./my_results/16_humann_pathabundance_relab_merged.tsv \
--file_name pathabun_relab
```

```sh
less -S ./my_results/16_humann_pathabundance_relab_merged.tsv
```

<br/>

#### - File split: stratified and unstratified table 생성
- Gene families 나 Pathway abundance 에 contribute 하는 taxa 정보를 포함/제외 하는 명령어 
1. Gene families
```sh
humann_split_stratified_table \
--input ./Outputs/06_humann_renorm/humann_genefamilies_relab_merged.tsv \
--output ./my_results/17_genefam_split
```
2. Pathways coverage
```sh
humann_split_stratified_table \
--input ./Outputs/06_humann_renorm/humann_pathcoverage_merged.tsv \
--output ./my_results/18_pathcov_split
```
3. Pathway abundance
```sh
humann_split_stratified_table \
--input ./Outputs/06_humann_renorm/humann_pathabundance_relab_merged.tsv \
--output ./my_results/19_pathabun_split
```



- 이렇게 합쳐진 파일은 [MaAsLin]("https://huttenhower.sph.harvard.edu/maaslin2/maaslin") 같은 통계분석 툴의 input file 로 사용 가능


<br/>

### - Downstream analysis
- 관심있는 phenotype 과 metagenome 에 대한 association 분석을 하기위해서는, 마이크로바이옴데이터에 특화되어 있는 다양한 통계 툴 이용가능
  - [MaAsLin]("https://huttenhower.sph.harvard.edu/maaslin/")
  - [ANCOM]("https://github.com/FrederickHuangLin/ANCOM-Code-Archive")
  - [HAllA]("https://huttenhower.sph.harvard.edu/halla/")
  - [ALDEX]("https://github.com/ggloor/ALDEx2_dev")
- HUMAnN 에서도 Kruskal-Wallis analysis 정도의 간단한 분석은 가능함
- phenotype 데이터와 메타지놈 데이터를 통합한 pcl 파일을 만들어서 이용
```sh
less -S ./Outputs/06_humann_renorm/sex_pathabund.pcl
```
```sh
humann_associate --input ./Outputs/06_humann_renorm/sex_pathabund.pcl \
--last-metadatum sex \
--focal-metadatum sex \
--focal-type categorical \
--output ./my_results/20_sex_path_stats.txt \
--fdr 1
```
<br>

### - HUMAnN bar plot 그리기
```sh
humann_barplot \
--input ./Outputs/06_humann_renorm/sex_pathabund.pcl \
--last-metadata sex \
--focal-feature PWY-7221 \
--focal-metadata sex \
--sort none \
--exclude-unclassified \
--output ./my_results/21_humann_barplot_PWY-7221
```
- sort: none, sum, dominant, braycurtis, braycurtis_m, metadata 중 선택 가능
- plotting 에 대한 기타 상세옵션은 humann_barplot --help 에서 확인가능

</br></br>

## 수고하셨습니다! </br>
### 본 스크립트는 잘 보관하셔서 차후 Metagenomics 분석하실 때 요긴하게 사용하시길 바랍니다~:poop::bug::sweat_smile::grin::heart:

---
강의자료 제작: <span style="color:blue">성균관대학교 삼성융합의과학원/강북삼성병원 **김한나**</span> (문의: 147942@hanmail.net) 
 
---









