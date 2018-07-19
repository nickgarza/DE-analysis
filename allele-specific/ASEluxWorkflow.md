# ASElux Workflow

ASElux is a novel allelic count aligner from creator Zong Miao of UCLA. It is currently the fastest allele specific aligner with RNA-Seq data. It is written in C++ and requires a the GCC compiler to be loaded.

The pipeline for the package creates a static index based off of a reference genome during the building stage as well as a dynamic index specific to the SNP file in reference to a passed .fasta or .fastq file

## Acquire ASElux from from GitHub

 To install, clone the repo from: <https://github.com/abl0719/ASElux.git> into a directory you wish to work in using CMD ``` git clone https://github.com/abl0719/ASElux.git .```

## Sample Run 
A sample can be run by going into the demo folder and running the ```script.sh``` file which should appear as:

```{bash}
ASElux="../ASElux"

#build the static index for ASElux
$ASElux build --gtf annotation.gtf --ref genome.fa --out demo

#align the fasta files
$ASElux align --fa --pe --readLen 50 --index demo --vcf snps.vcf --seqFiles read1.fa read2.fa --out demo_rc.txt
```

The path for  
```{bash}
ASElux="../ASElux"
```
should be re-routed to  
```{bash}
ASElux="../bin/ASElux"
```
in order to access the correct binary file to run ASElux.



After executing the ```script.sh``` file the output should appear as:
```
$/.script.sh
Building static index...
Anootation finished...
annotation created...
loading index: demo
Building dynamic index... 
Start alignment... 
Processed 25944 reads
average speed is: infreads/s
```
In the current ```demo``` directory, new files ```demo.annotation```, ```demo_genome.sa```, ```demo_gene.sa```, and ```demo_rc.txt``` should be generated. ```demo_rc.txt``` should contain the output from the total run and look like:
```
test_82952      51      51      0.5
test_79577      98      98      0.5
test_77887      94      94      0.5
test_74588      100     100     0.5
test_71009      50      50      0.5
test_74333      69      69      0.5
test_59544      98      98      0.5
test_37145      99      99      0.5
test_59371      110     110     0.5
test_56150      105     105     0.5
test_59533      97      97      0.5
test_50858      50      50      0.5
test_27118      55      55      0.5
test_72846      95      95      0.5
test_22557      97      97      0.5
test_22480      91      91      0.5
test_1962       103     103     0.5
test_22345      103     103     0.5
test_73005      112     112     0.5
test_59628      100     100     0.5
test_8469       41      41      0.5
test_20560      99      99      0.5
test_8278       47      47      0.5
```
From package description: "The output file of ASElux has four columns: SNP ID, reference allele count, alternative allele count,and proportion of the reference allele. The SNPs are not sorted in any specific order."

#### Some Notes
The parameters for both building stage of the program and the aligning stage are detailed on the github page linked above. The package is constantly updating to work with a variety of input types and the full description can be found in the link.

## Usage in SMC Research


For our data, we have pooled fastq PRO-Seq data from three separate donors. Our data is in reference to hg38 and thus requires the annotation file to build the static index. The proper files we used were an hg38 fasta file downloaded from: <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/> and the gtf annotation Homo_sapiens.GRCh38.87.gtf.

#### Note
As is, the gtf annotation file is lists CHRM only by #, but we need to have listings by chr# to have proper correspondence with our ```hg38.fa``` file during the static index build. To convert only requires a short perl command: ```perl -ne 'if (/^#/){print $_;}else{print "chr".$_;}' Homo_sapiens.GRCh38.87.gtf > Homo_sapiens.GRCh38.87.chr.gtf``` 

## Building the Static Index

To build the static index using our ```hg38.fa``` file and ```Homo_sapiens.GRCh38.87.chr.gtf``` I created a script (we'll call it ```build.sh```) that enabled a shortcut to the ```ASElux``` executable. It looked something like this:
```
ASElux="../ASElux-master/bin/ASElux"
#build the static index for ASElux
$ASElux build --gtf Homo_sapiens.GRCh38.87.chr.gtf --ref hg38.fa --out index
```
To make executable, run ```$chmod u+x build.sh```. Executing using ```$./build.sh``` should produce:
```
Building static index...
Annotation finished...
annotation created...
```
New files should have been created with names ```index.annotation```, ```index_genome.sa```, ```index_gene.sa```. This process should take the longest as we have a lot of data to index. You now have a static index for alignment and need not build it again.

## Building Dynamic Index and Aligning Reads
Next, we build the individual specific dynamic index and align our data to our static index. For this step I created a script (we'll call it ```align-sample.sh```) which looked something like this:
```
ASElux="../ASElux-master/bin/ASElux"
#align the fastq files
$ASElux align --fq --se --readLen 75 --index index --vcf individual_chr.vcf --seqFiles pooled.fastq --out individual-run_rc.txt
```
Make executable as before with ```chmod u+x align-sample.sh``` and run ```$/.align-sample.sh```.
#### Note
Some problems that could arise:
1. Ensuring proper read length  
To find the read length to use in the alignment stage, I used the FastQC tool to find the average read length of my fastq data
2. Correspondence of that fastq/fasta and the vcf file  
Just as when building the static index, we want to make sure that the first column in the vcf file is individual specific and has the same chr# formatting as our index and fastq file. If not, running the one liner ```perl -ne 'if (/^#/){print $_;}else{print "chr".$_;}' individual.vcf > individual_chr.vcf``` should do the trick.

#### Expected Output
Upon a successful run, we should see output 
```
loading index: index
Building dynamic index...
Start alignment...
Processed ########## reads
average speed is ###### reads/s
```
You should see a new file created named ```individual-run_rc.txt``` with columns for SNP ID, reference allele count, alternative allele count, and proportion of reference allele.
