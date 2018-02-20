## Objectif 
5 colonies de bactéries E.coli ont été séquencées. À partir des 5 fichiers fastq, réaliser un pipeline et identifier les mutations spécifiques à chaque colonie ( Les données ont été simulées avec wgsim ).

## Installation des dépendances
Pour cet exercice, vous avez besoin de : samtools, bwa,bcftools,snakemake. 

    sudo apt-get install samtools bwa bcftools snakemake

Dans la situation (très fréquente) ou vous n'êtes pas administrateur, installer [conda](https://conda.io/miniconda.html):

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x    Miniconda3-latest-Linux-x86_64.sh      
    ./Miniconda3-latest-Linux-x86_64.sh

Ajouter les dépôts bioconda pour bénéficier de l'ensemble des logiciels bioinformatiques.

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

### Crée un environnement appelé "master_big" 
        conda create -n master_big

### Activer l'environnement   
       source activate master_big

### Installer les applications    
      conda install bwa samtools snakemake bcftools 

### Télécharger les données du TP: 

    wget git@github.com:dridk/tp_snakemake.git

### Question : 
- Combien de bases dans genom/ecoli.fa ? 
- Combien de reads dans sample1.fastq ? 

## Création d'un pipeline brut
Nous allons créer un pipeline grâce à snakemake. Les étapes pour un échantillon (sample1.fastq) sont les suivantes : 

- Indexer le génome 
- Aligner un fastq sur le génome de référence  (*.fastq > *.sam)
- Convertir le SAM en BAM (*.sam > *.bam)
- Trier par position les reads alignés ( *.bam > *.sort.bam)
- Indexer le bam ( *.bam > *.bai )
- Caller les variants ( *.bam > *.vcf)

### indexer le génome ecoli.fa 

    bwa index genom/ecoli.fa 

Question : 
- Quelle est l'utilité d'un index ? 

### Alignement de sample1.fastq 

    bwa mem genom/ecoli.fa sample1.fastq > sample1.sam 

### Conversion du sam en bam 

    samtools view -b sample1.sam > sample1.bam

### Trie du bam par position 

    samtools sort sample1.bam > sample1.sort.bam 

### Index du bam 

    samtools index sample1.sort.bam 

### Visualiser votre alignement 

    samtools tview sample1.sort.bam     samtools tview sample1.sort.bam genom/ecoli.fa

### Variant mpillup 

    samtools mpileup -g -f genom/ecoli.fa sample1.sort.bam > sample1.bcf

### Variant calling 

    bcftools call -c -v sample1.bcf > sample1.vcf 

Question: 
- Quel(s) variant(s) avez-vous trouvé pour l'échantillon sample1 ? 

## Création du pipeline avec snakemake 

Crée un fichier Snakefile. Ce fichier contient des "règles" permettant de définir comment passer d'un fichier à un autre. Ces règles sont dans l'ordre que vous voulez. En demandant à snakemake un fichier, il définira lui-même la séquence de commande à exécuter pour produire ce fichier. 

### Règle d'alignement 

    rule alignement:
        input:
            "{sample}.fastq"
        output:
            "{sample}.sam"
        shell:
            "bwa mem genom/ecoli.fa {input} > {output}"

Tester votre règle sans l'exécuter avec la commande : 

    snakemake -np sample3.sam  
    snakemake -np sample4.sam 

### Autres règles 

Essayer de créer les autres règles jusqu'au vcf en partant des commandes défini plus haut. Si vous n'y arrivez pas, vous pouvez vous aider de la correction.

## Tester vos règles avant d'exécuter 

    snakemake -np sample1.vcf
    snakemake -np sample2.vcf
    snakemake -np sample3.vcf
    snakemake -np sample1.vcf sample2.vcf sample3.vcf sample4.vcf

## Executer sur 4 coeurs

    snakemake -p sample1.vcf sample2.vcf sample3.vcf sample4.vcf --cores 4

## Si un fichier est manquant ?

    rm sample2.sort.bam
    snakemake -p sample2.vcf

## Forcer l'exécution complète 

    snakemake -pF sample1.vcf 

## Combiner les vcfs 
Créer une dernière régle pour combiner l'ensemble des fichiers vcf 

    rule mergeAll : 
        input:
            "sample1.vcf",
            "sample2.vcf",
            "sample3.vcf",
            "sample4.vcf",
            "sample5.vcf"
            
        output:
            "allsample.vcf"
        shell:
            "cat {input} > {output}"

Afficher le graphe d'exécution. Vous aurez peut-être besoin de graphviz. 

    sudo apt install graphviz
    snakemake allsample.vcf --dag|dot|display 

Executer l'ensemble du pipeline : 

    snakemake -p allsample.vcf --cores 4 

### Question: 
Quelles sont les mutations retrouvées ? 
