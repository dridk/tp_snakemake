## Objectif 
5 colonies de bactéries E.coli ont été séquencées. À partir des 5 fichiers fastq, réaliser un pipeline et identifier les mutations spécifiques à chaque colonie ( Les données ont été simulées avec [wgsim](https://github.com/lh3/wgsim) ).

## Installation des dépendances
Pour cet exercice, vous avez besoin de : samtools, bwa,bcftools,snakemake. 

    sudo apt-get install samtools bwa bcftools snakemake

Dans la situation (très fréquente) ou vous n'êtes pas administrateur, installer vos dépendances par l'intermédiaire de [conda](https://conda.io/miniconda.html):

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh      
    ./Miniconda3-latest-Linux-x86_64.sh

Ajouter les dépôts bioconda pour avoir accès au catalogue des outils de bioinformatique

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

### Crée un environnement appelé "master_big" 
        conda create -n master_big

### Activer l'environnement   
       source activate master_big

### Installer les applications requises   
      conda install bwa samtools snakemake bcftools 

### Télécharger les données du TP: 
    git clone git@github.com:dridk/tp_snakemake.git

#### Question 1: 
- Combien de bases dans genom/ecoli.fa ? 
- Combien de reads dans sample1.fastq ? 

## Création d'un pipeline brut
Les étapes pour obtenir les variants de l'échantillon (sample1.fastq) sont les suivantes : 

- Indexer le génome avec bwa index
- Aligner un fastq sur le génome de référence avec bwa (*.fastq > *.sam)
- Convertir le SAM en BAM avec samtools (*.sam > *.bam)
- Trier par position les reads alignés avec samtools ( *.bam > *.sort.bam)
- Indexer le bam avec samtools ( *.bam > *.bai )
- Énumérer les bases à chaque position - pileup ( *.sort.bam > *.bcf)
- Détecter les variants significatifs - variant calling (*.bcf > *.vcf) 

### indexer le génome ecoli.fa 

    bwa index genom/ecoli.fa 

#### Question 2:
- Quelle est l'utilité d'un index ? 

### Alignement de sample1.fastq 

    bwa mem genom/ecoli.fa sample1.fastq > sample1.sam 

#### Question 3 : 
Afficher la position (colonne 4) et la sequence (colonne 10) des reads alignés.
Qu'observez vous au niveau des positions ?    

    samtools view -F4 sample1.sam |cut -f4,10

### Conversion du sam en bam 

    samtools view -b sample1.sam > sample1.bam

#### Question 4 : 
Quel est l’intérêt d'un fichier binaire ?   

### Trie du bam par position 

    samtools sort sample1.bam > sample1.sort.bam 

#### Question 5: 
Refaite la commande de la question 3 sur le fichier sample1.sort.bam généré. 

### Index du bam 

    samtools index sample1.sort.bam 

### Visualiser votre alignement (hors pipeline)

    samtools tview sample1.sort.bam
    samtools tview sample1.sort.bam genom/ecoli.fa

#### Question 6:
- Pour quelle raison observez vous autant de variation ? 
- Évaluer la profondeur ? 
- Évaluer la couverture ? 

### Variant mpileup
Pour chaque position nucléotidique du génome, cette commande dénombre les nucléotides observés sur les reads recouvrant cette position. 

    samtools mpileup -g -f genom/ecoli.fa sample1.sort.bam > sample1.bcf

### Variant calling 
Cette commande détecte les vrais variants du bruit de fond grâce à un modèle statistique. 

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

#### Question 6: Tester votre règle
-p : Afficher les commandes
-n : Ne pas les exécuter 

    snakemake -np sample3.sam  
    snakemake -np sample4.sam 

### Question 7 : Créer toutes les règles jusqu'au VCF

Essayer de créer les autres règles jusqu'au vcf en partant des commandes défini plus haut. Si vous n'y arrivez pas, vous pouvez vous aider, de la [doc officielle](https://snakemake.readthedocs.io/en/stable/) et en dernier recours de la [correction](https://github.com/dridk/tp_snakemake/blob/master/Snakefile.correction).    
Vous pouvez alors demander à snakemake comment générer un fichier donnée en argument:

    snakemake -np sample1.vcf
    snakemake -np sample2.vcf
    snakemake -np sample3.vcf
    snakemake -np sample1.vcf sample2.vcf sample3.vcf sample4.vcf

## Executer sur 4 coeurs

    snakemake -p sample1.vcf sample2.vcf sample3.vcf sample4.vcf --cores 4

## Si un fichier est manquant ?

    rm sample2.sort.bam
    snakemake -p sample2.vcf

## Forcer l'exécution de toutes les règles  

    snakemake -pF sample1.vcf 

## Combiner les vcfs 
Créer une dernière règle pour combiner l'ensemble des fichiers vcf 

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
