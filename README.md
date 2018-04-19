# Objectif 
5 colonies de bactéries E.coli ont été séquencées en NGS.    
Après avoir aligné ces 5 fichiers fastq sur un génome de référence, vous devez identifier les mutations spécifiques à chaque colonie ( Les données ont été simulées avec [wgsim](https://github.com/lh3/wgsim) pour permettre un calcul rapide sur une machine standard ).     
Pour cela, vous devez réaliser un pipeline snakemake décrivant l'ensemble des étapes pour passer d'un fichier [fastq](https://fr.wikipedia.org/wiki/FASTQ) à un fichier [vcf](https://en.wikipedia.org/wiki/Variant_Call_Format).    
Les étapes de votre pipeline sont les suivantes :    

- Indexer le génome de référence avec **bwa** index
- Aligner le fastq sur le génome de référence avec **bwa** (*.fastq > *.sam)
- Convertir le SAM en BAM avec **samtools** (*.sam > *.bam)
- Trier par position les reads alignés avec **samtools** ( *.bam > *.sort.bam)
- Indexer le bam avec **samtools** ( *.bam > *.bai )
- Énumérer les bases séquencées à chaque position **samtools pileup** ( *.sort.bam > *.bcf)
- Détecter les variants significatifs avec **bcftools call** (*.bcf > *.vcf) 
- Compresser les vcf avec **bgzip** (*.vcf > vcf.gz)
- Indexer les vcf avec **tabix** (*.vcf.gz > *.vcf.gz.tbi)

# Installation des dépendances
Pour cet exercice, vous avez besoin de : samtools, bwa, bcftools, tabix, snakemake.

## Installation via conda
Dans la situation (très fréquente) ou vous n'êtes pas administrateur, installer vos dépendances par l'intermédiaire de [conda](https://conda.io/miniconda.html).     
Pour installer conda :    

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    chmod +x Miniconda3-latest-Linux-x86_64.sh      
    ./Miniconda3-latest-Linux-x86_64.sh

Ajouter les dépôts bioconda pour avoir accès au catalogue des outils de bioinformatique:

    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

Crée un environnement appelé "master_big":

        conda create -n master_big

Activer/Désactiver l'environnement: 

       source activate master_big
       source deactivate

Installer les applications dans l'environnement:  

      conda install bwa samtools snakemake bcftools tabix

Télécharger les données du TP: 

    git clone git@github.com:dridk/tp_snakemake.git

### Question 1: 
- Combien de bases dans genom/ecoli.fa ? 
- Combien de reads dans sample1.fastq ? 

# Ligne de commande des étapes 
Les commandes de chaque étape sont décrite ci-dessous.    
Dans un premier temps, essayer manuellement d’exécuter chaque commande à partir du fichier *sample1.fastq*.

Indexer le génome *ecoli.fa*: 

    bwa index genom/ecoli.fa 

### Question 2:
- Pourquoi indexer le génome ? 

Alignement de *sample1.fastq*: 

    bwa mem genom/ecoli.fa sample1.fastq > sample1.sam 

Conversion du SAM en BAM: 

    samtools view -bS sample1.sam > sample1.bam

### Question 3:
- Afficher le contenu du fichier SAM et du fichier BAM avec **less** 
- Quel est la différence entre un fichier SAM et BAM ?  
- Afficher la position (colonne 4) et la séquence (colonne 10) des reads alignés. `` samtools view -F4 sample1.bam |cut -f4,10 ``
Qu'observez vous au niveau des positions ? 

Trier le bam par position: 

    samtools sort sample1.bam > sample1.sort.bam 

### Question 4: 
- Refaire la commande `` samtools view -F4 sample1.sort.bam |cut -f4,10 ``? 
- Qu'observez vous ? 

Indexer le bam:

    samtools index sample1.sort.bam 

### Question 5:
Visualiser votre alignement avec [IGV](http://software.broadinstitute.org/software/igv/) ou samtools:    
    `` samtools tview sample1.sort.bam ``      
    `` samtools tview sample1.sort.bam genom/ecoli.fa `` 

- Pour quelle raison observez vous autant de variation ? 
- Évaluer la profondeur ? 
- Évaluer la couverture ? 

Faire le mpileup :   
Pour chaque position nucléotidique du génome, cette commande dénombre les nucléotides observés sur les reads recouvrant cette position. 

    samtools mpileup -g -f genom/ecoli.fa sample1.sort.bam > sample1.bcf

Faire le variant calling :    
Cette commande détecte les vrais variants du bruit de fond grâce à un modèle statistique. 

    bcftools call -c -v sample1.bcf > sample1.vcf 

Compresser et indexer le fichier 

    bgzip sample1.vcf     # Produit un vcf.gz
    tabix sample1.vcf.gz  # Produit un vcf.gz.tbi

### Question 6:
- Afficher le fichier *sample1.vcf*.
- Quel(s) variant(s) avez-vous trouvé pour l'échantillon *sample1* ? 

## Création du pipeline avec snakemake 

Créer un fichier *Snakefile*. Ce fichier contient des "*règles*" ou "*rule*" permettant de définir comment passer d'un fichier à un autre. Ces règles sont dans l'ordre que vous voulez. En demandant à **snakemake** de produire, par exemple le fichier *sample2.bam*, il trouvera et exécutera lui même l'ensemble des règles pour produire ce fichier.    

Exemple avec la règle de conversion d'un sam en bam:

    rule sam2bam:
        input:
            "{sample}.sam"
        output:
            "{sample}.bam"
        shell:
            "samtools view -b {input} > {output}"

Cette règle peut se lire ainsi :    
" Pour produire le fichier *{qqch}.bam*, j'ai besoin du fichier *{qqch}.sam*. Si le fichier *{qqch}.sam* est absent trouver la règle pour le produire. Et ainsi de suite. 

### Question 7:
- Créer la règle d'alignement vu plus haut.  
- Tester votre règle avec la commande suivante (-n: ne rien faire  -p afficher les commandes)    
    ``snakemake -np sample3.sam``     
    ``snakemake -np sample4.sam ``


### Question 8: 

Essayer de créer les autres règles jusqu'au *{sample}.vcf.gz* à l'aide des commandes définies plus haut. Si vous n'y arrivez pas, vous pouvez vous aider, de la [doc officielle](https://snakemake.readthedocs.io/en/stable/) et en dernier recours de la [correction](https://github.com/dridk/tp_snakemake/blob/master/Snakefile.correction).    
Tester alors votre pipeline :     

    snakemake -np sample1.vcf.gz
    snakemake -np sample2.vcf.gz
    snakemake -np sample3.vcf.gz
    snakemake -np sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz 

Exécuter votre pipeline sur 4 coeurs:

    snakemake -p sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz --cores 4

Forcer l'exécution complète de tout le pipeline (-F) : 

    snakemake -pF sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz --cores 4

### Question 9:
- Modifier n'importe quel fichier et ré-exécuter snakemake. Que se passe t-il ? 

Créer une dernière règle pour combiner l'ensemble des fichiers vcf.gz dans un seul fichier allsample.vcf.gz.

    rule mergeAll : 
        input:
            "sample1.vcf.gz",
            "sample2.vcf.gz",
            "sample3.vcf.gz",
            "sample4.vcf.gz",
            "sample5.vcf.gz"
            
        output:
            "allsample.vcf.gz"
        shell:
            "bcftools merge {input}|bgzip> {output}"

Afficher le graphe d'exécution. Vous aurez peut-être besoin de graphviz. 

    conda install graphviz
    snakemake allsample.vcf.gz --dag|dot|display 

![Graphe du pipeline](https://github.com/dridk/tp_snakemake/blob/master/graph.png)

Exécuter l'ensemble du pipeline : 

    snakemake -p --cores 4 allsample.vcf.gz

#### Question 10: 
Quelles sont les mutations retrouvées dans les 5 colonies bactériennes ? 
