# Programme d'alignement de reads à un génome de référence

## Burrows-wheeler transform

Le jupyter notebook **cours_burrows-wheeler.ipynb** présente les bases de l'algorithme burrows-wheeler transform ainsi que les différents algorithmes mis en place pour la générer (par matrice ou space efficient) et pour réaliser l'alignement de reads à un génome de référence.

Les fichiers **alignements-classique.csv** & **alignements-space_efficient.csv** renseignent de l'alignement des reads **READSsars_cov_2_1e6.fasta** avec la séquence de référence **Hu-1.fasta**. Fichiers obtenus respectivement avec l'algorithme de burrows-wheeler classique (matrice) et space efficient. 

Les fichiers **Hu-1.fasta** & **READSsars_cov_2_1e6.fasta** du fait de leur taille n'ont pu être déposés sur le github. Pour tester le programme veuiller les déposer dans le fichier **/data**.

## Quick start

1. Clone du répertoire github

> Lien HTTPS

```
git clone https://github.com/Pierre-damase/Alignement-seq-reads.git
```

> Lien SSH

```
git clone git@github.com:Pierre-damase/Alignement-seq-reads.git
```

2. Créer l'environnement conda à partir du fichier bwt.yml (recommandé pour le jupyter notebook)

```
conda env create --file bwt.yml
```

3. Activer l'environnement conda

```
conda activate bwt
```

## Exécuter le programme

- Jupyter notebook

```
jupyter lab
```

- Module bwt

```
python -m bwt -ref data/Hu-1.fasta -reads data/READSsars_cov_2_1e6.fasta
```

Réalise l'alignement des reads de READSsars_cov_2_1e6.fasta avec le fichier de référence Hu-1.fasta en utilisant l'algorithme sapce efficient. Ces fichiers se trouve dans le dossier /data.


## Arguments

```
python -m bwt -ref <REF>.fasta -reads <READS>.fasta -bwt <ALGO>
```

1. Arguments obligatoires

> REF.fasta

La séquence de référence au format .fasta

> READS.fasta

Les reads au format .fasta

2. Argument optionnel
    
> ALGO
    
L'algorithme de la transformée de burrows-wheller, soit classique (matrice) ou space efficient (défaut)
    
## Auteurs

IMBERT Pierre

PENARD Esthel
