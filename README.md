# Programme d'alignement de reads à un génome de référence

## Burrows-wheeler transform

Utilisation de la transformée de burrows-wheeler.

Le jupyter notebook **cours_burrows-wheeler.ipynb** présente les bases de l'algorithme burrows-wheeler transform ainsi que les différents algorithmes mis en place pour la générer (par matrice ou space efficient) et pour réaliser l'alignement de reads à un génome de référence.

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
python -m bwt -ref <REF>.fasta -reads <READS>.fasta -bwt <ALGO>
```

1. Aguments obligatoires

> <REF>.fasta

La séquence de référence au format .fasta

> <READS>.fasta

Les reads au format .fasta

2. Agument optionnel
    
> <ALGO>
    
L'algorithme de la transformée de burrows-wheller, soit classique (matrice) ou space efficient (défaut)
    
## Auteurs

IMBERT Pierre

PENARD Esthel
