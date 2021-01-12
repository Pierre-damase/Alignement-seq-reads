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

## Exécuter le programme

Exécuter la commande

```
python -m bwt <REF>.fasta <READS>.fasta
```

- Aguments obligatoires

> <REF>.fasta

La séquence de référence au format .fasta

> <READS>.fasta

Les reads au format .fasta

## Auteurs

IMBERT Pierre

PENARD Esthel
