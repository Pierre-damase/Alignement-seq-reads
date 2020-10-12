<span style="color: #8b1538; font-size: 35px;">**Génomique**</span>

> But

Comment aligner des read sur un génome de référence de manière efficace ?

# Introduction

- Génome humain: 3 400 Mpb - de l'ordre du milliard
- Bactérie type escherichia coli: 4.64 Mpb - de l'ordre du million
- Virus type grippe: 0.013 Mpb

On ne peut donc pas aligner les read sur l'ensemble d'une génome.

<br>

# Solution de base - index sur les kmer

> Principe

- Création du dictionnaire de kmer, i.e l'**index**
- Rechechercher dans l'index

A chaque kmer est associé un ensemble de position dans le génome de réference. Il est donc plus aisé de faire l'alignement avec les read ensuite.

> Problème

- Kmer de taille 1: 4
- Kmer de taille 2: 4^2
- Kmer de taille n: 4^n 

Très rapidement le nombre de kmer possible explose.

<br>

# Transformés de Burrows-wheeler (94)

A l'origine, index pour la compression.

## Exemple

Référence: A T A T C G T

1. Générer l'ensemble des permutations de la séquence de référence

Le \$ a été choisi car c'est le caréctère juste avant le A - ordre alphabétique.

\begin{matrix}
\$ & A & T & A & T & C & G & T \\
 T &\$ & A & T & A & T & C & G \\
 G & T &\$ & A & T & A & T & C \\
 C & G & T &\$ & A & T & A & T \\
 T & C & G & T &\$ & A & T & A \\
 A & T & C & G & T &\$ & A & T \\
 T & A & T & C & G & T &\$ & A \\
 A & T & A & T & C & G & T &\$
\end{matrix}
 
 2. Trier par ordre alphabétique l'ensemble des permutations obtenues
 
\begin{matrix}
\$ & A & T & A & T & C & G & T \\
 A & T & A & T & C & G & T &\$ \\
 A & T & C & G & T &\$ & A & T \\
 C & G & T &\$ & A & T & A & T \\
 G & T &\$ & A & T & A & T & C \\
 T &\$ & A & T & A & T & C & G \\
 T & A & T & C & G & T &\$ & A \\
 T & C & G & T &\$ & A & T & A
\end{matrix}
 
 3. Sélectionner la dernière colonne - appelé la transformé de Burrows-wheeler
 
 <span style="color: #8b1538; font-size: 18px;">T \$ T T C G G A A</span>
  
 ## Propriété
 
 - Cette transformé a tendance a regroupé des patternes/mots équivalents ensemble (appelé **motifs**) - compression
 - A partir de la dernière colonne on peut générer la première (la réciproque est fausse) - toujours stocker la dernière
 - Les éléments d'un même motif sont toujours dans le même ordre quelque soit la colonne que l'on regarde