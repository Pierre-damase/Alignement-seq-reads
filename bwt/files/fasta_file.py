"""
Ce module permet de lire les deux fichiers fasta passés en paramètre.

  - Le fichier <reference>.fasta de la séquence de référence
  - Le fichier <reads>.fasta des reads

Le module permet de:
  - Lire les fichiers fasta
  - Pour la référence
      Renvoyer la séquence au format str
  - Pour les reads
      Renvoyer un dictionnaire clé: reads id & value: la séquence
"""

import sys


def read_reference(fasta):
    """Lit le fichier <reference>.fasta passé en paramètre.

    Parameter
    ---------
    fasta: str
        le nom du fichier <reference>.fasta

    Return
    ------
    reference: str
        la séquence de référence
    """
    with open(fasta, "r") as filin:
        reference = ""
        for line in filin:
            if not line.startswith(">"):
                reference += line.strip()
    return reference


def read_reads(fasta):
    """Lit le fichier <reads>.fasta passé en paramètre.

    Parameter
    ---------
    fasta: str
        le nom du fichier <reads>.fasta

    Return
    ------
    reads: dico
        dictionnaire des reads - clé: reads id & value: la séquence
    """
    with open(fasta, "r") as filin:
        reads = {}
        for line in filin:
            if line.startswith('>'):
                reads_id = line.split('|')[0].split('>')[1]  # extraire l'id
                reads[reads_id] = ""
            else:
                reads[reads_id] += line.strip()
    return reads


if __name__ == "__main__":
    sys.exit()  # Aucune action souhaitée 
