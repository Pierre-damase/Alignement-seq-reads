"""
Ce module permet de stocker les alignements obtenus dans un fichier texte.

Sous la  fome:

    Read1: Nombre de matches | positions de départ du.es match.es
    Read2: ...
"""

import csv
import sys


def write_transforme(transformee):
    """
    Sauvegarde la transformée dans un fichier .txt.

    Parameter
    ---------
    transformee: str
        la séquence transformée
    """
    with open("transformee.txt", "w") as filout:
        pass


def write_fasta(alignements):
    """
    Écrit dans un fichier .csv.

      - L'id du read
      - Le nombre de matches trouvé
      - La position de chaque match

    Parameter
    ---------
    alignements: dico
        dictionnaire de l'alignement des reads avec la séquence
          - clé: l'id de chaque read
          - value: nombre de matches (int) & positions (list)
    """
    with open("alignements.csv", "w") as filout:
        fields = ["Read", "Matched", "Positions"]
        f_writer = csv.DictWriter(filout, fieldnames=fields)
        f_writer.writeheader()
        for read, value in alignements.items():
            data = {
                "Read": read,
                "Matched": value['Matched'],
                "Positions": value['Positions']
            }
            f_writer.writerow(data)


if __name__ == "__main__":
    sys.exit()  # Aucun action souhaitée
