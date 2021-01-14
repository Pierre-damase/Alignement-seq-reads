"""
Cette méthode permet de déterminer les arguments requis pour le programme.

La validité des arguments est également vérifié.
"""

import argparse
import os
import sys


def valide_file(parser, fichier):
    """
    Test la validité des fichiers passés en paramètre.

      - format du fichier
      - .fasta attendu
      - fichier présent

    Parameter
    ---------
    parser
    fichier: str
        le nom du fichier
    """
    if not fichier.endswith(".fasta") or not os.path.exists(fichier):
        parser.error('Fichier {} non valide - format ou path (voir aide -h)'
                     .format(fichier))
    return fichier


def arguments():
    """
    Détermine les arguments requis pour le programme.

    Return
    ------
    args.reference: str
        la séquence de référence - .fasta
    args.reads:
        les reads - .fasta
    """
    parser = argparse.ArgumentParser()

    # Arguments nécessaires
    group = parser.add_argument_group('fasta', "les .fasta pour l'analyse")
    group.add_argument(
        '-ref', dest="ref", required="True",
        type=lambda x: valide_file(parser, x),
        help="la séquence de référence au format .fasta - située dans le fichier data"
    )
    group.add_argument(
        '-reads', dest="reads", required="True",
        type=lambda x: valide_file(parser, x),
        help="les reads au format .fasta - situés dans le fichier data"
    )

    parser.add_argument('-bwt', dest="bwt", default="space_efficient",
                        choices=['classique', 'space_efficient'],
                        help="""L'algorithme de la transformée de burrows-wheeler à appliquer: soit
                        classique (matrice), soit la version space efficient (par défaut)""")

    return parser.parse_args()


if __name__ == "__main__":
    sys.exit()  # aucune action souhaitée
