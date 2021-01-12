"""
Cette méthode permet de déterminer les arguments requis pour le programme.

La validité des arguments est également vérifié.
"""

import argparse
import os


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
        'ref', type=lambda x: valide_file(parser, x),
        help="la séquence de référence au format .fasta"
    )
    group.add_argument(
        'reads', type=lambda x: valide_file(parser, x),
        help="les reads au format .fasta"
    )

    args = parser.parse_args()

    return args.ref, args.reads


if __name__ == "__main__":
    sys.exit()  # aucune action souhaitée
