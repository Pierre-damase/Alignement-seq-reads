"""
Ce module permet d'aligner des reads à une séquence de référence.

La transformée de Burrows-wheeler est utilisée pour ce faire.

Usage
-----

    $ python -m bwt <REF> <READS>

    avec REF la séquence de référence au format .fasta
         READS les reads au format .fasta

    ps: être dans le dossier génomique
"""

from bwt.argument.definir_argument import arguments
from bwt.burrows_wheeler import burrows_wheeler_transform as bwt
from bwt.files import fasta_file as fasta


def main():
    # Définition des arguments du programme
    ref, reads = arguments()

    # Lecture des fichiers fasta
    reference = fasta.read_reference(ref)
    reads = fasta.read_reads(reads)

    ## Burrows-wheeler transform ##

    # La transformée de burrows-wheeler
    transformee = bwt.burrows_wheeler_transform(reference)

    # La séquence d'origine déterminée à partir de la transformée
    initial = bwt.inverse_bwt(transformee)

    # La tally table générée à partir de la transformée
    tally = bwt.tally_table(transformee)

    # Génère le nombre de caractères plus petits, et ceux pour chaque caractères
    inferieurs = bwt.nb_inferieurs(transformee)

    # Détermine le rang des caractères de la transformée
    ranks = bwt.generer_ranks(transformee)

    # La position des nucléotides de la transformée dans la séquence d'origine
    positions = bwt.generer_position(transformee, inferieurs, ranks)

    # L'alignement des reads avec la séquence de référence
    for id_read, read in reads.items():
        bornes = bwt.map(read, inferieurs, tally)
        alignement = bwt.afficher_alignement(reference, read, positions, bornes)

        print("Sequence: {}".format(reference))
        print("   Reads: {}".format(alignement), end="\n\n")


if __name__ == "__main__":
    main()
