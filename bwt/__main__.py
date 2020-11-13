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

import sys

from bwt.argument.definir_argument import arguments
from bwt.burrows_wheeler import burrows_wheeler_transform as bwt
from bwt.files import read_fasta_files as fasta
from bwt.files import write_file as save


def generer_bwt(reference):
    ## Burrows-wheeler transform ##

    # La transformée de burrows-wheeler
    # transformee = bwt.burrows_wheeler_transform(reference)
    transformee = bwt.bwt_space_efficient(reference)

    # La séquence d'origine déterminée à partir de la transformée
    # initial = bwt.inverse_bwt(transformee)

    # La tally table générée à partir de la transformée
    tally = bwt.tally_table(transformee)

    # Génère le nombre de caractères plus petits, et ceux pour chaque caractères
    inferieurs = bwt.nb_inferieurs(transformee)

    # Détermine le rang des caractères de la transformée
    ranks = bwt.generer_ranks(transformee)

    # La position des nucléotides de la transformée dans la séquence d'origine
    positions = bwt.generer_position(transformee, inferieurs, ranks)

    return tally, inferieurs, positions


def generer_alignements(alignements, reference, reads):
    # L'alignement des reads avec la séquence de référence

    tally, inferieurs, positions = generer_bwt(reference)

    for id_read, read in reads.items():
        print("\rRead {}".format(id_read), end="")
        bornes = bwt.map(read, inferieurs, tally)
        pos, matches = bwt.determiner_positions(positions, bornes)

        if not id_read in alignements.keys():
            alignements[id_read] = {'Matched': 0, 'Positions': []}

        alignements[id_read]['Matched'] += matches
        alignements[id_read]['Positions'].extend(pos)

    print("\nAlignements over\n")

    return alignements


def inverse_seq(seq):
    seq_inverse = ""
    dico = {'A': "T", 'T': "A", 'C': "G", 'G': "C"}

    for nucleotide in seq:
        seq_inverse += dico[nucleotide]

    return seq_inverse


def main():
    # Définition des arguments du programme
    ref, reads = arguments()

    # Lecture des fichiers fasta
    reference = fasta.read_reference(ref)
    reads = fasta.read_reads(reads)

    # Séquence inverse
    reference_inverse = inverse_seq(reference)

    # Alignements
    alignements = {}
    alignements = generer_alignements(alignements, reference, reads)
    alignements = generer_alignements(alignements, reference_inverse, reads)

    save.write_fasta(alignements)


if __name__ == "__main__":
    main()
