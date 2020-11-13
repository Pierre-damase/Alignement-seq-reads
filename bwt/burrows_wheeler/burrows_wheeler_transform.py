"""
Algorithme de burrows-wheeler.
"""

from collections import Counter
import copy
import sys


def burrows_wheeler_transform(seq):
    """
    Cette méthode permet d'appliquer l'algorithme de burrows-wheeter, i.e:

        - générer l'ensemble des permutations de la séquence
        - trier par ordre alphabétique l'ensemble des permutations obtenues
        - extraire la dernière colonne

    Parameter
    ---------
    seq: str
        la séquence de référence

    Return
    ------
    bwt: str
        la séquence extraite, i.e la transformé de burrows-wheeler
    """
    permutations = ['$' + seq]  # la matrice de burrows-wheeler
    for i in range(len(seq)):
        # Réupérer le dernier nucléotide de la séquence d'indice i
        last_element = permutations[i][-1]

        # Concaténer ce dernier élément à la séquence d'indice i (sans prendre
        # en compte le dernier élément) pour générer la séquence d'indice i+1
        tmp = last_element + permutations[i][0:-1]
        permutations.append(tmp)

    permutations.sort()  # trie par ordre alphabétique

    # extraire la dernière colonne
    bwt = ""
    for permut in permutations:
        bwt += permut[-1]  # => bwt = bwt + permut[-1]

    return bwt


def inverse_bwt(bwt):
    """
    Cette méthode permet de revenir à la seq de ref à partir de la transformée.

        - Ajouter dans inverse la transformé en tant que 1ère colonne
        - Trier par ordre alphabétique les lignes de inverse

    Parameter
    ---------
    bwt: str
        la séquence transformée

    Return
    ------
    seq: str
        la séquence d'origine
    """
    inverse = [""] * len(bwt)

    for _ in range(len(bwt)):
        # Ajoute la transformé en tant que 1ère colonne
        tmp = [bwt[i] + inverse[i] for i in range(len(bwt))]

        # Trie par ordre alphabétique les "ligne" la matrice
        inverse = sorted(tmp)

    # Extrait la séquence initiale
    seq = [row for row in inverse if row.startswith('$')][0]
    return seq.strip('$')


def tally_table(bwt):
    """
    Génère une tally table à partir de la transformée

    Parameter
    ---------
    bwt: str
        la séquence transformée

    Return
    ------
    tally: list
        la tally table stocké sous la forme d'une liste de dictionnaire
    """
    dico, tally = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}, []
    for nucleotide in bwt:
        dico[nucleotide] += 1  # Une occurence d'un nucleotide donné
        tally.append(copy.deepcopy(dico))  # Ajoute dans la tally table
    return tally


def nb_inferieurs(bwt):
    """
    Génère le nombre d'occurences de caractères plus petits dans la séquence.

    Et ceux pour chaque caractères.

    Parameter
    ---------
    bwt: str
        la séquence transformée

    Return
    ------
    inferieurs: dict
        un dictionnaire des occurences
    """
    inferieurs = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}
    count = Counter(sorted(bwt))
    nb = 0
    for nucleotide in inferieurs.keys():
        inferieurs[nucleotide] = nb
        nb += count[nucleotide]
    return inferieurs


def generer_ranks(bwt):
    """
    Détermine le rang de chaque caractères de la transformée

    Parameter
    ---------
    bwt: str
        la séquence transformée

    Return
    ------
    ranks: list
        liste des rangs
    """
    dico, ranks = {'$': 0, 'A': 0, 'C': 0, 'G': 0, 'T': 0}, []
    for nucleotide in bwt:
        ranks.append(dico[nucleotide])
        dico[nucleotide] += 1
    return ranks


def generer_position(bwt, inferieurs, ranks):
    """
    Détermine la position des nucléotides de la transformée ds la seq d'origine.

    Une position de -1 est affectée au $.

    Parameter
    ---------
    bwt: str
        la séquence transformée
    inferieurs: dict
        un dictionnaire des occurences
    ranks: list
        liste des rangs

    Return
    ------
    pos: list
        liste des positions
    """
    pos = [-1] * len(bwt)
    seq = ""
    x = len(bwt)-2
    p = 0
    while bwt[p] != "$":
        pos[p] = x
        p = inferieurs[bwt[p]] + ranks[p]
        x -= 1
    return pos


def map(read, inferieurs, tally):
    """
    Permet d'aligner le read avec le génome de référence.

    Parameter
    ---------
    read: str
        le read à aligner
    inferieurs: dict
        un dictionnaire des occurences
    tally: list
        la tally table stocké sous la forme d'une liste de dictionnaire

    Return
    ------
    borne_inf: int
        borne inféreur du vecteur positions
    borne_supp: int
        borne supérieur du vecteur positions
    """
    i = len(read) - 1
    borne_inf = inferieurs[read[i]]
    borne_supp = inferieurs[read[i]] + tally[-1][read[i]] - 1
    while i >= 1 and borne_inf <= borne_supp:
        i -= 1
        borne_inf = inferieurs[read[i]] + tally[borne_inf-1][read[i]]
        borne_supp = inferieurs[read[i]] + tally[borne_supp][read[i]] - 1
    return borne_inf, borne_supp


def determiner_alignement(ref, read, positions, bornes):
    """
    Renvoie l'alignement d'un read donné avec la séquence de référence.

    méthode a eviter sur les sequences de taille importante.

    Parameter
    ---------
    ref: str
        la séquence de référence
    read: str
        le read à aligner
    positions: list
        liste des positions
    bornes: tuple
        les bornes du vecteur positions

    Return
    ------
    alignement: str
        l'alignement du reads avec la séquence de référence
    matches: int
        le nombre de matches
    """
    # Déterminer les positions de départ de chaque motifs dans la séquence
    start_motif = []
    for i in range(bornes[0], bornes[1]+1):
        start_motif.append(positions[i]+1)
    start_motif.sort()
    matches = len(start_motif)

    # Déterminer l'alignement
    alignement = "." * len(ref)
    for pos in start_motif:
        alignement = alignement[:pos] + read + alignement[pos+len(read):]

    return alignement, matches


def determiner_positions(positions, bornes):
    """
    Détermine pour chaque read le nombre de matches.

    Ainsi que les positions de départ de chaque match.

    Parameter
    ---------
    positions: list
        liste des positions
    bornes: tuple
        les bornes du vecteur positions

    Return
    ------
    start_motif: list
        les positions de départ de chaque read dans la séquence
    matches: int
        le nombre d'alignement trouvé
    """
    # Déterminer les positions de départ de chaque motifs dans la séquence
    start_motif = []
    for i in range(bornes[0], bornes[1]+1):
        start_motif.append(positions[i]+1)
    start_motif.sort()
    matches = len(start_motif)

    return start_motif, matches


def find_all(nucleotide, seq):
    """
    Générateur qui renvoie l'index de chaque nucleotide donné dans la séquence.

    Parameter
    ---------
    seq: str
        la séquence de référence
    """
    index = seq.find(nucleotide)
    while index != -1:
        yield index
        index = seq.find(nucleotide, index+1)


def compare(pos1, pos2, seq):
    """

    """
    cpt = 1
    while cpt < len(seq):
        #print("{} & {}".format(seq[pos1:pos1+cpt], seq[pos2:pos2+cpt]))
        if seq[pos1:pos1+cpt] == seq[pos2:pos2+cpt]:
            cpt += 1
        elif seq[pos1:pos1+cpt] > seq[pos2:pos2+cpt]:
            return False
        else:
            return True


def bwt_space_efficient(seq):
    """
    Cette méthode permet d'appliquer l'algorithme de burrows-wheeter tansform.

    Parameter
    ---------
    seq: str
        la séquence de référence

    Return
    ------
    bwt: str
        la séquence extraite, i.e la transformé de burrows-wheeler
    """
    seq = "$" + copy.deepcopy(seq)
    dico = {'A': [], 'C': [], 'G': [], 'T': []}

    for nucleotide in dico:
        print("\rTraitement {}".format(nucleotide), end="")
        for pos in find_all(nucleotide, seq):
            if not dico[nucleotide]:
                dico[nucleotide].append(pos)
            else:
                # Comparer chaque nucleotide donné entre eux
                cpt = 0
                while cpt < len(dico[nucleotide]) and compare(dico[nucleotide][cpt], pos, seq):
                    cpt += 1
                dico[nucleotide].insert(cpt, pos)

    # Générer la transformée
    bwt = seq[-1]
    for value in dico.values():
        for pos in value:
            bwt += seq[pos-1]

    print("\rBurrows-wheeler transform over")

    return bwt


if __name__ == "__main__":
    sys.exit()  # aucune action souhaitée
