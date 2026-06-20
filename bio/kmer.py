"""
Lot 1 - Analyse de granularite (k-mers).

Un k-mer est simplement une sous-sequence de longueur fixe k, obtenue en
faisant glisser une fenetre le long d'un read.

Exemple avec la sequence "ACGTACG" et k=3 :
    ACG, CGT, GTA, TAC, ACG

Pourquoi compter les k-mers ? Si un k-mer apparait tres souvent dans
l'ensemble des reads, c'est probablement une vraie portion du genome.
S'il apparait une seule fois (ou tres rarement), c'est probablement une
erreur de lecture introduite par le sequenceur. L'histogramme de frequence
des k-mers permet de visualiser cette separation et de choisir un seuil
("solide" vs "bruit") qui sera reutilise au Lot 3.
"""
from collections import Counter

VALID_BASES = set("ACGT")


def split_into_kmers(sequence: str, k: int) -> list[str]:
    """Decoupe une sequence en tous ses k-mers chevauchants."""
    if k <= 0:
        raise ValueError("k doit etre strictement positif")
    if len(sequence) < k:
        return []
    return [sequence[i:i + k] for i in range(len(sequence) - k + 1)]


def count_kmers(sequences: list[str], k: int, ignore_ambiguous: bool = True) -> Counter:
    """
    Compte l'occurrence de chaque k-mer distinct sur l'ensemble des sequences.

    ignore_ambiguous : si True, on ignore les k-mers contenant une base autre
    que A/C/G/T (par exemple un 'N' = base indeterminee par le sequenceur).
    """
    counts: Counter = Counter()
    for seq in sequences:
        for kmer in split_into_kmers(seq, k):
            if ignore_ambiguous and not set(kmer) <= VALID_BASES:
                continue
            counts[kmer] += 1
    return counts


def frequency_histogram(kmer_counts: Counter) -> dict[int, int]:
    """
    Construit l'histogramme de frequence : pour chaque "nombre d'occurrences"
    (1 fois, 2 fois, 3 fois, ...), combien de k-mers DISTINCTS apparaissent
    exactement ce nombre de fois.

    C'est la fameuse courbe bimodale : un pic a basse frequence (erreurs)
    et un pic a plus haute frequence (vrai signal genomique).
    """
    histogram: Counter = Counter()
    for occurrences in kmer_counts.values():
        histogram[occurrences] += 1
    return dict(sorted(histogram.items()))


def suggest_solid_threshold(histogram: dict[int, int]) -> int:
    """
    Heuristique simple pour proposer automatiquement un seuil "solide" :
    on cherche le premier creux (minimum local) de l'histogramme apres le
    pic des erreurs, qui separe statistiquement le bruit du vrai signal.
    Si aucun creux net n'est trouve, on retombe sur un seuil par defaut de 2.
    """
    freqs = sorted(histogram.keys())
    if len(freqs) < 3:
        return 2
    counts = [histogram[f] for f in freqs]
    for i in range(1, len(counts) - 1):
        if counts[i] < counts[i - 1] and counts[i] <= counts[i + 1]:
            return freqs[i] + 1  # seuil = juste apres le creux
    return 2


def solid_kmers(kmer_counts: Counter, threshold: int) -> set[str]:
    """Renvoie l'ensemble des k-mers consideres comme 'solides' (frequence >= seuil)."""
    return {kmer for kmer, c in kmer_counts.items() if c >= threshold}
