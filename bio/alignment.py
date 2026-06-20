"""
Lot 2 - Alignement par programmation dynamique.

On cherche la Plus Longue Sous-Sequence Commune (LCS) entre deux reads.
"Sous-sequence" ne veut pas dire "sous-chaine contigue" : les lettres
trouvees peuvent etre espacees, mais doivent rester dans le meme ordre
dans les deux sequences.

Exemple : LCS("ACGTAC", "AGTC") = "AGTC" (longueur 4)

Principe de la programmation dynamique :
On construit une matrice dp ou dp[i][j] = longueur de la LCS entre les
i premiers caracteres de seq1 et les j premiers caracteres de seq2.

    - si seq1[i-1] == seq2[j-1] : dp[i][j] = dp[i-1][j-1] + 1  (on "matche")
    - sinon                     : dp[i][j] = max(dp[i-1][j], dp[i][j-1])

Complexite : O(n*m) en temps ET en memoire (n, m = longueurs des deux
sequences), car on remplit/stocke une matrice (n+1) x (m+1).

Pour "remonter" la matrice et reconstruire un alignement lisible (et pas
seulement la longueur), on fait un backtracking depuis dp[n][m] :
on suit le chemin qui a produit chaque valeur, ce qui nous dit, a chaque
position, si on a un "match" (les deux sequences avancent ensemble) ou un
"gap" (une seule des deux sequences avance, l'autre affiche un '-').
"""
from dataclasses import dataclass


@dataclass
class AlignmentResult:
    score: int               # longueur de la LCS = nombre de bases qui matchent
    seq1_id: str
    seq2_id: str
    aligned_seq1: str         # avec '-' pour les gaps
    aligned_seq2: str         # avec '-' pour les gaps
    match_line: str           # '|' sous chaque position qui matche, ' ' sinon
    overlap_start_seq1: int    # position (0-indexee) du 1er caractere matche dans seq1
    overlap_start_seq2: int    # position (0-indexee) du 1er caractere matche dans seq2
    length_seq1: int
    length_seq2: int

    def pretty_print(self, width: int = 80) -> str:
        """Rendu texte multi-lignes, decoupe en blocs de `width` caracteres."""
        blocks = []
        for i in range(0, len(self.aligned_seq1), width):
            s1 = self.aligned_seq1[i:i + width]
            mid = self.match_line[i:i + width]
            s2 = self.aligned_seq2[i:i + width]
            blocks.append(f"{self.seq1_id:>12}: {s1}\n{'':>12}  {mid}\n{self.seq2_id:>12}: {s2}")
        return "\n\n".join(blocks)


def lcs_align(seq1: str, seq2: str, seq1_id: str = "read1", seq2_id: str = "read2") -> AlignmentResult:
    """Calcule la LCS entre seq1 et seq2 et reconstruit un alignement lisible."""
    n, m = len(seq1), len(seq2)

    # --- 1. Remplissage de la matrice DP : O(n*m) en temps et en memoire ---
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    score = dp[n][m]

    # --- 2. Backtracking : on remonte de dp[n][m] vers dp[0][0] ---
    aligned1_rev: list[str] = []
    aligned2_rev: list[str] = []
    match_rev: list[str] = []
    first_match_i, first_match_j = None, None  # position du dernier match trouve en remontant

    i, j = n, m
    while i > 0 and j > 0:
        if seq1[i - 1] == seq2[j - 1] and dp[i][j] == dp[i - 1][j - 1] + 1:
            aligned1_rev.append(seq1[i - 1])
            aligned2_rev.append(seq2[j - 1])
            match_rev.append("|")
            first_match_i, first_match_j = i - 1, j - 1  # mis a jour a chaque match, donc
            # a la fin de la boucle ce sera le match le plus en amont (= debut du chevauchement)
            i -= 1
            j -= 1
        elif dp[i - 1][j] >= dp[i][j - 1]:
            aligned1_rev.append(seq1[i - 1])
            aligned2_rev.append("-")
            match_rev.append(" ")
            i -= 1
        else:
            aligned1_rev.append("-")
            aligned2_rev.append(seq2[j - 1])
            match_rev.append(" ")
            j -= 1

    # caracteres restants en debut de sequence (hors matrice, avant i=0 ou j=0)
    while i > 0:
        aligned1_rev.append(seq1[i - 1]); aligned2_rev.append("-"); match_rev.append(" ")
        i -= 1
    while j > 0:
        aligned1_rev.append("-"); aligned2_rev.append(seq2[j - 1]); match_rev.append(" ")
        j -= 1

    overlap_start_seq1 = first_match_i if first_match_i is not None else 0
    overlap_start_seq2 = first_match_j if first_match_j is not None else 0

    return AlignmentResult(
        score=score,
        seq1_id=seq1_id,
        seq2_id=seq2_id,
        aligned_seq1="".join(reversed(aligned1_rev)),
        aligned_seq2="".join(reversed(aligned2_rev)),
        match_line="".join(reversed(match_rev)),
        overlap_start_seq1=overlap_start_seq1,
        overlap_start_seq2=overlap_start_seq2,
        length_seq1=n,
        length_seq2=m,
    )
