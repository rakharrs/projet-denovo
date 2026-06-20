"""
Lot 3 - Algorithme de traversee "on-the-fly" (approche type Minia 2).

Idee centrale : on ne construit JAMAIS le Graphe de de Bruijn explicitement
(pas de liste de noeuds/aretes en memoire). A la place, a chaque etape, on
TESTE si un k-mer existe en interrogeant le Filtre de Bloom. Le graphe reste
"conceptuel" : ses aretes ne sont decouvertes qu'au moment ou on en a besoin.

Algorithme, pour un k-mer de depart (seed) :
    1. On regarde les 4 extensions possibles en 3' (on ajoute A, C, G ou T
       a la fin du dernier k-mer du contig courant, puis on retire la
       premiere base pour garder une longueur k -> c'est le k-mer suivant).
    2. On teste chaque extension dans le Bloom Filter.
       - 0 extension valide  -> fin de contig (cul-de-sac).
       - 1 extension valide  -> on avance, le contig s'allonge d'une base.
       - 2+ extensions valides -> embranchement (bifurcation) : on arrete
         ce contig ici et on note les branches possibles (dans cette
         version pedagogique, on n'explore pas recursivement chaque branche
         pour rester lisible, mais l'info est rapportee).

Important : comme le Bloom Filter peut renvoyer des faux positifs, une
extension testee "presente" peut en realite ne jamais avoir existe dans les
donnees d'origine -> ca cree un "chemin fantome" qui degrade la qualite de
l'assemblage (cf. bio/bloom.py et le rapport d'analyse critique).
"""
from dataclasses import dataclass, field

from .bloom import BloomFilter

BASES = "ACGT"


@dataclass
class Contig:
    sequence: str
    seed_kmer: str
    ended_reason: str  # "dead_end" | "branch" | "cycle_detected" | "max_length"
    branch_options: list[str] = field(default_factory=list)  # k-mers candidats si embranchement


def extend_forward(seed_kmer: str, bloom: BloomFilter, used_kmers: set[str],
                    max_length: int = 100_000) -> Contig:
    """Etend un contig en 3' a partir d'un k-mer de depart, jusqu'a cul-de-sac
    ou embranchement. Marque les k-mers visites dans `used_kmers` (cette
    structure ne sert qu'au bookkeeping de la traversee, ce n'est PAS le
    graphe de de Bruijn explicite)."""
    k = len(seed_kmer)
    contig = seed_kmer
    used_kmers.add(seed_kmer)

    while len(contig) - k < max_length:
        last_kmer = contig[-k:]
        candidates = []
        for base in BASES:
            candidate_kmer = last_kmer[1:] + base
            if candidate_kmer in bloom:
                candidates.append((base, candidate_kmer))

        if len(candidates) == 0:
            return Contig(contig, seed_kmer, "dead_end")

        if len(candidates) == 1:
            base, candidate_kmer = candidates[0]
            if candidate_kmer in used_kmers:
                return Contig(contig, seed_kmer, "cycle_detected")
            contig += base
            used_kmers.add(candidate_kmer)
            continue

        # len(candidates) >= 2 : embranchement, on s'arrete et on rapporte les options
        return Contig(contig, seed_kmer, "branch", branch_options=[c for _, c in candidates])

    return Contig(contig, seed_kmer, "max_length")


def has_predecessor(kmer: str, bloom: BloomFilter) -> bool:
    """
    Teste si `kmer` a un 'predecesseur' valide dans le graphe conceptuel,
    c-a-d s'il existe une base B telle que (B + kmer[:-1]) soit elle-meme un
    k-mer solide. Si oui, ce k-mer n'est pas un vrai point de DEPART de
    contig : il sera atteint naturellement en etendant depuis son
    predecesseur. Cette heuristique (utilisee par les assembleurs reels type
    Minia/SPAdes) evite de fragmenter artificiellement un contig simplement
    parce que l'ordre d'iteration des k-mers solides est arbitraire.
    """
    prefix = kmer[:-1]
    for base in BASES:
        if (base + prefix) in bloom:
            return True
    return False


def assemble(solid_kmers: set[str], k: int, bloom: BloomFilter,
             max_length: int = 100_000) -> list[Contig]:
    """
    Lance la traversee en partant prioritairement des k-mers qui n'ont pas
    de predecesseur (= vrais debuts de chemin dans le graphe conceptuel),
    puis traite les k-mers solides restants (ilots non atteints, boucles)
    pour ne perdre aucune donnee. Renvoie la liste des contigs obtenus.
    """
    used: set[str] = set()
    contigs: list[Contig] = []

    starting_points = [km for km in solid_kmers if not has_predecessor(km, bloom)]
    remaining = [km for km in solid_kmers if km not in set(starting_points)]

    for kmer in starting_points + remaining:
        if kmer in used:
            continue
        contig = extend_forward(kmer, bloom, used, max_length=max_length)
        contigs.append(contig)
    contigs.sort(key=lambda c: len(c.sequence), reverse=True)
    return contigs
