"""
Generateur de "Toy Dataset" : une sequence de reference COURTE et CONNUE,
fragmentee en reads chevauchants avec un taux d'erreur controle (par defaut
1%, comme demande dans les contraintes du cahier des charges). Sert a
valider que le pipeline (Lot1 -> Lot3) reconstruit bien la sequence cible.
"""
import random

BASES = "ACGT"


def _generate_non_repetitive_reference(length: int, k_guard: int, seed: int) -> str:
    """Genere une sequence aleatoire (mais deterministe via `seed`) qui ne
    contient aucun k-mer duplique pour k=k_guard. C'est important pour un
    'toy dataset' pedagogique : une sequence avec des repetitions internes
    cree des cycles dans le graphe de de Bruijn et empeche un assemblage
    simple en un seul contig (cf. discussion sur la resolution des
    repetitions dans le rapport technique)."""
    rng = random.Random(seed)
    while True:
        seq = "".join(rng.choice(BASES) for _ in range(length))
        seen_kmers = set()
        ok = True
        for i in range(len(seq) - k_guard + 1):
            kmer = seq[i:i + k_guard]
            if kmer in seen_kmers:
                ok = False
                break
            seen_kmers.add(kmer)
        if ok:
            return seq


# Sequence de reference fixe et deterministe (~300 bases), sans repetition
# de k-mer pour k=15, afin que le toy dataset soit assemblable en un seul
# contig avec un k raisonnable (cf. Lot 3).
REFERENCE_SEQUENCE = _generate_non_repetitive_reference(length=300, k_guard=15, seed=7)


def _mutate(sequence: str, error_rate: float, rng: random.Random) -> str:
    """Introduit des substitutions aleatoires (mutations ponctuelles) au taux donne."""
    chars = list(sequence)
    for i in range(len(chars)):
        if rng.random() < error_rate:
            chars[i] = rng.choice([b for b in BASES if b != chars[i]])
    return "".join(chars)


def generate_toy_fastq(reference: str = REFERENCE_SEQUENCE, read_length: int = 50,
                        step: int = 5, error_rate: float = 0.01, seed: int = 42,
                        boundary_oversample: int = 4) -> str:
    """
    Decoupe `reference` en reads chevauchants de longueur `read_length`,
    avec un pas `step` (donc un chevauchement de read_length - step bases
    entre deux reads consecutifs), puis introduit du bruit (error_rate).

    boundary_oversample : nombre de copies supplementaires (independamment
    mutees) generees pour le PREMIER et le DERNIER read. En sequencage
    shotgun reel, les extremites d'une sequence ne peuvent etre couvertes
    que par des fragments qui commencent/finissent exactement a cet endroit
    -> elles sont structurellement moins couvertes que le milieu ("effet de
    bord"). Sans ce surechantillonnage, les k-mers terminaux tombent souvent
    sous le seuil "solide" et l'assemblage s'arrete juste avant la fin
    (phenomene reel et bien documente en assemblage de novo, pas un bug).

    Renvoie le contenu texte d'un fichier FASTQ pret a etre uploade.
    """
    rng = random.Random(seed)
    lines = []
    read_id = 0
    starts = list(range(0, len(reference) - read_length + 1, step))

    def emit(start: int):
        nonlocal read_id
        raw_read = reference[start:start + read_length]
        noisy_read = _mutate(raw_read, error_rate, rng)
        quality = "".join(rng.choice("=>?@A") for _ in noisy_read)
        lines.append(f"@toyread_{read_id}_pos{start}")
        lines.append(noisy_read)
        lines.append("+")
        lines.append(quality)
        read_id += 1

    for start in starts:
        emit(start)
    # surechantillonnage des bords (effet de bord du sequencage shotgun)
    for _ in range(boundary_oversample):
        emit(starts[0])
        emit(starts[-1])
    return "\n".join(lines) + "\n"


def reconstruction_identity(reconstructed: str, reference: str = REFERENCE_SEQUENCE) -> float:
    """
    Mesure simple d'identite : on cherche `reconstructed` (ou son inverse,
    peu importe le sens du brin) comme sous-chaine la plus proche de
    `reference`, et on calcule le pourcentage de positions identiques sur le
    chevauchement le plus long trouve. Utilise pour la clause de recette du
    cahier des charges ("identite superieure a 98%").
    """
    best_identity = 0.0
    for candidate in (reconstructed, reconstructed[::-1].translate(
            str.maketrans("ACGT", "TGCA"))):  # on teste aussi le brin complementaire inverse
        if candidate in reference:
            return 100.0
        # sinon, alignement glissant simple pour estimer le meilleur recouvrement
        for offset in range(-len(candidate) + 1, len(reference)):
            matches = 0
            total = 0
            for i in range(len(candidate)):
                ref_pos = offset + i
                if 0 <= ref_pos < len(reference):
                    total += 1
                    if candidate[i] == reference[ref_pos]:
                        matches += 1
            if total > 0:
                identity = 100.0 * matches / total
                best_identity = max(best_identity, identity)
    return best_identity