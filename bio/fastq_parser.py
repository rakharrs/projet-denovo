
from dataclasses import dataclass


@dataclass
class Read:
    """Un fragment unique issu du sequencage."""
    identifier: str
    sequence: str
    quality: str | None = None  # absent si le read vient d'un FASTA

    def __len__(self) -> int:
        return len(self.sequence)

    def mean_quality(self) -> float | None:
        """Score de qualite moyen (echelle Phred, base ASCII 33)."""
        if not self.quality:
            return None
        scores = [ord(c) - 33 for c in self.quality]
        return sum(scores) / len(scores)


def parse_fastq(text: str) -> list[Read]:
    """Parse le contenu brut d'un fichier FASTQ en liste de Read."""
    lines = [l for l in text.splitlines() if l != ""]
    reads = []
    i = 0
    while i < len(lines):
        if not lines[i].startswith("@"):
            raise ValueError(
                f"Ligne {i + 1} invalide : un bloc FASTQ doit commencer par '@', "
                f"trouve : {lines[i][:30]!r}"
            )
        identifier = lines[i][1:].strip()
        sequence = lines[i + 1].strip().upper()
        # lines[i+2] doit etre "+"
        if not lines[i + 2].startswith("+"):
            raise ValueError(f"Ligne {i + 3} invalide : attendu '+', trouve {lines[i + 2][:10]!r}")
        quality = lines[i + 3].strip()
        if len(quality) != len(sequence):
            raise ValueError(
                f"Read {identifier!r} : la sequence ({len(sequence)} caracteres) et la "
                f"qualite ({len(quality)} caracteres) n'ont pas la meme longueur"
            )
        reads.append(Read(identifier, sequence, quality))
        i += 4
    return reads


def parse_fasta(text: str) -> list[Read]:
    """Parse le contenu brut d'un fichier FASTA en liste de Read (sans qualite)."""
    reads = []
    identifier = None
    seq_chunks: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if identifier is not None:
                reads.append(Read(identifier, "".join(seq_chunks).upper()))
            identifier = line[1:].strip()
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if identifier is not None:
        reads.append(Read(identifier, "".join(seq_chunks).upper()))
    return reads


def parse_reads(text: str, filename: str) -> list[Read]:
    """Detecte le format (FASTQ ou FASTA) a partir de l'extension et parse."""
    lower = filename.lower()
    if lower.endswith((".fastq", ".fq")):
        return parse_fastq(text)
    if lower.endswith((".fasta", ".fa", ".fna")):
        return parse_fasta(text)
    # fallback : on devine a partir du premier caractere non vide
    stripped = text.lstrip()
    if stripped.startswith("@"):
        return parse_fastq(text)
    if stripped.startswith(">"):
        return parse_fasta(text)
    raise ValueError("Format non reconnu : le fichier doit etre un .fastq/.fq ou .fasta/.fa")


def reads_to_fasta(reads: list[Read]) -> str:
    """Convertit une liste de Read en texte FASTA (Lot 1 : conversion selective)."""
    lines = []
    for r in reads:
        lines.append(f">{r.identifier}")
        # on retourne a la ligne tous les 70 caracteres, convention FASTA classique
        for i in range(0, len(r.sequence), 70):
            lines.append(r.sequence[i:i + 70])
    return "\n".join(lines) + "\n"
