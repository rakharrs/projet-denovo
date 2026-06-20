
import hashlib
import math


class BloomFilter:
    def __init__(self, size_bits: int, num_hashes: int):
        if size_bits <= 0 or num_hashes <= 0:
            raise ValueError("size_bits et num_hashes doivent etre positifs")
        self.size_bits = size_bits
        self.num_hashes = num_hashes
        self.bit_array = bytearray((size_bits + 7) // 8)  # 1 bit reel par bit logique
        self.inserted_count = 0  # n : nombre d'insertions effectuees (pour les stats)

    def _get_bit(self, index: int) -> int:
        byte_index, bit_offset = divmod(index, 8)
        return (self.bit_array[byte_index] >> bit_offset) & 1

    def _set_bit(self, index: int) -> None:
        byte_index, bit_offset = divmod(index, 8)
        self.bit_array[byte_index] |= (1 << bit_offset)

    def _hash_positions(self, item: str):
        """Genere k positions dans le tableau de bits a partir de k 'sels' differents.
        On simule k fonctions de hachage independantes avec une seule fonction
        cryptographique (sha256) en faisant varier un sel i"""
        item_bytes = item.encode("utf-8")
        for i in range(self.num_hashes):
            digest = hashlib.sha256(item_bytes + i.to_bytes(4, "big")).digest()
            position = int.from_bytes(digest[:8], "big") % self.size_bits
            yield position

    def add(self, item: str) -> None:
        for pos in self._hash_positions(item):
            self._set_bit(pos)
        self.inserted_count += 1

    def __contains__(self, item: str) -> bool:
        return all(self._get_bit(pos) == 1 for pos in self._hash_positions(item))







    def estimated_false_positive_rate(self) -> float:
        """p ≈ (1 - e^(-k*n/m)) ^ k, avec n = inserted_count."""
        k, n, m = self.num_hashes, self.inserted_count, self.size_bits
        if n == 0:
            return 0.0
        return (1 - math.exp(-k * n / m)) ** k

    def fill_ratio(self) -> float:
        """Proportion de bits a 1 dans le tableau (indicateur visuel de saturation)."""
        ones = sum(bin(byte).count("1") for byte in self.bit_array)
        return ones / self.size_bits

    @staticmethod
    def optimal_num_hashes(size_bits: int, expected_items: int) -> int:
        """k_opt = (m/n) * ln(2), arrondi et borne a au moins 1."""
        if expected_items <= 0:
            return 1
        return max(1, round((size_bits / expected_items) * math.log(2)))

    @staticmethod
    def optimal_size_bits(expected_items: int, target_fp_rate: float) -> int:
        """m = -(n * ln(p)) / (ln(2))^2, pour un taux de faux positifs cible p."""
        if expected_items <= 0 or not (0 < target_fp_rate < 1):
            raise ValueError("expected_items > 0 et 0 < target_fp_rate < 1 requis")
        return max(8, math.ceil(-(expected_items * math.log(target_fp_rate)) / (math.log(2) ** 2)))












def memory_comparison(num_distinct_kmers: int, kmer_length: int, bloom_size_bits: int) -> dict:
    """
    Compare la memoire occupee par le Bloom Filter face a un dictionnaire/set
    Python standard qui stockerait explicitement chaque k-mer distinct.

    Pour le dictionnaire : on estime le cout par cle a ~ (overhead Python str
    + structure du set), une approximation courante est 49 octets d'overhead
    par objet str + 1 octet par caractere, plus l'overhead de la table de
    hachage du set lui-meme (~50-70 octets par entree en pratique CPython).
    """
    bytes_per_kmer_in_set = 49 + kmer_length + 60  # approximation CPython realiste
    dict_bytes = num_distinct_kmers * bytes_per_kmer_in_set
    bloom_bytes = bloom_size_bits / 8
    return {
        "dict_set_bytes": dict_bytes,
        "dict_set_mb": dict_bytes / (1024 ** 2),
        "bloom_bytes": bloom_bytes,
        "bloom_mb": bloom_bytes / (1024 ** 2),
        "ratio": dict_bytes / bloom_bytes if bloom_bytes else float("inf"),
    }
