"""
Service layer for bioinformatics pipeline.
Orchestrates all the bio modules and provides high-level functionality.
"""

from dataclasses import dataclass, field
from typing import Optional
from bio.fastq_parser import parse_reads
from bio.kmer import count_kmers, frequency_histogram
from bio.alignment import lcs_align
from bio.assembler import assemble
from bio.bloom import BloomFilter


@dataclass
class PipelineConfig:
    """Configuration parameters for the pipeline."""
    k: int = 31
    bloom_size_bits: int = 10_000
    bloom_num_hashes: int = 7
    solid_threshold: int = 2


@dataclass
class QCResult:
    """Quality control results."""
    num_reads: int
    total_bases: int
    min_length: int
    max_length: int
    avg_length: float
    total_kmers: int
    unique_kmers: int
    singleton_kmers: int
    rare_kmers: int
    noise_ratio: float
    histogram: dict  # frequency -> count


@dataclass
class AssemblyResult:
    """Assembly results."""
    num_contigs: int
    num_solid_kmers: int
    contigs: list = field(default_factory=list)
    bloom_false_positive_estimate: float = 0.0


@dataclass
class AlignmentResult:
    """Alignment results."""
    score: int
    aligned_seq1: str
    aligned_seq2: str
    match_line: str
    overlap_start_seq1: int
    overlap_start_seq2: int


class BioinformaticsPipeline:
    """Main pipeline orchestrator."""
    
    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config or PipelineConfig()
        self.reads = []
        self.sequences = []
        self.kmer_counts = None
        self.solid_kmers = set()
        self.bloom = None
        self.contigs = []
        
    def load_reads(self, filename: str, content: str) -> bool:
        """Load and parse a read file in FASTQ/FASTA format."""
        try:
            self.reads = parse_reads(content, filename)
            self.sequences = [r.sequence for r in self.reads]
            return len(self.reads) > 0
        except Exception as e:
            raise ValueError(f"Failed to parse input file: {str(e)}")
    
    def run_qc(self) -> QCResult:
        """Run quality control analysis."""
        if not self.sequences:
            raise ValueError("No sequences loaded")
        
        # K-mer analysis
        self.kmer_counts = count_kmers(self.sequences, self.config.k)
        histogram = frequency_histogram(self.kmer_counts)
        
        # Statistics
        lengths = [len(seq) for seq in self.sequences]
        total_bases = sum(lengths)
        
        singleton_count = sum(1 for count in self.kmer_counts.values() if count == 1)
        rare_count = sum(1 for count in self.kmer_counts.values() if count <= 2)
        noise_ratio = 100.0 * singleton_count / len(self.kmer_counts) if self.kmer_counts else 0.0

        return QCResult(
            num_reads=len(self.sequences),
            total_bases=total_bases,
            min_length=min(lengths),
            max_length=max(lengths),
            avg_length=total_bases / len(self.sequences),
            total_kmers=sum(self.kmer_counts.values()),
            unique_kmers=len(self.kmer_counts),
            singleton_kmers=singleton_count,
            rare_kmers=rare_count,
            noise_ratio=noise_ratio,
            histogram=dict(sorted(histogram.items()))
        )
    
    def run_assembly(self) -> AssemblyResult:
        """Run assembly pipeline."""
        if not self.kmer_counts:
            raise ValueError("Run QC first")
        
        # Select solid k-mers
        self.solid_kmers = {
            km for km, count in self.kmer_counts.items() 
            if count >= self.config.solid_threshold
        }
        
        if not self.solid_kmers:
            raise ValueError(f"No solid k-mers found with threshold={self.config.solid_threshold}")
        
        # Build Bloom Filter
        self.bloom = BloomFilter(
            size_bits=self.config.bloom_size_bits,
            num_hashes=self.config.bloom_num_hashes
        )
        for kmer in self.solid_kmers:
            self.bloom.add(kmer)
        
        # Assemble
        self.contigs = assemble(self.solid_kmers, self.config.k, self.bloom)
        
        # Estimate false positive rate (theoretical)
        n = len(self.solid_kmers)
        m = self.config.bloom_size_bits
        k = self.config.bloom_num_hashes
        fp_rate = (1 - (1 - 1/m) ** (k*n)) ** k  # theoretical formula
        
        return AssemblyResult(
            num_contigs=len(self.contigs),
            num_solid_kmers=len(self.solid_kmers),
            contigs=[
                {
                    'sequence': c.sequence,
                    'length': len(c.sequence),
                    'seed_kmer': c.seed_kmer,
                    'ended_reason': c.ended_reason,
                    'branch_options': c.branch_options
                }
                for c in self.contigs
            ],
            bloom_false_positive_estimate=fp_rate * 100
        )
    
    def align_sequences(self, seq1: str, seq2: str) -> AlignmentResult:
        """Align two sequences."""
        try:
            result = lcs_align(seq1.upper(), seq2.upper(), "Seq1", "Seq2")
            return AlignmentResult(
                score=result.score,
                aligned_seq1=result.aligned_seq1,
                aligned_seq2=result.aligned_seq2,
                match_line=result.match_line,
                overlap_start_seq1=result.overlap_start_seq1,
                overlap_start_seq2=result.overlap_start_seq2
            )
        except Exception as e:
            raise ValueError(f"Alignment failed: {str(e)}")
    
    def get_sequences_for_alignment(self) -> dict:
        """Get available sequences (reads + contigs) for alignment."""
        sequences = {}
        
        # Add reads
        for i, read in enumerate(self.reads[:20]):  # Limit to first 20
            sequences[f'read_{i}'] = {
                'label': f'Read {i} (len={len(read.sequence)})',
                'sequence': read.sequence
            }
        
        # Add contigs
        for i, contig in enumerate(self.contigs):
            sequences[f'contig_{i}'] = {
                'label': f'Contig {i} (len={len(contig.sequence)})',
                'sequence': contig.sequence
            }
        
        return sequences
