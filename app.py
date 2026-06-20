"""
Flask application
"""

import os
import json
from flask import Flask, render_template, request, jsonify, send_file
import io

from bio.toy_dataset import REFERENCE_SEQUENCE, generate_toy_fastq, reconstruction_identity
from bio.fastq_parser import parse_reads, reads_to_fasta
from bio_service import BioinformaticsPipeline, PipelineConfig

app = Flask(__name__)
app.secret_key = 'asg-2026-secret-key'

# Configuration
UPLOAD_FOLDER = os.path.join(os.path.dirname(__file__), 'uploads')
ALLOWED_EXTENSIONS = {'fastq', 'fq', 'fasta', 'fa', 'fna', 'txt'}
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# ============ MAIN PIPELINE PAGE ============

@app.route('/')
def index():
    """Main pipeline interface."""
    return render_template('pipeline.html')

@app.route('/api/pipeline', methods=['POST'])
def run_pipeline():
    """
    Main API endpoint: runs the complete pipeline.
    Expects: FASTQ file + parameters (k, bloom_size, bloom_hashes, solid_threshold)
    Returns: JSON with QC results, assembly results, and contigs for alignment.
    """
    try:
        # Get parameters
        k = request.form.get('k', 31, type=int)
        bloom_size = request.form.get('bloom_size', 10000, type=int)
        bloom_hashes = request.form.get('bloom_hashes', 7, type=int)
        solid_threshold = request.form.get('solid_threshold', 2, type=int)
        
        # Validate parameters
        if k < 1 or k > 127:
            return jsonify({'error': 'k doit être entre 1 et 127'}), 400
        if solid_threshold < 1:
            return jsonify({'error': 'Le seuil de solidité doit être >= 1'}), 400
        
        # Check file
        if 'file' not in request.files:
            return jsonify({'error': 'Aucun fichier sélectionné'}), 400
        
        file = request.files['file']
        if file.filename == '' or not allowed_file(file.filename):
            return jsonify({'error': 'Format de fichier invalide'}), 400
        
        # Read file
        content = file.read().decode('utf-8', errors='ignore')
        
        # Initialize pipeline
        config = PipelineConfig(
            k=k,
            bloom_size_bits=bloom_size,
            bloom_num_hashes=bloom_hashes,
            solid_threshold=solid_threshold
        )
        pipeline = BioinformaticsPipeline(config)
        
        # Load read file (FASTQ or FASTA)
        pipeline.load_reads(file.filename, content)
        
        # Run QC
        qc_result = pipeline.run_qc()
        
        # Run assembly
        assembly_result = pipeline.run_assembly()
        
        # Get sequences for alignment
        sequences_for_alignment = pipeline.get_sequences_for_alignment()
        
        return jsonify({
            'success': True,
            'qc': {
                'num_reads': qc_result.num_reads,
                'total_bases': qc_result.total_bases,
                'min_length': qc_result.min_length,
                'max_length': qc_result.max_length,
                'avg_length': f"{qc_result.avg_length:.1f}",
                'total_kmers': qc_result.total_kmers,
                'unique_kmers': qc_result.unique_kmers,
                'singleton_kmers': qc_result.singleton_kmers,
                'rare_kmers': qc_result.rare_kmers,
                'noise_ratio': qc_result.noise_ratio,
                'histogram': qc_result.histogram
            },
            'assembly': {
                'num_contigs': assembly_result.num_contigs,
                'num_solid_kmers': assembly_result.num_solid_kmers,
                'bloom_fp_estimate': f"{assembly_result.bloom_false_positive_estimate:.2f}%",
                'contigs': assembly_result.contigs
            },
            'alignment': {
                'sequences': sequences_for_alignment
            }
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/align', methods=['POST'])
def align():
    """
    Alignment endpoint: aligns two sequences.
    Expects: seq1, seq2 (JSON body)
    Returns: alignment result
    """
    try:
        data = request.get_json()
        seq1 = data.get('seq1', '').strip()
        seq2 = data.get('seq2', '').strip()
        
        if not seq1 or not seq2:
            return jsonify({'error': 'Deux séquences requises'}), 400
        
        if len(seq1) > 5000 or len(seq2) > 5000:
            return jsonify({'error': 'Séquences trop longues (max 5000 bases)'}), 400
        
        # Validate
        valid_bases = set('ACGT')
        if not all(c in valid_bases for c in seq1.upper()):
            return jsonify({'error': 'Seq1 contient des bases invalides'}), 400
        if not all(c in valid_bases for c in seq2.upper()):
            return jsonify({'error': 'Seq2 contient des bases invalides'}), 400
        
        # Align
        pipeline = BioinformaticsPipeline()
        result = pipeline.align_sequences(seq1, seq2)
        
        return jsonify({
            'success': True,
            'score': result.score,
            'aligned_seq1': result.aligned_seq1,
            'aligned_seq2': result.aligned_seq2,
            'match_line': result.match_line,
            'overlap_start_seq1': result.overlap_start_seq1,
            'overlap_start_seq2': result.overlap_start_seq2,
            'percent_identity': f"{(result.score / max(len(seq1), len(seq2)) * 100):.1f}%"
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/toy_dataset')
def get_toy_dataset():
    """Download toy dataset for demonstration."""
    try:
        fastq_content = generate_toy_fastq(REFERENCE_SEQUENCE)
        return send_file(
            io.BytesIO(fastq_content.encode('utf-8')),
            mimetype='text/plain',
            as_attachment=True,
            download_name='toy_dataset.fastq'
        )
    except Exception as e:
        return str(e), 500

@app.route('/api/convert_fasta', methods=['POST'])
def convert_fasta():
    """Convert uploaded FASTQ/FASTA reads into a FASTA download."""
    try:
        if 'file' not in request.files:
            return jsonify({'error': 'Aucun fichier sélectionné'}), 400
        file = request.files['file']
        if file.filename == '' or not allowed_file(file.filename):
            return jsonify({'error': 'Format de fichier invalide'}), 400
        content = file.read().decode('utf-8', errors='ignore')
        reads = parse_reads(content, file.filename)
        fasta_content = reads_to_fasta(reads)
        return send_file(
            io.BytesIO(fasta_content.encode('utf-8')),
            mimetype='text/plain',
            as_attachment=True,
            download_name=f'{os.path.splitext(file.filename)[0]}.fasta'
        )
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/demo')
def demo():
    """
    Run demo on toy dataset and return results.
    """
    try:
        # Generate toy dataset
        fastq_content = generate_toy_fastq(REFERENCE_SEQUENCE)
        
        # Initialize pipeline
        config = PipelineConfig(k=31, solid_threshold=2)
        pipeline = BioinformaticsPipeline(config)
        pipeline.load_reads('toy_dataset.fastq', fastq_content)
        
        # Run analysis
        qc_result = pipeline.run_qc()
        assembly_result = pipeline.run_assembly()
        sequences_for_alignment = pipeline.get_sequences_for_alignment()
        
        # Compute identity of largest contig
        identity = 0.0
        if pipeline.contigs:
            largest = pipeline.contigs[0]
            identity = reconstruction_identity(largest.sequence)
        
        return jsonify({
            'success': True,
            'reference_length': len(REFERENCE_SEQUENCE),
            'qc': {
                'num_reads': qc_result.num_reads,
                'total_bases': qc_result.total_bases,
                'min_length': qc_result.min_length,
                'max_length': qc_result.max_length,
                'avg_length': f"{qc_result.avg_length:.1f}",
                'total_kmers': qc_result.total_kmers,
                'unique_kmers': qc_result.unique_kmers,
                'singleton_kmers': qc_result.singleton_kmers,
                'rare_kmers': qc_result.rare_kmers,
                'noise_ratio': qc_result.noise_ratio,
                'histogram': qc_result.histogram
            },
            'assembly': {
                'num_contigs': assembly_result.num_contigs,
                'num_solid_kmers': assembly_result.num_solid_kmers,
                'bloom_fp_estimate': f"{assembly_result.bloom_false_positive_estimate:.2f}%",
                'contigs': assembly_result.contigs
            },
            'validation': {
                'largest_contig_length': len(pipeline.contigs[0].sequence) if pipeline.contigs else 0,
                'identity': f"{identity:.1f}%",
                'identity_pass': identity > 98
            },
            'alignment': {
                'sequences': sequences_for_alignment
            }
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(debug=True)

