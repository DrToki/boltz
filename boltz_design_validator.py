#!/usr/bin/env python3
"""
Boltz Protein Design Validator
Simple pipeline for validating designed binders against reference complex structures.

Usage:
    python boltz_design_validator.py reference.pdb --target-chain A --binder-chain B --binder-sequences sequences.txt
"""

import argparse
import logging
from pathlib import Path
from typing import List
import json
import sys
from boltz_design_pipeline.boltz_pipeline import BoltzDesignPipeline
from boltz_design_pipeline.structure_utils import ValidationResult, PDBComplexParser

def load_binder_sequences(sequences_file: Path) -> List[str]:
    """Load binder sequences from file (one per line)"""
    with sequences_file.open() as f:
        sequences = [line.strip() for line in f if line.strip()]
    return sequences

def write_results_summary(results: List[ValidationResult], output_file: Path):
    """Write results to JSON file with metrics summary, separating successful and failed results"""
    # Separate successful and failed results
    successful_results = [r for r in results if r.prediction_success]
    failed_results = [r for r in results if not r.prediction_success]
    
    if not results:
        summary = {
            "validation_results": [],
            "failed_predictions": [],
            "summary_stats": {
                "total_binders": 0,
                "successful_predictions": 0,
                "failed_predictions": 0,
                "avg_complex_plddt": 0.0,
                "avg_complex_rmsd": 0.0,
                "avg_complex_rmsd_aligned": 0.0,
                "high_confidence_count": 0,
                "good_rmsd_count": 0
            }
        }
    else:
        # Calculate summary stats from successful results only
        plddt_scores = [r.confidence_scores.get("complex_plddt", 0) for r in successful_results if r.confidence_scores]
        rmsd_scores = [r.rmsd_metrics.get("complex_rmsd", 0) for r in successful_results if r.rmsd_metrics.get("complex_rmsd", -1) >= 0]
        rmsd_aligned_scores = [r.rmsd_metrics.get("complex_rmsd_aligned", 0) for r in successful_results if r.rmsd_metrics.get("complex_rmsd_aligned", -1) >= 0]
        
        summary = {
            "validation_results": [],
            "failed_predictions": [],
            "summary_stats": {
                "total_binders": len(results),
                "successful_predictions": len(successful_results),
                "failed_predictions": len(failed_results),
                "avg_complex_plddt": sum(plddt_scores) / len(plddt_scores) if plddt_scores else 0.0,
                "avg_complex_rmsd": sum(rmsd_scores) / len(rmsd_scores) if rmsd_scores else 0.0,
                "avg_complex_rmsd_aligned": sum(rmsd_aligned_scores) / len(rmsd_aligned_scores) if rmsd_aligned_scores else 0.0,
                "high_confidence_count": sum(1 for score in plddt_scores if score > 0.8),
                "good_rmsd_count": sum(1 for score in rmsd_scores if score < 2.0)
            }
        }
    
    # Add successful results
    for i, result in enumerate(successful_results):
        summary["validation_results"].append({
            "binder_id": f"binder_{i:03d}",
            "sequence": result.binder_sequence,
            "confidence_scores": result.confidence_scores,
            "rmsd_metrics": result.rmsd_metrics,
            "affinity_scores": result.affinity_scores,
            "structure_file": str(result.predicted_structure) if result.predicted_structure else ""
        })
    
    # Add failed results separately
    for i, result in enumerate(failed_results):
        summary["failed_predictions"].append({
            "binder_id": f"failed_{i:03d}",
            "sequence": result.binder_sequence,
            "error_message": result.error_message,
            "prediction_success": False
        })
    
    with output_file.open('w') as f:
        json.dump(summary, f, indent=2)

def list_pdb_chains(pdb_path: Path, parser: PDBComplexParser):
    """List all chains and sequences in PDB file for user selection"""
    chains_info = parser.list_chains_in_pdb(pdb_path)
    
    print(f"\nChains found in {pdb_path.name}:")
    print("-" * 60)
    for chain_id, sequence in chains_info.items():
        print(f"Chain {chain_id}: {len(sequence):3d} residues - {sequence[:50]}{'...' if len(sequence) > 50 else ''}")
    print("-" * 60)
    
    return chains_info

def main():
    parser = argparse.ArgumentParser(
        description="Validate protein binder designs using Boltz",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  List chains:     python boltz_design_validator.py complex.pdb --list-chains
  Auto-detect:     python boltz_design_validator.py complex.pdb --auto-chains --sequences "MKLLVV..."
  Specify chains:  python boltz_design_validator.py complex.pdb --target-chain A --binder-chain B
  With sequences:  python boltz_design_validator.py complex.pdb --target-chain A --binder-chain B --sequences "MKLLVV..."
        """)
    parser.add_argument("reference_pdb", type=Path, help="Reference PDB file of target-binder complex")
    parser.add_argument("--target-chain", help="Target chain ID in reference PDB")
    parser.add_argument("--binder-chain", help="Binder chain ID in reference PDB") 
    
    # Sequence input options (mutually exclusive is clearer but not enforced for flexibility)
    parser.add_argument("--binder-sequences", type=Path, help="File with binder sequences (one per line)")
    parser.add_argument("--sequences", nargs="+", help="Binder sequences as command line arguments")
    
    parser.add_argument("--list-chains", action="store_true", help="List all chains in PDB and exit")
    parser.add_argument("--auto-chains", action="store_true", help="Auto-detect target/binder chains (largest two)")
    parser.add_argument("--output-dir", type=Path, default=Path("./validation_results"), help="Output directory")
    parser.add_argument("--work-dir", type=Path, default=Path("./boltz_work"), help="Working directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Validate PDB file exists
    if not args.reference_pdb.exists():
        print(f"Error: PDB file {args.reference_pdb} not found")
        return 1
    
    # Create single PDB parser instance for reuse
    pdb_parser = PDBComplexParser()
    
    # Handle chain listing mode
    if args.list_chains:
        try:
            list_pdb_chains(args.reference_pdb, pdb_parser)
            print("\nUse the chain IDs with --target-chain and --binder-chain arguments")
            return 0
        except Exception as e:
            print(f"Error reading PDB file: {e}")
            return 1
    
    # Handle chain detection
    if args.auto_chains or (not args.target_chain or not args.binder_chain):
        if args.auto_chains:
            # Auto-detect chains
            try:
                target_chain, binder_chain = pdb_parser.auto_detect_chains(args.reference_pdb)
                print(f"Auto-detected chains: Target={target_chain}, Binder={binder_chain}")
            except Exception as e:
                print(f"Error auto-detecting chains: {e}")
                return 1
        else:
            # Show available chains and prompt user
            try:
                chains_info = list_pdb_chains(args.reference_pdb, pdb_parser)
                
                if not args.target_chain:
                    print(f"\nPlease specify --target-chain from available chains: {list(chains_info.keys())}")
                if not args.binder_chain:
                    print(f"Please specify --binder-chain from available chains: {list(chains_info.keys())}")
                
                return 1
            except Exception as e:
                print(f"Error reading PDB file: {e}")
                return 1
    else:
        target_chain = args.target_chain
        binder_chain = args.binder_chain
    
    # CRITICAL: Validate that target and binder chains are different
    if target_chain == binder_chain:
        print(f"Error: Target chain ({target_chain}) and binder chain ({binder_chain}) must be different")
        return 1
    
    # Setup logging
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Load binder sequences
    if args.binder_sequences:
        try:
            binder_sequences = load_binder_sequences(args.binder_sequences)
        except Exception as e:
            print(f"Error loading binder sequences from {args.binder_sequences}: {e}")
            return 1
    elif args.sequences:
        binder_sequences = args.sequences
    else:
        # Auto-extract binder sequence from PDB if no sequences provided
        try:
            chains_info = pdb_parser.list_chains_in_pdb(args.reference_pdb)
            
            if binder_chain not in chains_info:
                print(f"Error: Binder chain '{binder_chain}' not found in PDB file")
                print(f"Available chains: {list(chains_info.keys())}")
                return 1
                
            binder_sequence = chains_info[binder_chain]
            binder_sequences = [binder_sequence]
            print(f"Auto-extracted binder sequence from chain {binder_chain}: {len(binder_sequence)} residues")
            
        except Exception as e:
            print(f"Error extracting binder sequence from PDB: {e}")
            return 1
    
    # Validate sequences
    if not binder_sequences:
        print("Error: No binder sequences provided")
        return 1
    
    for i, seq in enumerate(binder_sequences):
        if not seq or not seq.strip():
            print(f"Error: Empty sequence at position {i+1}")
            return 1
        if len(seq.strip()) < 5:
            print(f"Warning: Very short sequence at position {i+1}: {len(seq.strip())} residues")
    
    print(f"Validating {len(binder_sequences)} binder designs against {args.reference_pdb}")
    print(f"Target chain: {target_chain}, Binder chain: {binder_chain}")
    
    # Strip whitespace from sequences
    binder_sequences = [seq.strip() for seq in binder_sequences]
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Run validation pipeline
    pipeline = BoltzDesignPipeline(work_dir=args.work_dir)
    
    try:
        results = pipeline.validate_binder_designs(
            reference_pdb=args.reference_pdb,
            target_chain=target_chain,
            binder_chain=binder_chain,
            binder_sequences=binder_sequences
        )
        
        # Write results
        results_file = args.output_dir / "validation_results.json"
        write_results_summary(results, results_file)
        
        # Print summary
        print(f"\nValidation Complete!")
        print(f"Results written to: {results_file}")
        
        # Calculate and display summary stats
        successful_results = [r for r in results if r.prediction_success]
        failed_results = [r for r in results if not r.prediction_success]
        
        print(f"\nSummary:")
        print(f"Total binders processed: {len(results)}")
        print(f"Successful predictions: {len(successful_results)}")
        print(f"Failed predictions: {len(failed_results)}")
        
        if failed_results:
            print(f"\nFailed binders:")
            for result in failed_results:
                print(f"  - {result.binder_sequence[:20]}... Error: {result.error_message}")
        
        if successful_results:
            # Calculate metrics for successful results only
            plddt_scores = [r.confidence_scores.get('complex_plddt', 0) for r in successful_results if r.confidence_scores]
            rmsd_scores = [r.rmsd_metrics.get('complex_rmsd', 0) for r in successful_results if r.rmsd_metrics.get('complex_rmsd', -1) >= 0]
            rmsd_aligned_scores = [r.rmsd_metrics.get('complex_rmsd_aligned', 0) for r in successful_results if r.rmsd_metrics.get('complex_rmsd_aligned', -1) >= 0]
            
            if plddt_scores:
                avg_plddt = sum(plddt_scores) / len(plddt_scores)
                print(f"\nSuccess Metrics:")
                print(f"Average complex pLDDT: {avg_plddt:.3f}")
                
            if rmsd_scores:
                avg_rmsd = sum(rmsd_scores) / len(rmsd_scores)
                print(f"Average complex RMSD (unaligned): {avg_rmsd:.3f} Å")
                
            if rmsd_aligned_scores:
                avg_rmsd_aligned = sum(rmsd_aligned_scores) / len(rmsd_aligned_scores)
                print(f"Average complex RMSD (aligned): {avg_rmsd_aligned:.3f} Å")
            
            # Top candidates
            if plddt_scores:
                sorted_results = sorted(successful_results, key=lambda x: x.confidence_scores.get("complex_plddt", 0), reverse=True)
                print(f"\nTop 3 candidates by pLDDT:")
                for i, result in enumerate(sorted_results[:3]):
                    rmsd_unaligned = result.rmsd_metrics.get('complex_rmsd', 0)
                    rmsd_aligned = result.rmsd_metrics.get('complex_rmsd_aligned', 0)
                    affinity = result.affinity_scores.get('affinity_probability_binary', 0)
                    print(f"{i+1}. pLDDT: {result.confidence_scores.get('complex_plddt', 0):.3f}, "
                          f"RMSD: {rmsd_unaligned:.3f}/{rmsd_aligned:.3f} Å (unaligned/aligned), "
                          f"Affinity: {affinity:.3f}")
        else:
            print("No successful predictions generated")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())