# Implementation Plan: Boltz Protein Design Validator

## Simple Pipeline Approach

### Core Concept: Single Pipeline Script

Input: PDB file of target-binder complex + designed binder sequences  
Output: Predictions with confidence scores and RMSD metrics against reference

### Step 1: PDB Structure Handler (0.5 days)
Location: `boltz_design_pipeline/structure_utils.py`

```python
from pathlib import Path
from typing import Tuple, Dict, List
from dataclasses import dataclass
import numpy as np
from Bio import PDB
from Bio.PDB import PDBIO, Select

@dataclass
class ComplexInfo:
    """Information extracted from input PDB complex"""
    target_chain_id: str
    binder_chain_id: str
    target_sequence: str
    binder_sequence: str
    target_coords: np.ndarray
    binder_coords: np.ndarray
    complex_coords: np.ndarray
    pdb_path: Path

@dataclass  
class ValidationResult:
    """Complete validation result with metrics"""
    binder_sequence: str
    predicted_structure: Path
    confidence_scores: Dict[str, float]
    rmsd_metrics: Dict[str, float]
    affinity_scores: Dict[str, float]

class PDBComplexParser:
    """Parse PDB files to extract target-binder information"""
    
    def __init__(self):
        self.parser = PDB.PDBParser(QUIET=True)
        
    def list_chains_in_pdb(self, pdb_path: Path) -> Dict[str, str]:
        """List all chains in PDB with their sequences for chain selection"""
        structure = self.parser.get_structure("complex", pdb_path)
        chains_info = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.id
                sequence, _ = self._extract_chain_info(structure, chain_id)
                chains_info[chain_id] = sequence
                
        return chains_info
    
    def parse_complex(self, pdb_path: Path, target_chain: str, binder_chain: str) -> ComplexInfo:
        """Extract sequences and coordinates from PDB complex"""
        structure = self.parser.get_structure("complex", pdb_path)
        
        # Extract chain sequences and coordinates
        target_seq, target_coords = self._extract_chain_info(structure, target_chain)
        binder_seq, binder_coords = self._extract_chain_info(structure, binder_chain)
        
        # Get full complex coordinates
        complex_coords = self._get_complex_coords(structure, [target_chain, binder_chain])
        
        return ComplexInfo(
            target_chain_id=target_chain,
            binder_chain_id=binder_chain,
            target_sequence=target_seq,
            binder_sequence=binder_seq,
            target_coords=target_coords,
            binder_coords=binder_coords,
            complex_coords=complex_coords,
            pdb_path=pdb_path
        )
    
    def auto_detect_chains(self, pdb_path: Path) -> Tuple[str, str]:
        """Auto-detect target and binder chains (largest two chains)"""
        chains_info = self.list_chains_in_pdb(pdb_path)
        
        # Sort chains by sequence length (descending)
        sorted_chains = sorted(chains_info.items(), key=lambda x: len(x[1]), reverse=True)
        
        if len(sorted_chains) < 2:
            raise ValueError(f"PDB file must contain at least 2 protein chains, found {len(sorted_chains)}")
        
        target_chain = sorted_chains[0][0]  # Longest chain as target
        binder_chain = sorted_chains[1][0]  # Second longest as binder
        
        return target_chain, binder_chain
    
    def extract_target_structure(self, complex_info: ComplexInfo, output_path: Path) -> Path:
        """Extract target chain as separate PDB for template use"""
        structure = self.parser.get_structure("complex", complex_info.pdb_path)
        
        class TargetSelector(Select):
            def accept_chain(self, chain):
                return chain.id == complex_info.target_chain_id
                
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(output_path), TargetSelector())
        return output_path

class RMSDCalculator:
    """Calculate RMSD between predicted and reference structures"""
    
    @staticmethod
    def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """Calculate RMSD between two coordinate arrays"""
        if len(coords1) != len(coords2):
            raise ValueError("Coordinate arrays must have same length")
        return np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))
    
    def calculate_complex_rmsd(self, pred_structure: Path, ref_complex: ComplexInfo) -> Dict[str, float]:
        """Calculate RMSD metrics for predicted vs reference complex"""
        pred_coords = self._extract_coords_from_cif(pred_structure)
        
        return {
            "complex_rmsd": self.calculate_rmsd(pred_coords, ref_complex.complex_coords),
            "target_rmsd": self.calculate_rmsd(pred_coords[:len(ref_complex.target_coords)], ref_complex.target_coords),
            "binder_rmsd": self.calculate_rmsd(pred_coords[len(ref_complex.target_coords):], ref_complex.binder_coords)
        }
```

### Step 2: Simple Boltz Pipeline (1 day)
Location: `boltz_design_pipeline/boltz_pipeline.py`

```python
from pathlib import Path
from typing import List, Dict
import subprocess
import yaml
import json
import logging
from .structure_utils import ComplexInfo, ValidationResult, PDBComplexParser, RMSDCalculator

class BoltzDesignPipeline:
    """Simple pipeline for protein design validation"""
    
    def __init__(self, work_dir: Path = Path("./boltz_work")):
        self.work_dir = work_dir
        self.work_dir.mkdir(exist_ok=True)
        self.pdb_parser = PDBComplexParser()
        self.rmsd_calc = RMSDCalculator()
        
    def validate_binder_designs(
        self, 
        reference_pdb: Path,
        target_chain: str,
        binder_chain: str, 
        binder_sequences: List[str]
    ) -> List[ValidationResult]:
        """Complete pipeline: PDB input -> Boltz predictions -> RMSD metrics"""
        
        # Step 1: Parse reference complex
        logging.info(f"Parsing reference complex: {reference_pdb}")
        ref_complex = self.pdb_parser.parse_complex(reference_pdb, target_chain, binder_chain)
        
        # Step 2: Extract target structure for template
        target_template = self.work_dir / "target_template.pdb" 
        self.pdb_parser.extract_target_structure(ref_complex, target_template)
        
        # Step 3: Convert PDB to mmCIF for Boltz template
        template_cif = self._convert_pdb_to_cif(target_template)
        
        # Step 4: Run predictions for each binder
        results = []
        for i, binder_seq in enumerate(binder_sequences):
            logging.info(f"Processing binder {i+1}/{len(binder_sequences)}")
            
            # Create Boltz config
            config = self._create_binder_config(
                ref_complex.target_sequence, 
                binder_seq, 
                template_cif
            )
            
            # Run prediction
            pred_structure = self._run_boltz_prediction(config, f"binder_{i:03d}")
            
            # Calculate metrics
            confidence_scores = self._extract_confidence_scores(f"binder_{i:03d}")
            affinity_scores = self._extract_affinity_scores(f"binder_{i:03d}")
            rmsd_metrics = self.rmsd_calc.calculate_complex_rmsd(pred_structure, ref_complex)
            
            results.append(ValidationResult(
                binder_sequence=binder_seq,
                predicted_structure=pred_structure,
                confidence_scores=confidence_scores,
                rmsd_metrics=rmsd_metrics,
                affinity_scores=affinity_scores
            ))
            
        return results
    
    def _create_binder_config(self, target_seq: str, binder_seq: str, template_path: Path) -> Dict:
        """Create Boltz YAML configuration for binder prediction"""
        return {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "target",
                        "sequence": target_seq,
                        "msa": "empty"  # Use template instead
                    }
                },
                {
                    "protein": {
                        "id": "binder",
                        "sequence": binder_seq, 
                        "msa": "empty"  # De novo binder
                    }
                }
            ],
            "templates": [
                {
                    "cif": str(template_path),
                    "chain_id": "target"
                }
            ],
            "properties": [
                {
                    "affinity": {
                        "binder": "binder"
                    }
                }
            ]
        }
    
    def _run_boltz_prediction(self, config: Dict, output_name: str) -> Path:
        """Run Boltz prediction and return structure path"""
        config_path = self.work_dir / f"{output_name}.yaml"
        
        with config_path.open('w') as f:
            yaml.dump(config, f)
        
        cmd = [
            "boltz", "predict", str(config_path),
            "--out_dir", str(self.work_dir),
            "--model", "boltz2",
            "--diffusion_samples", "1"  # Fast for screening
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Find generated structure
        pred_dir = self.work_dir / f"boltz_results_{output_name}" / "predictions" / output_name
        structure_files = list(pred_dir.glob("*_model_0.cif"))
        
        if not structure_files:
            raise RuntimeError(f"No structure generated for {output_name}")
            
        return structure_files[0]
```

### Step 3: Main Pipeline Script (0.5 days)
Location: `boltz_design_validator.py` (Main executable script)

```python
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
from boltz_design_pipeline.boltz_pipeline import BoltzDesignPipeline
from boltz_design_pipeline.structure_utils import ValidationResult

def load_binder_sequences(sequences_file: Path) -> List[str]:
    """Load binder sequences from file (one per line)"""
    with sequences_file.open() as f:
        sequences = [line.strip() for line in f if line.strip()]
    return sequences

def write_results_summary(results: List[ValidationResult], output_file: Path):
    """Write results to JSON file with metrics summary"""
    summary = {
        "validation_results": [],
        "summary_stats": {
            "total_binders": len(results),
            "avg_complex_plddt": sum(r.confidence_scores.get("complex_plddt", 0) for r in results) / len(results),
            "avg_complex_rmsd": sum(r.rmsd_metrics.get("complex_rmsd", 0) for r in results) / len(results),
            "high_confidence_count": sum(1 for r in results if r.confidence_scores.get("complex_plddt", 0) > 0.8),
            "good_rmsd_count": sum(1 for r in results if r.rmsd_metrics.get("complex_rmsd", 999) < 2.0)
        }
    }
    
    for i, result in enumerate(results):
        summary["validation_results"].append({
            "binder_id": f"binder_{i:03d}",
            "sequence": result.binder_sequence,
            "confidence_scores": result.confidence_scores,
            "rmsd_metrics": result.rmsd_metrics,
            "affinity_scores": result.affinity_scores,
            "structure_file": str(result.predicted_structure)
        })
    
    with output_file.open('w') as f:
        json.dump(summary, f, indent=2)

def list_pdb_chains(pdb_path: Path):
    """List all chains and sequences in PDB file for user selection"""
    from boltz_design_pipeline.structure_utils import PDBComplexParser
    
    parser = PDBComplexParser()
    chains_info = parser.list_chains_in_pdb(pdb_path)
    
    print(f"\nChains found in {pdb_path.name}:")
    print("-" * 60)
    for chain_id, sequence in chains_info.items():
        print(f"Chain {chain_id}: {len(sequence):3d} residues - {sequence[:50]}{'...' if len(sequence) > 50 else ''}")
    print("-" * 60)
    
    return chains_info

def main():
    parser = argparse.ArgumentParser(description="Validate protein binder designs using Boltz")
    parser.add_argument("reference_pdb", type=Path, help="Reference PDB file of target-binder complex")
    parser.add_argument("--target-chain", help="Target chain ID in reference PDB (auto-detect if not provided)")
    parser.add_argument("--binder-chain", help="Binder chain ID in reference PDB (auto-detect if not provided)") 
    parser.add_argument("--binder-sequences", type=Path, help="File with binder sequences (one per line)")
    parser.add_argument("--sequences", nargs="+", help="Binder sequences as command line arguments")
    parser.add_argument("--list-chains", action="store_true", help="List all chains in PDB and exit")
    parser.add_argument("--auto-chains", action="store_true", help="Auto-detect target/binder chains (largest two)")
    parser.add_argument("--output-dir", type=Path, default=Path("./validation_results"), help="Output directory")
    parser.add_argument("--work-dir", type=Path, default=Path("./boltz_work"), help="Working directory")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    
    args = parser.parse_args()
    
    # Handle chain listing mode
    if args.list_chains:
        list_pdb_chains(args.reference_pdb)
        print("\nUse the chain IDs with --target-chain and --binder-chain arguments")
        return 0
    
    # Handle chain detection
    from boltz_design_pipeline.structure_utils import PDBComplexParser
    
    if args.auto_chains or (not args.target_chain or not args.binder_chain):
        pdb_parser = PDBComplexParser()
        
        if args.auto_chains:
            # Auto-detect chains
            target_chain, binder_chain = pdb_parser.auto_detect_chains(args.reference_pdb)
            print(f"Auto-detected chains: Target={target_chain}, Binder={binder_chain}")
        else:
            # Show available chains and prompt user
            chains_info = list_pdb_chains(args.reference_pdb)
            
            if not args.target_chain:
                print(f"\nPlease specify --target-chain from available chains: {list(chains_info.keys())}")
            if not args.binder_chain:
                print(f"Please specify --binder-chain from available chains: {list(chains_info.keys())}")
            
            return 1
    else:
        target_chain = args.target_chain
        binder_chain = args.binder_chain
    
    # Setup logging
    log_level = logging.INFO if args.verbose else logging.WARNING
    logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
    
    # Load binder sequences
    if args.binder_sequences:
        binder_sequences = load_binder_sequences(args.binder_sequences)
    elif args.sequences:
        binder_sequences = args.sequences
    else:
        parser.error("Must provide either --binder-sequences file or --sequences")
    
    print(f"Validating {len(binder_sequences)} binder designs against {args.reference_pdb}")
    print(f"Target chain: {target_chain}, Binder chain: {binder_chain}")
    
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
        print(f"Average complex pLDDT: {sum(r.confidence_scores.get('complex_plddt', 0) for r in results) / len(results):.3f}")
        print(f"Average complex RMSD: {sum(r.rmsd_metrics.get('complex_rmsd', 0) for r in results) / len(results):.3f} Å")
        
        # Top candidates
        sorted_results = sorted(results, key=lambda x: x.confidence_scores.get("complex_plddt", 0), reverse=True)
        print(f"\nTop 3 candidates by pLDDT:")
        for i, result in enumerate(sorted_results[:3]):
            print(f"{i+1}. pLDDT: {result.confidence_scores.get('complex_plddt', 0):.3f}, "
                  f"RMSD: {result.rmsd_metrics.get('complex_rmsd', 0):.3f} Å, "
                  f"Affinity: {result.affinity_scores.get('affinity_probability_binary', 0):.3f}")
        
    except Exception as e:
        logging.error(f"Pipeline failed: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())
```

## Implementation Details

### Step 4: Missing Implementation Methods (0.5 days)

Key methods to implement in the pipeline classes:

```python
# In PDBComplexParser class
def _extract_chain_info(self, structure, chain_id: str) -> Tuple[str, np.ndarray]:
    """Extract sequence and coordinates from chain"""
    chain = structure[0][chain_id]
    sequence = ""
    coords = []
    
    for residue in chain:
        if PDB.is_aa(residue):
            sequence += PDB.Polypeptide.three_to_one(residue.get_resname())
            if 'CA' in residue:
                coords.append(residue['CA'].get_coord())
    
    return sequence, np.array(coords)

def _get_complex_coords(self, structure, chain_ids: List[str]) -> np.ndarray:
    """Get coordinates for entire complex"""
    all_coords = []
    for chain_id in chain_ids:
        _, coords = self._extract_chain_info(structure, chain_id)
        all_coords.append(coords)
    return np.concatenate(all_coords)

# In BoltzDesignPipeline class  
def _convert_pdb_to_cif(self, pdb_path: Path) -> Path:
    """Convert PDB to mmCIF format for Boltz template"""
    cif_path = pdb_path.with_suffix('.cif')
    # Use existing tools like BioPython or external converters
    from Bio.PDB import MMCIFIO
    structure = self.pdb_parser.parser.get_structure("temp", pdb_path)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))
    return cif_path

def _extract_confidence_scores(self, output_name: str) -> Dict[str, float]:
    """Extract confidence metrics from Boltz output"""
    pred_dir = self.work_dir / f"boltz_results_{output_name}" / "predictions" / output_name
    confidence_files = list(pred_dir.glob("confidence_*.json"))
    
    if confidence_files:
        with confidence_files[0].open() as f:
            data = json.load(f)
        return {
            "complex_plddt": data.get("complex_plddt", 0.0),
            "interface_plddt": data.get("complex_iplddt", 0.0),
            "ptm": data.get("ptm", 0.0),
            "iptm": data.get("iptm", 0.0)
        }
    return {}

def _extract_affinity_scores(self, output_name: str) -> Dict[str, float]:
    """Extract affinity predictions from Boltz output"""
    pred_dir = self.work_dir / f"boltz_results_{output_name}" / "predictions" / output_name
    affinity_files = list(pred_dir.glob("affinity_*.json"))
    
    if affinity_files:
        with affinity_files[0].open() as f:
            data = json.load(f)
        return {
            "affinity_probability_binary": data.get("affinity_probability_binary", 0.0),
            "affinity_pred_value": data.get("affinity_pred_value", 0.0)
        }
    return {}

# In RMSDCalculator class
def _extract_coords_from_cif(self, cif_path: Path) -> np.ndarray:
    """Extract coordinates from CIF file"""
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("pred", cif_path)
    
    coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if 'CA' in residue:
                    coords.append(residue['CA'].get_coord())
    
    return np.array(coords)
```

## Simple Usage Examples

### Sequence Extraction and Chain Detection
```bash
# List all chains and sequences in PDB file
python boltz_design_validator.py reference_complex.pdb --list-chains

# Output:
# Chains found in reference_complex.pdb:
# ------------------------------------------------------------
# Chain A: 150 residues - MKLLVGVVEQWKRQLEDGRTLSADIQQWLAKQSQELKRTEEQ...
# Chain B:  85 residues - AKLSILPWGHCDEFGHIKLMNPQRSTVWXYZ...
# Chain C:  12 residues - MGSSHHHHHHSS
# ------------------------------------------------------------

# Auto-detect target/binder chains (largest two)
python boltz_design_validator.py reference_complex.pdb \
    --auto-chains \
    --sequences "MKLLVGVVEQWKRQ" "AKLSILPWGHC"

# Manual chain specification (if auto-detection incorrect)
python boltz_design_validator.py reference_complex.pdb \
    --target-chain A --binder-chain B \
    --sequences "MKLLVGVVEQWKRQ" "AKLSILPWGHC"
```

### Basic Usage
```bash
# Validate single binder with automatic chain detection
python boltz_design_validator.py reference_complex.pdb \
    --auto-chains \
    --sequences "MKLLVGVVEQWKRQ" "AKLSILPWGHC"

# Validate from file with manual chain selection
python boltz_design_validator.py reference_complex.pdb \
    --target-chain A --binder-chain B \
    --binder-sequences designed_binders.txt

# If chains not specified, tool will prompt
python boltz_design_validator.py reference_complex.pdb \
    --sequences "MKLLVGVVEQWKRQ"
# Output: Lists available chains and asks for --target-chain and --binder-chain
```

### Expected Output
```json
{
  "validation_results": [
    {
      "binder_id": "binder_000",
      "sequence": "MKLLVGVVEQWKRQ",
      "confidence_scores": {
        "complex_plddt": 0.85,
        "interface_plddt": 0.78,
        "ptm": 0.82,
        "iptm": 0.79
      },
      "rmsd_metrics": {
        "complex_rmsd": 1.2,
        "target_rmsd": 0.8,
        "binder_rmsd": 2.1
      },
      "affinity_scores": {
        "affinity_probability_binary": 0.92,
        "affinity_pred_value": -1.5
      }
    }
  ],
  "summary_stats": {
    "total_binders": 10,
    "avg_complex_plddt": 0.75,
    "avg_complex_rmsd": 2.3,
    "high_confidence_count": 6,
    "good_rmsd_count": 4
  }
}
```

## File Structure (Simplified)
```
boltz_design_pipeline/
├── __init__.py
├── structure_utils.py      # PDB parsing, RMSD calculation
└── boltz_pipeline.py       # Main pipeline logic

boltz_design_validator.py   # Main executable script
sequences.txt               # Input binder sequences
validation_results.json     # Output results
```

## Implementation Timeline (Simplified)
- **Day 1**: Structure utilities and PDB handling
- **Day 2**: Boltz pipeline integration
- **Day 3**: Main script and testing

## Key Design Decisions

1. **Single PDB Input**: Takes reference complex, extracts target automatically
2. **Template Reuse**: Uses target from reference as template for all predictions  
3. **RMSD Benchmarking**: Direct comparison against reference structure
4. **Simple Interface**: Command-line script with JSON output
5. **No Boltz Changes**: Pure wrapper approach

This simplified approach focuses on the core requirement: validating designed binders against a reference complex structure with confidence and RMSD metrics.