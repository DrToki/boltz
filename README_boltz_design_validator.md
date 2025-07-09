# Boltz Design Validator

A Python pipeline for validating protein binder designs against reference complex structures using the Boltz prediction engine.

## Purpose

This tool validates designed protein binders by:
- Extracting target and binder sequences from reference PDB files
- Generating structure predictions using Boltz
- Calculating binding affinity, confidence scores, and RMSD metrics
- Providing comprehensive validation reports

## Requirements

- Python 3.8+
- Boltz prediction engine installed
- Required Python packages:
  - `boltz_design_pipeline` (custom module)
  - `argparse`, `pathlib`, `json`, `logging`

## Installation

```bash
# Ensure Boltz is installed and accessible
# Install required dependencies
pip install -r requirements_design_validator.txt
```

## Usage

### Basic Usage Modes

#### 1. List Available Chains
```bash
python boltz_design_validator.py reference.pdb --list-chains
```

#### 2. Auto-detect Chains with Custom Sequences
```bash
python boltz_design_validator.py reference.pdb --auto-chains --sequences "MKLLVVTG..." "AEQKLISEED..."
```

#### 3. Specify Chains with Sequence File
```bash
python boltz_design_validator.py reference.pdb --target-chain A --binder-chain B --binder-sequences sequences.txt
```

#### 4. Use PDB Binder Sequence (NEW)
```bash
python boltz_design_validator.py reference.pdb --target-chain A --binder-chain B
```
*Automatically extracts binder sequence from chain B in the PDB file*

### Command Line Arguments

| Argument | Description |
|----------|-------------|
| `reference_pdb` | Reference PDB file of target-binder complex |
| `--target-chain` | Target chain ID in reference PDB |
| `--binder-chain` | Binder chain ID in reference PDB |
| `--binder-sequences` | File with binder sequences (one per line) |
| `--sequences` | Binder sequences as command line arguments |
| `--list-chains` | List all chains in PDB and exit |
| `--auto-chains` | Auto-detect target/binder chains (largest two) |
| `--output-dir` | Output directory (default: ./validation_results) |
| `--work-dir` | Working directory (default: ./boltz_work) |
| `--verbose` | Verbose logging |

## Output

The pipeline generates:

### 1. JSON Results File (`validation_results.json`)
```json
{
  "validation_results": [
    {
      "binder_id": "binder_001",
      "sequence": "MKLLVVTG...",
      "confidence_scores": {
        "complex_plddt": 0.85
      },
      "rmsd_metrics": {
        "complex_rmsd": 1.2,
        "complex_rmsd_aligned": 0.8
      },
      "affinity_scores": {
        "affinity_probability_binary": 0.92
      },
      "structure_file": "path/to/predicted_structure.pdb"
    }
  ],
  "failed_predictions": [...],
  "summary_stats": {
    "total_binders": 10,
    "successful_predictions": 8,
    "failed_predictions": 2,
    "avg_complex_plddt": 0.82,
    "avg_complex_rmsd": 1.5,
    "high_confidence_count": 6,
    "good_rmsd_count": 7
  }
}
```

### 2. Console Summary
- Success/failure counts
- Average confidence scores and RMSD values
- Top candidates ranked by pLDDT
- Failed predictions with error messages

## Examples

### Validate Multiple Designed Binders
```bash
# Create sequences.txt with one sequence per line
echo "MKLLVVTGDQYADSVKGRFTISRDYSKNTLYLQMNSLR" > sequences.txt
echo "AEQKLISEEDLKAVEEAHSSLMAKMEQLLMGWDDEDE" >> sequences.txt

python boltz_design_validator.py complex.pdb --target-chain A --binder-chain B --binder-sequences sequences.txt
```

### Quick Validation with Auto-detection
```bash
python boltz_design_validator.py complex.pdb --auto-chains --sequences "MKLLVVTGDQYADSVKGRFTISRDYSKNTLYLQMNSLR"
```

### Validate Against Reference Binder
```bash
python boltz_design_validator.py complex.pdb --target-chain A --binder-chain B
```

## Metrics Explained

- **pLDDT**: Confidence score (0-1), higher is better
- **RMSD**: Root Mean Square Deviation in Ångstroms, lower is better
- **Affinity Probability**: Predicted binding probability (0-1), higher is better
- **High Confidence**: pLDDT > 0.8
- **Good RMSD**: RMSD < 2.0 Å

## Troubleshooting

### Common Issues
1. **Chain not found**: Use `--list-chains` to see available chains
2. **No sequences provided**: Specify sequences via file, command line, or use PDB extraction
3. **Prediction failures**: Check sequence validity and Boltz installation
4. **Memory issues**: Reduce number of sequences or increase system memory

### Error Messages
- `Chain X not found in PDB file`: Verify chain ID exists in structure
- `Error auto-detecting chains`: PDB may have unusual structure
- `No binder sequences provided`: Must specify sequences or use auto-extraction

## Performance Notes

- Processing time depends on sequence length and Boltz model complexity
- Large sequence files may require significant memory
- Results are cached in the work directory for faster re-runs