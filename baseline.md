# Boltz Codebase Analysis - Baseline Documentation

## Overview

Boltz is a family of deep learning models for biomolecular interaction prediction, with two main versions:
- **Boltz-1**: Initial model focused on structure prediction
- **Boltz-2**: Enhanced model with binding affinity prediction and template support

## Codebase Structure

### Core Architecture
```
src/boltz/
├── main.py                    # Main CLI interface and prediction workflow
├── data/                      # Data processing pipeline
│   ├── parse/                # Input parsing (YAML, FASTA, mmCIF)
│   ├── feature/              # Feature extraction and preparation
│   ├── tokenize/             # Structure tokenization for model input
│   ├── module/               # PyTorch Lightning data modules
│   └── types.py              # Core data structures and types
├── model/                     # Neural network architecture
│   ├── models/               # Boltz1 and Boltz2 model definitions
│   ├── modules/              # Model components (diffusion, transformers)
│   ├── layers/               # Building blocks (attention, etc.)
│   └── loss/                 # Loss functions
└── examples/                  # Example input files
```

### Key Data Structures

**Core Types:**
- `Structure/StructureV2`: 3D molecular structure representation
- `Record`: Metadata and chain information
- `Target`: Complete input specification with structure, sequences, templates
- `Tokenized`: Model-ready representation with tokens and bonds

**Template System:**
- `TemplateInfo`: Template alignment and mapping information
- Templates stored as `StructureV2` objects with tokenized representations
- Sequence alignment using BioPython PairwiseAligner

## Input Format Analysis

### YAML Input Capabilities
```yaml
sequences:          # Molecular sequences (protein, DNA, RNA, ligands)
constraints:        # Covalent bonds, pocket contacts, distance constraints  
templates:          # Structural templates (Boltz-2 only, proteins only)
properties:         # Binding affinity prediction
```

### Template Features (Boltz-2 Only)
- **Input**: mmCIF files with optional chain mapping
- **Processing**: Automatic sequence alignment for template-query mapping
- **Integration**: 3D coordinate priors and local geometry guidance
- **Limitations**: Protein chains only, requires reasonable sequence similarity

## Model Capabilities

### Boltz-1
- Structure prediction for proteins, nucleic acids, small molecules
- MSA-driven predictions using evolutionary information
- Diffusion-based sampling with confidence estimation

### Boltz-2 (Enhanced)
- All Boltz-1 capabilities plus:
- **Template Support**: Structural priors from known structures
- **Binding Affinity**: Two prediction modes (classification + regression)
- **Enhanced Architecture**: Improved transformer modules and training

### Key Prediction Parameters
```python
# Structure prediction settings
recycling_steps: 3           # Model iterations (default)
sampling_steps: 200          # Diffusion sampling steps  
diffusion_samples: 1         # Number of structure samples
step_scale: 1.5              # Boltz-2 sampling temperature

# Affinity prediction (Boltz-2)
affinity_probability_binary  # Binder vs non-binder (0-1)
affinity_pred_value         # log(IC50) in μM units
```

## Template System Deep Dive

### Template Processing Pipeline
1. **Parsing**: Load mmCIF files and extract structure information
2. **Alignment**: Sequence alignment between query and template chains
3. **Tokenization**: Convert template structure to model tokens
4. **Feature Generation**: Extract 3D coordinates, frames, validity masks
5. **Model Integration**: Process through dedicated transformer modules

### Template Data Flow
```
mmCIF → Structure → Alignment → TemplateInfo → Tokenized Templates → Model Features
```

### Template Limitations
- **Protein Only**: No support for DNA/RNA or ligand templates
- **Boltz-2 Only**: Not available in Boltz-1
- **Sequence Dependent**: Requires reasonable sequence similarity for alignment
- **Single Domain**: Best suited for single-domain template matching

## MSA Processing

### MSA Requirements
- **Default**: MSA required for protein chains
- **Auto-generation**: `--use_msa_server` flag enables mmseqs2 server
- **Custom MSA**: Support for pre-computed .a3m files
- **No MSA Mode**: `msa: empty` for single sequence (discouraged)

### MSA Data Structures
- `MSA`: Sequences, deletions, residue information
- Paired and unpaired sequence handling
- Taxonomy and clustering information

## Protein Design Integration Strategy

### Identified Integration Points

1. **Template-Based Target Reuse**
   - Generate high-quality target structure with full MSA
   - Use target structure as template for binder predictions
   - Set target `msa: empty` to skip MSA processing in screening

2. **De Novo Binder Handling**
   - Use `msa: empty` for designed binder sequences
   - Rely on protein language model embeddings
   - No MSA server calls for synthetic sequences

3. **Batch Processing Optimization**
   - Pre-process target template once
   - Reuse across multiple binder predictions
   - Parallel processing of binder variants

### Recommended Workflow

**Phase 1: Target Preparation**
```yaml
# target_reference.yaml
sequences:
  - protein:
      id: target
      sequence: "TARGET_SEQUENCE"
      msa: "path/to/target.a3m"  # Full evolutionary MSA
```

**Phase 2: Binder Validation**
```yaml
# binder_validation.yaml  
sequences:
  - protein:
      id: target
      sequence: "TARGET_SEQUENCE" 
      msa: empty                  # Skip MSA, use template
  - protein:
      id: binder
      sequence: "DESIGNED_BINDER"
      msa: empty                  # De novo sequence
templates:
  - cif: "target_reference.cif"
    chain_id: target
properties:
  - affinity:
      binder: binder
```

## Performance Considerations

### Computational Bottlenecks
1. **MSA Generation**: Most expensive step for novel sequences
2. **Template Processing**: One-time cost, reusable across predictions
3. **Diffusion Sampling**: Scales with sample count and steps

### Optimization Strategies
1. **Template Caching**: Pre-process and reuse target templates
2. **MSA Skipping**: Use `msa: empty` for known structures with templates
3. **Parallel Screening**: Process multiple binders simultaneously
4. **Reduced Sampling**: Single diffusion sample for initial screening

## Output Analysis

### Structure Quality Metrics
- `pLDDT`: Per-residue confidence (0-1, higher better)
- `PAE`: Predicted Aligned Error (Ångströms, lower better)  
- `ipTM`: Interface Template Modeling score
- `complex_plddt`: Overall complex confidence

### Binding Assessment  
- `affinity_probability_binary`: Binder classification (0-1)
- `affinity_pred_value`: Binding strength as log(IC50) μM
- Interface contact analysis and stability metrics

## Integration Recommendations

### Minimal Code Changes Required
- **No Core Modifications**: Use existing YAML interface
- **Wrapper Functions**: Build validation pipeline on top of CLI
- **Template Management**: Implement target structure reuse logic
- **Batch Processing**: Orchestrate multiple boltz predict calls

### Proposed Architecture
```python
class BoltzDesignValidator:
    """Protein design validation using Boltz"""
    
    def prepare_target_template(self, target_seq, msa_path):
        """Generate high-quality target structure"""
        
    def validate_binder(self, binder_seq, target_template):
        """Validate single binder design"""
        
    def batch_validate(self, binder_sequences):
        """Parallel validation of multiple binders"""
```

## Conclusion

Boltz provides a robust foundation for protein design validation through:
1. **Template System**: Efficient target structure reuse
2. **MSA Flexibility**: Handle both evolutionary and de novo sequences  
3. **Affinity Prediction**: Quantitative binding assessment
4. **Batch Processing**: Scalable screening workflows

The template-based approach enables AF2-style optimization by using high-quality target structures as priors while maintaining prediction speed for designed binder evaluation.