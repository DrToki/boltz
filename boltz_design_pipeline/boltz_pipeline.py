from pathlib import Path
from typing import List, Dict
import subprocess
import yaml
import json
import logging
from Bio.PDB import MMCIFIO
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
            
            try:
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
                    affinity_scores=affinity_scores,
                    prediction_success=True,
                    error_message=""
                ))
                
            except Exception as e:
                logging.error(f"Failed to process binder {i}: {e}")
                # Add failed result with clear failure markers
                results.append(ValidationResult(
                    binder_sequence=binder_seq,
                    predicted_structure=Path(""),
                    confidence_scores={},
                    rmsd_metrics={},
                    affinity_scores={},
                    prediction_success=False,
                    error_message=str(e)
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
        
        logging.info(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        
        # Find generated structure
        pred_dir = self.work_dir / f"boltz_results_{output_name}" / "predictions" / output_name
        structure_files = list(pred_dir.glob("*_model_0.cif"))
        
        if not structure_files:
            raise RuntimeError(f"No structure generated for {output_name}")
            
        return structure_files[0]
    
    def _convert_pdb_to_cif(self, pdb_path: Path) -> Path:
        """Convert PDB to mmCIF format for Boltz template"""
        cif_path = pdb_path.with_suffix('.cif')
        
        # Use BioPython to convert
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
        return {
            "complex_plddt": 0.0,
            "interface_plddt": 0.0,
            "ptm": 0.0,
            "iptm": 0.0
        }

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
        return {
            "affinity_probability_binary": 0.0,
            "affinity_pred_value": 0.0
        }