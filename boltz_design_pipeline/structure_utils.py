from pathlib import Path
from typing import Tuple, Dict, List
from dataclasses import dataclass
import numpy as np
from Bio import PDB
from Bio.PDB import PDBIO, Select, MMCIFParser

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
    prediction_success: bool = True
    error_message: str = ""

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
                if sequence:  # Only include chains with valid sequences
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
    
    def _extract_chain_info(self, structure, chain_id: str) -> Tuple[str, np.ndarray]:
        """Extract sequence and coordinates from chain"""
        try:
            chain = structure[0][chain_id]
        except KeyError:
            raise ValueError(f"Chain {chain_id} not found in structure")
            
        sequence = ""
        coords = []
        
        for residue in chain:
            if PDB.is_aa(residue):
                try:
                    sequence += PDB.Polypeptide.three_to_one(residue.get_resname())
                    if 'CA' in residue:
                        coords.append(residue['CA'].get_coord())
                except KeyError:
                    # Skip non-standard residues that can't be converted
                    continue
        
        return sequence, np.array(coords)

    def _get_complex_coords(self, structure, chain_ids: List[str]) -> np.ndarray:
        """Get coordinates for entire complex"""
        all_coords = []
        for chain_id in chain_ids:
            _, coords = self._extract_chain_info(structure, chain_id)
            all_coords.append(coords)
        return np.concatenate(all_coords)

class RMSDCalculator:
    """Calculate RMSD between predicted and reference structures"""
    
    @staticmethod
    def calculate_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
        """Calculate RMSD between two coordinate arrays"""
        if len(coords1) != len(coords2):
            raise ValueError(f"Coordinate arrays must have same length: {len(coords1)} vs {len(coords2)}")
        return np.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))
    
    @staticmethod
    def kabsch_align(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, float]:
        """
        Align coords1 to coords2 using Kabsch algorithm and return aligned coords1 and RMSD
        
        Parameters
        ----------
        coords1 : np.ndarray
            Coordinates to be aligned (N, 3)
        coords2 : np.ndarray  
            Reference coordinates (N, 3)
            
        Returns
        -------
        Tuple[np.ndarray, float]
            Aligned coordinates and RMSD after alignment
        """
        if len(coords1) != len(coords2):
            raise ValueError(f"Coordinate arrays must have same length: {len(coords1)} vs {len(coords2)}")
        
        # Center coordinates
        centroid1 = np.mean(coords1, axis=0)
        centroid2 = np.mean(coords2, axis=0)
        
        coords1_centered = coords1 - centroid1
        coords2_centered = coords2 - centroid2
        
        # Compute covariance matrix
        H = coords1_centered.T @ coords2_centered
        
        # SVD decomposition
        U, S, Vt = np.linalg.svd(H)
        
        # Compute rotation matrix
        R = Vt.T @ U.T
        
        # Ensure proper rotation (det(R) = 1)
        if np.linalg.det(R) < 0:
            Vt[-1, :] *= -1
            R = Vt.T @ U.T
        
        # Apply rotation and translation
        coords1_aligned = (coords1_centered @ R.T) + centroid2
        
        # Calculate RMSD after alignment
        rmsd = np.sqrt(np.mean(np.sum((coords1_aligned - coords2)**2, axis=1)))
        
        return coords1_aligned, rmsd
    
    def calculate_complex_rmsd(self, pred_structure: Path, ref_complex: ComplexInfo) -> Dict[str, float]:
        """Calculate RMSD metrics for predicted vs reference complex with structural alignment"""
        pred_coords = self._extract_coords_from_cif(pred_structure)
        
        # Calculate RMSD for different parts of the complex
        n_target = len(ref_complex.target_coords)
        n_binder = len(ref_complex.binder_coords)
        
        if len(pred_coords) != (n_target + n_binder):
            # Handle length mismatch - align what we can
            min_complex_len = min(len(pred_coords), len(ref_complex.complex_coords))
            
            if min_complex_len < 3:
                # Cannot align with fewer than 3 points
                return {
                    "complex_rmsd": -1.0,
                    "complex_rmsd_aligned": -1.0,
                    "target_rmsd": -1.0,
                    "target_rmsd_aligned": -1.0,
                    "binder_rmsd": -1.0,
                    "binder_rmsd_aligned": -1.0
                }
            
            # Align available coordinates
            pred_coords_trunc = pred_coords[:min_complex_len]
            ref_coords_trunc = ref_complex.complex_coords[:min_complex_len]
            
            _, complex_rmsd_aligned = self.kabsch_align(pred_coords_trunc, ref_coords_trunc)
            complex_rmsd_unaligned = self.calculate_rmsd(pred_coords_trunc, ref_coords_trunc)
            
            return {
                "complex_rmsd": complex_rmsd_unaligned,
                "complex_rmsd_aligned": complex_rmsd_aligned,
                "target_rmsd": -1.0,
                "target_rmsd_aligned": -1.0,
                "binder_rmsd": -1.0,
                "binder_rmsd_aligned": -1.0
            }
        
        # Full complex alignment
        _, complex_rmsd_aligned = self.kabsch_align(pred_coords, ref_complex.complex_coords)
        complex_rmsd_unaligned = self.calculate_rmsd(pred_coords, ref_complex.complex_coords)
        
        # Target chain alignment
        pred_target = pred_coords[:n_target]
        if len(pred_target) >= 3:
            _, target_rmsd_aligned = self.kabsch_align(pred_target, ref_complex.target_coords)
            target_rmsd_unaligned = self.calculate_rmsd(pred_target, ref_complex.target_coords)
        else:
            target_rmsd_aligned = target_rmsd_unaligned = -1.0
        
        # Binder chain alignment  
        pred_binder = pred_coords[n_target:]
        if len(pred_binder) >= 3:
            _, binder_rmsd_aligned = self.kabsch_align(pred_binder, ref_complex.binder_coords)
            binder_rmsd_unaligned = self.calculate_rmsd(pred_binder, ref_complex.binder_coords)
        else:
            binder_rmsd_aligned = binder_rmsd_unaligned = -1.0
        
        return {
            "complex_rmsd": complex_rmsd_unaligned,
            "complex_rmsd_aligned": complex_rmsd_aligned,
            "target_rmsd": target_rmsd_unaligned,
            "target_rmsd_aligned": target_rmsd_aligned,
            "binder_rmsd": binder_rmsd_unaligned,
            "binder_rmsd_aligned": binder_rmsd_aligned
        }
    
    def _extract_coords_from_cif(self, cif_path: Path) -> np.ndarray:
        """Extract coordinates from CIF file"""
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("pred", cif_path)
        
        coords = []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if 'CA' in residue:
                        coords.append(residue['CA'].get_coord())
        
        return np.array(coords)