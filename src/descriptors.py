from rdkit import Chem
from rdkit.Chem import Descriptors

def compute_descriptors(smiles):
    """
    Computes basic chemical descriptors using RDKit.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        dict: Dictionary of descriptors (MolWt, LogP, etc.)
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return {
            "MolWt": None,
            "LogP": None,
            "NumHDonors": None,
            "NumHAcceptors": None
        }

    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol)
    }
