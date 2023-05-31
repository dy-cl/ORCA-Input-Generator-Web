from rdkit import Chem
from rdkit.Chem import AllChem, AddHs

def to_xyz(mol_string, mol_input_format):
    """
    Convert molecular input to .xyz format using RDKit.

    Parameters:
    - molecule_string (str): Molecular input string.
    - mol_input_format (str): Initial format type of the molecular string.

    Returns:
    - None: If XYZ coordinates could not be generated.
    - xyz_content: Otherwise.
    """
    # Create RDKit molecule object
    mol = Chem.MolFromSmiles(mol_string)

    if mol is None:
        print("Error: Failed to read the molecule.")
        return None

    # Add implicit hydrogens
    mol = Chem.AddHs(mol)


    # Generate 3D coordinates
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)

    # Write to XYZ format
    xyz_content = Chem.MolToXYZBlock(mol)

    # Remove the first line (atom count)
    lines = xyz_content.strip().split('\n')
    xyz_content = '\n'.join(lines[1:])

    return xyz_content