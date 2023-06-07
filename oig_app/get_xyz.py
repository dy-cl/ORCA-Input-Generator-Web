from rdkit import Chem
from rdkit.Chem import AllChem, AddHs

def to_xyz(mol_string, mol_input_format):
    """
    Convert molecular input to .xyz format using RDKit.

    Parameters:
    - molecule_string (str): Molecular input string.
    - mol_input_format (str): Initial format type of the molecular string.
        Valid values: 'smiles', 'inchi'

    Returns:
    - None: If XYZ coordinates and InChI could not be generated.
    - xyz_content: Otherwise.
    """
    if mol_input_format == 'smiles':
        #Create RDKit molecule object from SMILES
        mol = Chem.MolFromSmiles(mol_string)
        if mol is None:
            print("Error: Failed to read the molecule from SMILES.")
            return None

    elif mol_input_format == 'inchi':
        #Create RDKit molecule object from InChI
        mol = Chem.MolFromInchi(mol_string)
        if mol is None:
            print("Error: Failed to read the molecule from InChI.")
            return None

    else:
        print("Error: Invalid input format. Valid formats are 'smiles' and 'inchi'.")
        return None

    #Add implicit hydrogens
    mol = Chem.AddHs(mol)

    #Generate 3D coordinates
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs = True, useBasicKnowledge = True, randomSeed = 10)

    #Optimize molecule
    AllChem.MMFFOptimizeMolecule(mol)

    #Write to XYZ format
    xyz_content = Chem.MolToXYZBlock(mol)

    #Remove the first line (atom count)
    lines = xyz_content.strip().split('\n')
    xyz_content = '\n'.join(lines[1:])

    return xyz_content