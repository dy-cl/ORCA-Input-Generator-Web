import io
import base64
from ase import Atoms
from ase.io import write
from ase.visualize import view

def draw_mol(xyz_data):
    xyz_data = xyz_data.strip()

    lines = xyz_data.split('\n')

    symbols = []
    positions = []
    for line in lines:
        symbol, x, y, z = line.split()
        symbols.append(symbol)
        positions.append((float(x), float(y), float(z)))

    atoms = Atoms(symbols=symbols, positions=positions)

    #Save the molecular structure as an image in memory
    image_stream = io.BytesIO()
    write(image_stream, atoms, format='png', scale = 100)
    image_stream.seek(0)

    #Convert the image to a base64-encoded string
    image_base64 = base64.b64encode(image_stream.read()).decode('utf-8')

    return image_base64