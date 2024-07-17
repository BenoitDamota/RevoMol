import io
import os

from PIL import Image
from rdkit.Chem import Draw
from rdkit.Chem.rdDepictor import Compute2DCoords

from evomol.representation import MolecularGraph, Molecule


def draw(
    molecule: Molecule,
    show: bool = True,
    size: int = 200,
    write_to_path: str = "",
) -> None:
    """
    Drawing the molecule
    """

    mol = molecule.get_representation(MolecularGraph).mol.GetMol()

    # Computing coordinates and making sure the properties are computed
    Compute2DCoords(mol)
    mol.UpdatePropertyCache()

    # Drawing the molecule
    dr = Draw.rdMolDraw2D.MolDraw2DCairo(size, size)
    opts = dr.drawOptions()

    # Transparent background if not writing to file
    if not write_to_path:
        opts.clearBackground = False

    dr.DrawMolecule(mol)
    dr.FinishDrawing()

    # Loading the molecule as a PIL object
    bytes_images = dr.GetDrawingText()
    image = Image.open(io.BytesIO(bytes_images))

    if show:
        image.show()

    if write_to_path:
        # Creating directories if they don't exist
        os.makedirs(os.path.dirname(write_to_path), exist_ok=True)

        # Writing image to disk
        image.save(write_to_path, "PNG")

    # return image


def draw_with_index(molecule: Molecule) -> None:

    mol = molecule.get_representation(MolecularGraph).mol.GetMol()

    d = Draw.rdMolDraw2D.MolDraw2DCairo(1000, 1000)  # or MolDraw2DSVG to get SVGs
    d.drawOptions().addStereoAnnotation = True
    d.drawOptions().addAtomIndices = True
    d.DrawMolecule(mol)
    d.FinishDrawing()
    d.WriteDrawingText(f"{molecule}.png")

    # mol.GetAtomWithIdx(2).SetProp("atomNote", "foo")
    # mol.GetBondWithIdx(0).SetProp("bondNote", "bar")
