from rdkit.Chem import Draw
from IPython.display import SVG, display
import matplotlib.pyplot as plt
from io import BytesIO
from PIL import Image
import cairosvg

import sys
import os


# Add the parent directory to the path to import the module evomol
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from evomol.representation import Molecule, SMILES, MolecularGraph

Molecule.id_representation_class = SMILES
Molecule.representations_class = [MolecularGraph]


def mol_to_svg(molecule: Molecule, size=500, draw_index=True):
    mol = molecule.get_representation(MolecularGraph).mol.GetMol()

    d = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)
    d.drawOptions().addStereoAnnotation = True
    if draw_index:
        d.drawOptions().addAtomIndices = True
    d.DrawMolecule(mol)
    d.FinishDrawing()

    svg = d.GetDrawingText()

    return svg


def draw_in_jupyter(svg: str) -> None:
    display(SVG(svg))


def draw_in_matplotlib(svg: str) -> None:
    # convert SVG to PNG using cairosvg
    png_data = cairosvg.svg2png(bytestring=svg.encode("utf-8"))

    # load the PNG image with PIL
    image = Image.open(BytesIO(png_data))

    # create a matplotlib figure and display the image
    fig, ax = plt.subplots(figsize=(image.width / 100, image.height / 100))
    ax.imshow(image)
    ax.axis("off")
    plt.show()


def draw_multiple_svgs_in_matplotlib(svgs: list) -> None:
    # convert each SVG to PNG using cairosvg
    png_datas = [
        Image.open(BytesIO(cairosvg.svg2png(bytestring=svg.encode("utf-8"))))
        for svg in svgs
    ]

    # set the number of columns and rows for the subplots
    num_images = len(png_datas)
    # 3 columns maximum
    cols = min(3, num_images)
    rows = (num_images + cols - 1) // cols

    # create a figure with subplots
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 5))

    # display each image in a subplot
    for i, (image, ax) in enumerate(zip(png_datas, axes.flat)):
        ax.imshow(image)
        ax.axis("off")

    # remove the axis for the empty subplots
    for j in range(i + 1, len(axes.flat)):
        axes.flat[j].axis("off")

    plt.tight_layout()
    plt.show()


molecules = [Molecule("CCO"), Molecule("CCC"), Molecule("CCN")]
svgs = [mol_to_svg(mol) for mol in molecules]
draw_multiple_svgs_in_matplotlib(svgs)
