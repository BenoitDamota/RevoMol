from rdkit.Chem import Draw
from IPython.display import SVG, display
import matplotlib.pyplot as plt
from io import BytesIO
from PIL import Image
import cairosvg

from evomol.representation import MolecularGraph, Molecule


def mol_to_svg(
    molecule: Molecule,
    size=500,
    draw_index=True,
    notes_on_atoms: dict[int, str] = None,
    notes_on_bonds: dict[int, str] = None,
) -> str:
    mol = molecule.get_representation(MolecularGraph).mol.GetMol()

    # other : Draw.rdMolDraw2D.MolDraw2DCairo
    d = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)
    d.drawOptions().addStereoAnnotation = True
    if draw_index:
        d.drawOptions().addAtomIndices = True

    if notes_on_atoms is not None:
        for atom_idx, note in notes_on_atoms.items():
            mol.GetAtomWithIdx(atom_idx).SetProp("atomNote", note)
    if notes_on_bonds is not None:
        for bond_idx, note in notes_on_bonds.items():
            mol.GetBondWithIdx(bond_idx).SetProp("bondNote", note)

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
