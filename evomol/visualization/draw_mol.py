"""
This module provides functions to draw molecules in Jupyter notebooks or in
matplotlib figures. The SVG format is used for the images because it provides
better quality than PNG images. The SVG images are converted to PNG using
cairosvg and then displayed in matplotlib figures.

You might want to look at the following link
http://www.rdkit.org/new_docs/Cookbook.html#drawing-molecules-jupyter
as it provides a way to draw molecules in Jupyter notebooks using RDKit.


"""

from io import BytesIO
from typing import Optional

import cairosvg
import matplotlib.pyplot as plt
from IPython.display import SVG, display
from PIL import Image
from rdkit.Chem import Draw

from evomol.representation import MolecularGraph, Molecule


def mol_to_svg(
    molecule: Molecule,
    size: int = 500,
    draw_index: bool = True,
    notes_on_atoms: Optional[dict[int, str]] = None,
    notes_on_bonds: Optional[dict[int, str]] = None,
) -> str:
    """Convert a molecule to an SVG image. You can add notes on atoms and bonds.

    Args:
        molecule (Molecule): The molecule to draw.
        size (int, optional): Size of the plot. Defaults to 500.
        draw_index (bool, optional): Draw the index of the atoms on the plot.
            Defaults to True.
        notes_on_atoms (dict[int, str], optional): notes to add on atoms.
            Format: {atom_index: note}. Defaults to None.
        notes_on_bonds (dict[int, str], optional): notes to add on bonds.
            Format: {bond_index: note}. Defaults to None.

    Returns:
        str: SVG image of the molecule.
    """
    mol = molecule.get_representation(MolecularGraph).mol.GetMol()

    # other way : Draw.rdMolDraw2D.MolDraw2DCairo
    d = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)

    d.drawOptions().addStereoAnnotation = True

    # add atom indices
    if draw_index:
        d.drawOptions().addAtomIndices = True

    # add notes on atoms and bonds
    if notes_on_atoms is not None:
        for atom_idx, note in notes_on_atoms.items():
            mol.GetAtomWithIdx(atom_idx).SetProp("atomNote", note)
    if notes_on_bonds is not None:
        for bond_idx, note in notes_on_bonds.items():
            mol.GetBondWithIdx(bond_idx).SetProp("bondNote", note)

    # draw the molecule
    d.DrawMolecule(mol)
    d.FinishDrawing()

    # get the SVG image and return it
    svg: str = d.GetDrawingText()
    return svg


def draw_in_jupyter(svg: str) -> None:
    """Display an SVG image in a Jupyter notebook. The quality is better than
    in matplotlib.

    Args:
        svg (str): SVG image to display.
    """
    display(SVG(svg))  # type: ignore[no-untyped-call]


def draw_in_matplotlib(svg: str, show: bool = True, save_to_path: str = "") -> None:
    """Display an SVG image in a matplotlib figure.

    Args:
        svg (str): SVG image to display.
        show (bool, optional): Show the image in the notebook or new window.
            Defaults to True.
        save_to_path (str, optional): Path to save the image.
            Defaults to "".
    """
    # convert SVG to PNG using cairosvg
    png_data = cairosvg.svg2png(bytestring=svg.encode("utf-8"))

    # load the PNG image with PIL
    image = Image.open(BytesIO(png_data))

    # create a matplotlib figure and display the image
    fig, ax = plt.subplots(figsize=(image.width / 100, image.height / 100))
    ax.imshow(image)
    ax.axis("off")

    if show:
        plt.show()

    if save_to_path:
        fig.savefig(save_to_path)


def draw_multiple_svgs_in_matplotlib(
    svgs: list[str],
    show: bool = True,
    save_to_path: str = "",
) -> None:
    """Display multiple SVG images in a matplotlib figure.

    Args:
        svgs (list): List of SVG images to display.
        show (bool, optional):  Show the image in the notebook or new window.
            Defaults to True.
        save_to_path (str, optional): Path to save the image.
            Defaults to "".
    """
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
    for image, ax in zip(png_datas, axes.flat):  # type: ignore[union-attr]
        ax.imshow(image)
        ax.axis("off")

    # remove the axis for the empty subplots
    for j in range(len(png_datas), len(axes.flat)):  # type: ignore[union-attr]
        axes.flat[j].axis("off")  # type: ignore[union-attr]

    plt.tight_layout()
    if show:
        plt.show()

    if save_to_path:
        fig.savefig(save_to_path)
