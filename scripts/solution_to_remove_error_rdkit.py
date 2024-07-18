"""
This script is a proposition to remove the error from the conversion of a
molecule to carbon atoms.

When one of the atoms in the molecule has a valence > 4, the conversion to
carbon atoms will fail (as the limit is 4). In this case, rdkit will print a
warning to stderr in the format:
[hh:mm:ss] Explicit valence for atom # x C, n, is greater than permitted
You can test the error with the SMILES "C1=C2C3=C[SH]134NN24"
"""

import ctypes
import io
import os
import sys
import tempfile
from contextlib import contextmanager
from typing import Optional

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

libc = ctypes.CDLL(None)
c_stderr = ctypes.c_void_p.in_dll(libc, "stderr")
c_flush = libc.fflush


@contextmanager
def stderr_redirector(stream):
    # The original fd stderr points to.
    original_stderr_fd = sys.stderr.fileno()

    def _redirect_stderr(to_fd):
        """Redirect stderr to the given file descriptor."""
        # Flush the C-level buffer stderr
        c_flush(c_stderr)
        # Flush and close sys.stderr - also closes the file descriptor (fd)
        sys.stderr.close()
        # Make original_stderr_fd point to the same file as to_fd
        os.dup2(to_fd, original_stderr_fd)
        # Create a new sys.stderr that points to the redirected fd
        sys.stderr = io.TextIOWrapper(os.fdopen(original_stderr_fd, "wb"))

    # Save a copy of the original stderr fd in saved_stderr_fd
    saved_stderr_fd = os.dup(original_stderr_fd)
    try:
        # Create a temporary file and redirect stderr to it
        tfile = tempfile.TemporaryFile(mode="w+b")
        _redirect_stderr(tfile.fileno())
        # Yield to caller, then redirect stderr back to the saved fd
        yield
        _redirect_stderr(saved_stderr_fd)
        # Copy contents of temporary file to the given stream
        tfile.flush()
        tfile.seek(0, io.SEEK_SET)
        stream.write(tfile.read())
    finally:
        tfile.close()
        os.close(saved_stderr_fd)


def try_convert_to_carbon(smiles: str) -> Optional[str]:
    """Try to convert the atoms in the molecule to carbon atoms.
    Should not work if the molecule contains P or S atoms with valence > 4.

    Args:
        smiles (str): SMILES representation of the molecule

    Returns:
        Optional[str]: SMILES representation of the molecule with all atoms
            converted to carbon atoms or None if it fails
    """
    f = io.BytesIO()
    with stderr_redirector(f):
        try:
            gscaf = MurckoScaffold.MakeScaffoldGeneric(
                MurckoScaffold.GetScaffoldForMol(Chem.MolFromSmiles(smiles))
            )
            smiles_scaffolds: str = Chem.MolToSmiles(gscaf)
            return smiles_scaffolds
        except Exception:
            pass

    output = f.getvalue().decode("utf-8")
    if output:
        print(f"Error with molecule {smiles}: {output}")
    return None


if __name__ == "__main__":

    try_convert_to_carbon("C1=C2C3=C[SH]134NN24")
