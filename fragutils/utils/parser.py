# Series of functions to parse input files
from fragutils.alysis.models import Object, Owner
from fragutils.utils.rdkit_utils import (
    _parse_ligand_sdf,
    _get_c_of_mass,
    RDKitPh4,
    RDKitAtom,
    _get_water_coords,
    _get_waters,
    _get_res,
    _get_res_rmsds,
    _parse_pdb,
)
from rdkit import Chem


def _get_c_of_mass_list(mols):
    c_of_mass_list = []
    for m in mols:
        c_of_mass_list.append(_get_c_of_mass(m))
    return c_of_mass_list


def parse_ligands(input_file, input_type="sdf"):
    mols = _parse_ligand_sdf(input_file=input_file)
    # Now return them with their name and centre of mass
    return _get_c_of_mass_list(mols)


def parse_ligand_ph4s(input_mols):
    """
    Function to return a series of ligand based pharmacophores.
    :param input_mols: the RDKit molecules
    :return: the molecule based pharmacophores
    """
    rdkit_ph4 = RDKitPh4()
    rdkit_atom = RDKitAtom()
    output_pharma_list = []
    for mol in input_mols:
        if not mol:
            pharma_list = []
        else:
            pharma_list = rdkit_ph4.generate_ph4_for_mol(rdmol=mol)
            atom_list = rdkit_atom.generate_atoms_for_mol(mol)
            x, y, z = _get_c_of_mass(rdmol=mol)
            c_of_m_feat = (x, y, z, "c_of_m")
            pharma_list.append(c_of_m_feat)
            pharma_list.extend(atom_list)
        output_pharma_list.append(pharma_list)
    return output_pharma_list


def parse_waters(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of waters - return waters.
    :param input_pdb: the input PDB files
    :param input_mol: the input molecule (to use as a reference around which to limit)
    :return: tuple threes of coordinates of the waters
    """
    output_water_list = []
    # First just get the waters from the file
    for input_pdb in input_pdbs:
        waters = _get_waters(open(input_pdb).readlines())
        output_water_list.append(_get_water_coords(waters))
    return output_water_list


def parse_residues(input_pdbs, input_mol=None, max_dist=10.0):
    """
    Function to parse a series of PDB files of proteins.
    :param input_pdb: the input PDB files - with identifiers
    :return: a dict (Key Residue -> value list of molecules)
    """
    owner_list = []
    res_dict = {}
    for input_pdb in input_pdbs:
        # Loop through the residues
        mol = _parse_pdb(input_pdb)
        this_res_dict = _get_res(mol)
        for key in this_res_dict:
            if key in res_dict:
                res_dict[key].append(res_dict[key])
            else:
                res_dict[key] = [res_dict[key]]
    for res in res_dict:
        rmsd_coords = _get_res_rmsds(res_dict[res])
        out_l = []
        res = Object(rmsd_coords, res)
        out_l.append(res)
        owner = Owner(out_l, input_pdb)
        owner_list.append(owner)
    return owner_list


def get_file(file_path, output_format, file_counter):
    if output_format == "smi":
        return Chem.SmilesWriter(file_path + "_" + str(file_counter) + ".smi")
    else:
        return Chem.SDWriter(file_path + "_" + str(file_counter) + ".sdf")


def parse_mols(input_file, input_format):
    if input_format == "smi":
        return Chem.SmilesMolSupplier(input_file)
    else:
        return Chem.SDMolSupplier(input_file)
