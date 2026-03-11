from rdkit import Chem
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.Chem import AllChem, RWMol, Atom
from fragutils.utils.network_utils import (
    rebuild_smi,
    get_ring_ring_splits,
    get_fragments,
    ret_comb_index,
    conv_at_xe,
)
import re

XE_PATT = r"[0-9]{3,}Xe"


def get_mol(mol_parser, mol_data):
    rd_mol = mol_parser(mol_data)
    return RWMol(AllChem.AddHs(rd_mol))


def decorate_smi(input_smi):
    """
    Decorate an input SMILES with a pseudeo molecule with all desirable changes.
    This can be added to the list of input molecules - and then we know if we get an Antimony back - it's not a real trans.
    :param input_smi: the input smiles
    :return: a list of SMILES around a given molecule.
    """
    # Add replacement groups (e.g. At) to all Ring positions in turn.
    return decorate_mol(get_mol, Chem.MolFromSmiles, input_smi)


def decorate_3d_mol(input_mol_file, get_indices=False):
    res_dict = decorate_mol(
        get_mol, Chem.MolFromMolBlock, input_mol_file, three_d_mol=True
    )
    out_dict = {}
    for res in res_dict:
        atom_pairs = res_dict[res]
        mol = AllChem.AddHs(Chem.MolFromMolBlock(input_mol_file), addCoords=True)
        conf = mol.GetConformer()
        for i, atom_pair in enumerate(atom_pairs):
            if get_indices:
                new_atom_pair = (atom_pair[0], atom_pair[1], "NA", True)
                out_dict[conv_at_xe(res) + "__" + str(i)] = [new_atom_pair]
            else:
                atom_one = conf.GetAtomPosition(atom_pair[0])
                atom_two = conf.GetAtomPosition(atom_pair[1])
                out_dict[conv_at_xe(res) + "__" + str(i)] = [
                    (atom_one.x, atom_one.y, atom_one.z),
                    (atom_two.x, atom_two.y, atom_two.z),
                ]
    return out_dict


def decorate_mol(get_new_mol, mol_parse, mol_data, three_d_mol=False):
    mol = get_new_mol(mol_parse, mol_data)
    patt = Chem.MolFromSmarts("[*;R]-;!@[H]")
    # Get the list of atom Indices to replace
    out_atom_repls = [x for x in mol.GetSubstructMatches(patt)]
    # Now replace with At - an produce a new mol everytime
    new_mols = {}
    for atom_pairs in out_atom_repls:
        atom = atom_pairs[1]
        rw_mol = get_new_mol(mol_parse, mol_data)
        rw_mol.ReplaceAtom(atom, Atom(85))
        newer_mol = rw_mol.GetMol()
        this_mol = Chem.MolToSmiles(
            Chem.MolFromSmiles(Chem.MolToSmiles(newer_mol, isomericSmiles=True)),
            isomericSmiles=True,
        )
        if three_d_mol:
            if this_mol in new_mols:
                new_mols[this_mol].append(atom_pairs)
            else:
                new_mols[this_mol] = [atom_pairs]
        else:
            new_mols[this_mol] = atom_pairs
    return new_mols


def deletion_linker_smi(input_smi, iso_labels=True):
    mol = Chem.MolFromSmiles(input_smi)
    return deletion_linker_mol(mol, iso_labels)


def deletion_linker_sd(input_mol, iso_labels=True):
    mol = Chem.MolFromMolBlock(input_mol)
    return deletion_linker_mol(mol, iso_labels)


def get_3d_vects_for_mol(input_mol, iso_labels=True):
    tot_dict = del_link_coord(input_mol, iso_labels=iso_labels)
    tot_dict["additions"] = decorate_3d_mol(input_mol)
    return tot_dict


def get_vect_indices_for_mol(input_mol):
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(Chem.MolFromMolBlock(input_mol)))
    AllChem.EmbedMolecule(mol)
    input_mol = Chem.MolToMolBlock(mol)
    tot_dict = del_link_coord(input_mol, get_indices=True, iso_labels=False)
    tot_dict["additions"] = decorate_3d_mol(input_mol, get_indices=True)
    return tot_dict


def convert_dict_to_mols(tot_dict):
    """
    :param tot_dict:
    :return:
    """
    mol_list = []
    for smiles in tot_dict:
        # Now generate the molecules for that
        mol = RWMol()
        atoms = tot_dict[smiles]
        print(atoms)
        for atom in atoms:
            atom = Atom(6)
            mol.AddAtom(atom)
        # for i in range(len(atoms)-1):
        #     mol.AddBond(i,i+1)
        mol = mol.GetMol()
        AllChem.EmbedMolecule(mol)
        conf = mol.GetConformer()
        for i, atom in enumerate(atoms):
            point_3d = Point3D(atom[0], atom[1], atom[2])
            conf.SetAtomPosition(i, point_3d)
        mol = conf.GetOwningMol()
        mol.SetProp("_Name", smiles)
        mol_list.append(mol)
    return mol_list


def del_link_coord(input_mol, get_indices=False, iso_labels=False):
    tot_vals = deletion_linker_sd(input_mol, iso_labels=iso_labels)
    deletions = [Chem.MolToSmiles(x, isomericSmiles=True) for x in tot_vals[0]]
    linkers = [Chem.MolToSmiles(x, isomericSmiles=True) for x in tot_vals[1]]
    ring = [Chem.MolToSmiles(x, isomericSmiles=True) for x in tot_vals[2]]

    if not iso_labels:
        iso_vals = deletion_linker_sd(input_mol, iso_labels=True)
    else:
        iso_vals = tot_vals
    deletions_iso = [Chem.MolToSmiles(x, isomericSmiles=True) for x in iso_vals[0]]
    linkers_iso = [Chem.MolToSmiles(x, isomericSmiles=True) for x in iso_vals[1]]
    ring_iso = [Chem.MolToSmiles(x, isomericSmiles=True) for x in iso_vals[2]]

    out_d = {"linkers": {}, "deletions": {}, "ring": {}}
    for i, x in enumerate(linkers):
        ret_ans = get_atom_coords(
            x, Chem.MolFromMolBlock(input_mol), get_indices, linkers_iso[i]
        )
        out_d["linkers"][ret_ans[0]] = ret_ans[1:]
    for i, x in enumerate(deletions):
        ret_ans = get_atom_coords(
            x, Chem.MolFromMolBlock(input_mol), get_indices, deletions_iso[i]
        )
        out_d["deletions"][ret_ans[0]] = ret_ans[1:]
    for i, x in enumerate(ring):
        ret_ans = get_atom_coords(
            x, Chem.MolFromMolBlock(input_mol), get_indices, ring_iso[i]
        )
        out_d["ring"][ret_ans[0]] = ret_ans[1:]
    return out_d


def deletion_linker_mol(mol, iso_labels=True):
    """
    Produce all the linker deletion SMILES and replace with Li
    :param input_smi: the input SMI
    :return:
    """
    nr = mol.GetRingInfo().NumRings()
    fragments = get_fragments(mol, iso_labels=iso_labels)
    out_mols = []
    linker_mol_list = []
    ring_repl_list = []
    ring_ring_splits = get_ring_ring_splits(mol, do_comb_index=True)
    if ring_ring_splits:
        for ring_ring_split in ring_ring_splits:
            rebuilt_smi = rebuild_smi(ring_ring_split, ring_ring=True)
            new_mol = Chem.MolFromSmiles(rebuilt_smi)
            if new_mol.GetRingInfo().NumRings() < nr:
                continue
            linker_mol_list.append(new_mol)
    for i in range(len(fragments)):
        new_list = []
        for j, item in enumerate(fragments):
            if i == j:
                continue
            new_list.append(item)
        rebuilt_smi = rebuild_smi(new_list, ring_ring=False)
        if "." in rebuilt_smi:
            new_mol = Chem.MolFromSmiles(rebuilt_smi)
            # Only consider linker replacements
            if new_mol.GetRingInfo().NumRings() < nr:
                ring_repl_list.append(new_mol)
            else:
                linker_mol_list.append(new_mol)
            continue
        new_mol = Chem.MolFromSmiles(rebuilt_smi)
        # If the resulting
        if new_mol.GetRingInfo().NumRings() < nr:
            continue
        out_mols.append(new_mol)
    return out_mols, linker_mol_list, ring_repl_list


def find_atom_pairs(smiles_input, get_indices, iso_smiles):
    """
    Find the indices of the atom pairs from a SMILES input
    :param smiles_input:
    :return:
    """
    matches = re.findall(XE_PATT, smiles_input)
    if iso_smiles:
        iso_matches = re.findall(XE_PATT, iso_smiles)
    else:
        iso_matches = None
        isotope = None
    ind_list = []
    for i, match in enumerate(matches):
        index = int(match[:-2])
        if iso_matches:
            isotope = int(iso_matches[i][:-2])
        indices = ret_comb_index(index, get_indices, isotope)
        ind_list.append(indices)
    return ind_list


def get_atom_coords(smiles_input, mol, get_indices=False, iso_smiles=None):
    """
    Get a
    :param smiles_input:
    :param mol:
    :return:
    """
    ind_list = find_atom_pairs(smiles_input, get_indices, iso_smiles)
    clean_smi = re.sub(XE_PATT, "Xe", smiles_input)
    out_list = [clean_smi]
    conf = mol.GetConformer()
    for atom_pair in ind_list:
        atom_one = conf.GetAtomPosition(atom_pair[0])
        atom_two = conf.GetAtomPosition(atom_pair[1])
        if get_indices:
            out_list.append((atom_pair[0], atom_pair[1], atom_pair[2], False))
        else:
            out_list.extend(
                [
                    (atom_one.x, atom_one.y, atom_one.z),
                    (atom_two.x, atom_two.y, atom_two.z),
                ]
            )
    return out_list


def link_li(rebuilt_smi):
    mol = Chem.MolFromSmiles(rebuilt_smi)
    mol = RWMol(mol)
    bons = [x[0] for x in mol.GetSubstructMatches(Chem.MolFromSmarts("[Xe]"))]
    mol.AddBond(bons[0], bons[1])
    return mol.GetMol()


def addition_smi(input_smi):
    smis = decorate_smi(input_smi)
    return [Chem.MolFromSmiles(x) for x in smis]


def get_ring_removals(smi):
    rw_mol = RWMol(Chem.MolFromSmiles(smi))
    rings = rw_mol.GetRingInfo().AtomRings()
    out_mols = {}
    for ring in rings:
        new_mol = Chem.MolFromSmiles(smi)
        for atom in ring:
            new_mol.GetAtomWithIdx(atom).SetAtomicNum(0)
        Chem.DeleteSubstructs(new_mol, Chem.MolFromSmarts("[#0]"))
        Chem.GetMolFrags(new_mol)
        out_mols[Chem.MolToSmiles(new_mol, isomericSmiles=True)] = ring
    return out_mols


def conv_smiles(additions, deletions, linkers, ring_removals):
    additions = [conv_at_xe(Chem.MolToSmiles(x)) for x in additions]
    deletions = [Chem.MolToSmiles(x) for x in deletions]
    linkers = [Chem.MolToSmiles(x) for x in linkers]
    ring_removals = [Chem.MolToSmiles(x) for x in ring_removals]
    linkers.extend(ring_removals)
    return additions, deletions, linkers, ring_removals


def get_add_del_link(smi, asSmiles=True):
    additions = addition_smi(smi)
    res = deletion_linker_smi(smi)
    linkers = res[1]
    ring_removals = res[2]
    deletions = res[0]
    if asSmiles:
        additions, deletions, linkers, ring_removals = conv_smiles(
            additions, deletions, linkers, ring_removals
        )
    return [additions, deletions, linkers]
