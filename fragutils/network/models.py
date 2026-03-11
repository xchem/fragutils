from rdkit import Chem
from fragutils.utils.network_utils import get_num_ring_atoms, get_type, simplified_graph


class NodeHolder(object):
    """
    A Class to hold a series of nodes
    """

    def __init__(self, iso_flag=True):
        self.node_list = set()
        self.edge_list = set()
        self.iso_flag = iso_flag

    def create_or_retrieve_node(self, child_smi):
        new_node = Node(Chem.MolFromSmiles(child_smi), self.iso_flag)
        if new_node not in self.node_list:
            self.node_list.add(new_node)
            return new_node, True
        return new_node, False

    def create_or_retrieve_edge(self, excluded_smi, child_smi, input_node, new_node):
        """
        :param excluded_smi:
        :param child_smi:
        :param input_node:
        :param new_node:
        :return:
        """
        new_edge = Edge(excluded_smi, child_smi, input_node, new_node, self.iso_flag)
        new_edge.NODES = [input_node, new_node]
        self.edge_list.add(new_edge)
        return new_edge

    def get_edges(self):
        """
        :return:
        """
        return set(self.edge_list)


class Node(object):
    """
    A Class to hold a Node on the graph
    """

    def __eq__(self, other):
        return self.HASH == other.HASH

    def __ne__(self, other):
        return self.HASH != other.HASH

    def __hash__(self):
        return self.HASH

    def __init__(self, input_mol=None, iso_flag=True):
        self.iso_flag = iso_flag
        if not input_mol:
            return
        if type(input_mol) == str:
            input_mol = Chem.MolFromSmiles(input_mol)
        self.SMILES = Chem.MolToSmiles(input_mol, isomericSmiles=self.iso_flag)
        self.HAC = input_mol.GetNumHeavyAtoms()
        self.RAC, split_indices = get_num_ring_atoms(input_mol)
        self.RING_SMILES = simplified_graph(self.SMILES, iso_flag=self.iso_flag)
        self.RDMol = input_mol
        self.EDGES = []
        self.HASH = hash(self.SMILES)

    def __str__(self):
        return " ".join(
            ["NODE", self.SMILES, str(self.HAC), str(self.RAC), self.RING_SMILES]
        )


class Edge(object):
    """
    A Class to hold an edge
    """

    def __eq__(self, other):
        return self.HASH == other.HASH

    def __ne__(self, other):
        return self.HASH != other.HASH

    def __hash__(self):
        return self.HASH

    def __init__(self, excluded_smi, rebuilt_smi, node_one, node_two, iso_flag=True):
        self.EXCLUDE_SMILES = Chem.MolToSmiles(
            Chem.MolFromSmiles(excluded_smi), isomericSmiles=False
        )
        self.EXCLUDE_TYPE = get_type(excluded_smi)
        self.REBUILT_SMILES = Chem.MolToSmiles(
            Chem.MolFromSmiles(rebuilt_smi), isomericSmiles=False
        )
        self.REBUILT_RING_SMILES = simplified_graph(rebuilt_smi, iso_flag=iso_flag)
        self.REBUILT_TYPE = get_type(rebuilt_smi)
        self.EXCLUDED_RING_SMILES = simplified_graph(excluded_smi, iso_flag=iso_flag)
        self.NODES = [node_one, node_two]
        self.REPR = " ".join(
            ["EDGE", self.NODES[0].SMILES, self.NODES[1].SMILES, self.get_label()]
        )
        self.HASH = hash(
            "||".join(
                ["EDGE", self.NODES[0].SMILES, self.NODES[1].SMILES, self.get_label()]
            )
        )

    def get_label(self):
        return "|".join(
            [
                self.EXCLUDE_TYPE,
                self.EXCLUDE_SMILES,
                self.EXCLUDED_RING_SMILES,
                self.REBUILT_TYPE,
                self.REBUILT_SMILES,
                self.REBUILT_RING_SMILES,
            ]
        )

    def __str__(self):
        return self.REPR


class Attr(object):
    """
    A Class of attributes
    """

    def __init__(self, smiles=None, property_list=None, input_str=None):
        if input_str:
            self.SMILES = input_str.split()[1]
            self.PROP_LIST = input_str.split()[2:]
        else:
            self.SMILES = smiles
            self.PROP_LIST = property_list
        if self.SMILES is None:
            raise ValueError("Att")

    def __str__(self):
        tot_list = ["ATTR", self.SMILES]
        tot_list.extend(self.PROP_LIST)
        return " ".join(tot_list)
