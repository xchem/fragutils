import random
import logging

from fragutils.utils.network_utils import write_results, get_driver, canon_input

logger = logging.getLogger(__name__)


class ReturnObject(object):
    def __init__(
        self,
        start_smi,
        end_smi,
        label,
        edge_count,
        change_frag,
        iso_label,
        compound_ids,
    ):
        """
        Build this object.
        :param start_smi:
        :param end_smi:
        :param label:
        :param iso_label:
        :param frag_type:
        :param edge_count:
        :param change_frag:
        :param compound_ids: - list of compound_ids for the end_smi.
        """
        self.start_smi = start_smi
        self.end_smi = end_smi
        self.label = label
        self.iso_label = iso_label
        self.frag_type = None
        self.edge_count = edge_count
        self.change_frag = change_frag
        self.compound_ids = compound_ids

    def __str__(self):
        out_list = [self.label, str(self.edge_count), self.frag_type]
        return "_".join(out_list)


def find_double_edge_query(smiles):
    return (
        "MATCH (sta:F2 {smiles:'%(smiles)s'})-[nm:FRAG]-(mid:F2)-[ne:FRAG]-(end:Mol)"
        " where abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1"
        " and sta.smiles <> end.smiles"
        " RETURN sta, nm, mid, ne, end"
        " order by split(nm.label, '|')[4],"
        " split(ne.label, '|')[2];" % {"smiles": smiles}
    )


def find_triple_edge_growth_query(
    smiles,
    heavy_atom_diff_min=6,
    heavy_atom_diff_max=10,
    mid_heavy_atom_diff_min=-1,
    mid_heavy_atom_diff_max=3,
):
    return (
        "MATCH (sta:F2 {smiles:'%(smiles)s'})-[nm:FRAG]-(mid_one:F2)-[ne:FRAG]-(mid:Mol)-[nm2:FRAG]-(mid_two:F2)-[ne2:FRAG]-(end:Mol)"
        " where end.hac-sta.hac > %(hacmin)d and end.hac-sta.hac <= %(hacmax)d"
        " and mid.hac-sta.hac > %(chacmin)d and mid.hac-sta.hac <= %(chacmax)d"
        " and sta.smiles <> mid.smiles and sta.smiles <> end.smiles "
        " WITH collect("
        "{"
        "end: end.smiles,"
        "mid: mid.smiles,"
        "frag_one: mid_one.smiles,"
        "frag_two: mid_two.smiles"
        "}"
        ") AS edges"
        " RETURN edges"
        % {
            "smiles": smiles,
            "hacmin": heavy_atom_diff_min,
            "hacmax": heavy_atom_diff_max,
            "chacmin": mid_heavy_atom_diff_min,
            "chacmax": mid_heavy_atom_diff_max,
        }
    )


def add_follow_ups_query(smiles):
    return (
        "MATCH (sta:F2 {smiles:'%(smiles)s'})-[nm:FRAG]-(mid:F2)-[ne:FRAG]-(end:Mol)"
        " where abs(sta.hac-end.hac) <= 3 and abs(sta.chac-end.chac) <= 1"
        " and sta.smiles <> end.smiles "
        " MERGE (end)-[:FOLLOW_UP]->(sta)" % {"smiles": smiles}
    )


def find_proximal_query(smiles):
    return (
        "match p = (n:F2{smiles:'%(smiles)s'})-[nm]-(m:Mol)"
        " where abs(n.hac-m.hac) <= 3 and abs(n.chac-m.chac) <= 1"
        " return n, nm, m"
        " order by split(nm.label, '|')[4];" % {"smiles": smiles}
    )


def get_type(r_group_form, sub_one, sub_two):
    if "." in r_group_form:
        if "C1" in sub_two:
            return "ring_linker"
        return "linker"
    if "C1" in sub_two:
        return "ring_replacement"
    return "replacement"


def define_double_edge_type(record):
    """
    Create return object from the graph query results returned for double edge jump molecules
    :param: query record
    :return: ReturnObject
    """
    mol_one = record["sta"]
    label = str(record["ne"]["label"].split("|")[4])
    iso_label = str(record["ne"]["label"].split("|")[5])
    change_frag = str(record["ne"]["label"].split("|")[2])
    mol_two = record["mid"]
    mol_three = record["end"]
    diff_one = mol_one["hac"] - mol_two["hac"]
    diff_two = mol_two["hac"] - mol_three["hac"]
    ret_obj = ReturnObject(
        mol_one["smiles"],
        mol_three["smiles"],
        label,
        2,
        change_frag,
        iso_label,
        mol_three["cmpd_ids"],
    )
    if "." in label:
        ret_obj.frag_type = "LINKER"
    elif diff_one >= 0 and diff_two >= 0:
        ret_obj.frag_type = "DELETION"
    elif diff_one <= 0 and diff_two <= 0:
        ret_obj.frag_type = "ADDITION"
    else:
        ret_obj.frag_type = "REPLACE"
    return ret_obj


def define_proximal_type(record):
    """
    Create return object from the graph query results returned for (naerby) proximal molecules
    :param: query record
    :return: ReturnObject
    """
    mol_one = record["n"]
    label = str(record["nm"]["label"].split("|")[4])
    iso_label = str(record["nm"]["label"].split("|")[5])
    change_frag = str(record["nm"]["label"].split("|")[2])
    mol_two = record["m"]
    ret_obj = ReturnObject(
        mol_one["smiles"],
        mol_two["smiles"],
        label,
        1,
        change_frag,
        iso_label,
        mol_two["cmpd_ids"],
    )
    if "." in label:
        ret_obj.frag_type = "LINKER"
    elif mol_one["hac"] - mol_two["hac"] > 0:
        ret_obj.frag_type = "DELETION"
    elif mol_one["hac"] - mol_two["hac"] < 0:
        ret_obj.frag_type = "ADDITION"
    else:
        ret_obj.frag_type = "REPLACE"
    return ret_obj


def organise(records, num_picks):
    """
    Create Output dictionary from the list of ReturnObject records
    :param: ReturnObject records
    :return: Output dictionary
    """
    out_d = {}
    smi_set = set()
    for rec in records:
        rec_key = str(rec)
        addition = {
            "change": rec.change_frag,
            "end": rec.end_smi,
            "compound_ids": rec.compound_ids,
        }
        if rec_key in out_d:
            out_d[rec_key]["addition"].append(addition)
        else:
            out_d[rec_key] = {"vector": rec.iso_label, "addition": [addition]}
        smi_set.add(rec.end_smi)
    if num_picks:
        max_per_hypothesis = num_picks / len(out_d)
    out_smi = []
    for rec in out_d:
        # TODO here is the logic as to ordering replacements
        if num_picks:
            random.shuffle(out_d[rec]["addition"])
            out_d[rec]["addition"] = out_d[rec]["addition"][:max_per_hypothesis]
        else:
            out_d[rec]["addition"] = out_d[rec]["addition"]
        out_smi.extend(out_d[rec])
    return out_d


def get_picks(smiles, num_picks, graph_url="neo4j", graph_auth="neo4j/neo4j"):
    smiles = canon_input(smiles)
    driver = get_driver(graph_url, graph_auth)
    with driver.session() as session:
        records = []
        for record in session.run(find_proximal_query(smiles)):
            ans = define_proximal_type(record)
            records.append(ans)
        for record in session.run(find_double_edge_query(smiles)):
            ans = define_double_edge_type(record)
            records.append(ans)
        for label in list(set([x.label for x in records])):
            # Linkers are meaningless
            if "." in label:
                continue
        if records:
            orga_dict = organise(records, num_picks)
            return orga_dict
        else:
            print("Nothing found for input: " + smiles)


def get_full_graph(
    smiles, graph_url="neo4j", graph_auth="neo4j/neo4j", isomericSmiles=True
):
    """
    For the given smiles and graph credentials, query the Neo4j database
    :param: smiles
    :param: graph_url - note that for local Neo4j graphs, this can be set to 'localhost'
    :param: graph_auth - username/password.
    :param: isomericSmiles True/False
    :return: Output dictionary - This is returned "as-is" from the Fragalysis-Backend api/graph query.
    """

    msg = f'get_full_graph("{smiles}", graph_url={graph_url}, isomericSmiles={isomericSmiles})'
    smiles = canon_input(smiles, isomericSmiles)
    msg += f' canon_smiles="{smiles}"'

    driver = get_driver(graph_url, graph_auth)
    records = []
    with driver.session() as session:
        debug_proximal_num_records = 0
        for record in session.run(find_proximal_query(smiles)):
            debug_proximal_num_records += 1
            ans = define_proximal_type(record)
            records.append(ans)
        msg += f" n_proximals={debug_proximal_num_records}"
        debug_double_edge_num_records = 0
        for record in session.run(find_double_edge_query(smiles)):
            debug_double_edge_num_records += 0
            ans = define_double_edge_type(record)
            records.append(ans)
        msg += f" n_double_edge={debug_double_edge_num_records}"
        for label in list(set([x.label for x in records])):
            # Linkers are meaningless
            if "." in label:
                continue

    if records:
        orga_dict = organise(records, None)
        orga_dict_len = len(orga_dict)

        msg += f" orga_dict_len={orga_dict_len}"
        logger.info(f"GRAPH -> {msg}")
        return orga_dict

    msg += f" (fall-through) ret=None"
    logger.info(f"GRAPH -> {msg}")
    return None


def custom_query(query, graph_url="neo4j", graph_auth="neo4j/neo4j"):
    driver = get_driver(graph_url, graph_auth)
    records = []
    with driver.session() as session:
        for record in session.run(query):
            records.append(record)
    return records


def write_picks(smiles, num_picks):
    img_dict = write_results(get_picks(smiles, num_picks))
    for key in img_dict:
        out_f = open(key + ".svg", "w")
        out_f.write(img_dict[key])
