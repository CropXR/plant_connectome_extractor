"""
Scripts to infer, expand, extract interaction networks from plantconnectome.

Author: Tijmen van Butselaar, Ben Noordijk
"""

from pathlib import Path
from typing import Annotated

import typer
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

from helpers import parse_input_list_in_file, keep_only_edges_with_relevant_nodes, \
    extract_from_plantconnectome, annotate_from_pmid

app = typer.Typer()


@app.command()
def expand_network_of_interest(full_edgelist: Path,
                               original_edges: Path,
                               out_path: Path) -> nx.DiGraph:
    """Given an original network with genes of interest, and a network with
    all connections (between many genes), this function will add all nodes
    that are intermediate (i.e. situated between two nodes in
    the original network).

    :param full_edgelist: Edge list dataframe that contains all edges
    :param original_edges: Path to edge list (tsv) that only contains genes
                              of interest as nodes
    :return: Graph object of expanded network
    """
    full_edgelist_df = pd.read_csv(full_edgelist, sep='\t')
    original_edges_df = pd.read_csv(original_edges, sep='\t')

    full_edgelist_df['Interaction'] = "interacts with"
    original_edges_df['added_as_intermediate'] = False
    potential_expansion_network = nx.from_pandas_edgelist(full_edgelist_df,
                                              source='Source',
                                              target='Target',
                                              edge_attr=True,
                                              create_using=nx.DiGraph)
    nx.set_node_attributes(potential_expansion_network, True, 'extra')
    original_graph = nx.from_pandas_edgelist(original_edges_df,
                                             source='Source',
                                             target='Target',
                                             edge_attr=True,
                                             create_using=nx.DiGraph)
    nx.set_node_attributes(original_graph, False, 'extra')
    nx.draw_circular(original_graph, with_labels=True)
    plt.show()
    edges_to_add = []
    original_nodes = set(original_graph.nodes)
    for node in potential_expansion_network.nodes:
        if node in original_nodes:
            # No need to add nodes we already know
            continue
        predecessors = potential_expansion_network.predecessors(node)
        successors = potential_expansion_network.successors(node)

        overlapping_predecessors = original_nodes.intersection(predecessors)
        overlapping_successors = original_nodes.intersection(successors)
        # Remove overlap; ignore genes that both regulate
        # and are regulated by the same node
        no_dupl_predecessor = overlapping_predecessors - overlapping_successors
        no_dupl_successor = overlapping_successors - overlapping_predecessors
        if (len(no_dupl_predecessor) == 0
                or len(no_dupl_successor) == 0):
            # Not situated between two nodes from the original network
            continue
        for predecessor in no_dupl_predecessor:
            edge_data = potential_expansion_network.get_edge_data(predecessor, node)
            edge_data['added_as_intermediate'] = True
            edges_to_add.append((predecessor, node, edge_data))
            print(f'Adding {predecessor}->{node}')
        for successor in no_dupl_successor:
            edge_data = potential_expansion_network.get_edge_data(node, successor)
            edge_data['added_as_intermediate'] = True
            edges_to_add.append((node, successor, edge_data))
            print(f'Adding {node}->{successor}')

    original_graph.add_edges_from(edges_to_add)
    # Add node attribute that shows if it was in the original graph
    nx.draw_circular(original_graph, with_labels=True)
    plt.show()

    # TODO look at data storing, try to reduce redundancy
    out_df = nx.to_pandas_edgelist(original_graph)
    out_df.to_csv(out_path, sep='\t')

    # nx.write_edgelist(expanded_graph, 'expanded_testcase_PC.tsv', data=['Interaction', 'Interaction Type', 'PMIDs', 'Pubmed ID', 'added_as_intermediate'] delimiter='\t')

    return original_graph


@app.command()
def create_proto_network(
        out_dir: Annotated[
            Path,
            typer.Option(help='Directory in which output files will be placed')
        ],
        genes_of_interest: Annotated[
            Path,
            typer.Option(...,
                          help='Path to file that lists genes you are interested in,'
                                ' separated by space')
        ],
        pp_of_interest:  Annotated[
            Path,
            typer.Option(...,
                help='Path to file that lists proteins/phenotypes you are interested in,'
                     ' separated by space')
        ],
        re_extract_data: Annotated[
            bool,
            typer.Option("--re-extract-data",
                help="If provided, re-extract all data from plantconnectome")
        ] = False):
    """Given a list of genes/molecules and process/phenotypes, extract all their connections
    from plant connectome using the API.

    This is info that can be used to create a new network, or expand an
    existing network."""
    ## Dictionary that maps genes and process/phenotype of interest to a type
    # Ensure everything is uppercase
    # Gene/Molecule
    genes_of_interest_list = parse_input_list_in_file(genes_of_interest)
    roi_dict = {gene: 'GM' for gene in genes_of_interest_list}
    # Process/phenotype
    pp_of_interest_list = parse_input_list_in_file(pp_of_interest)
    roi_dict.update({pp: 'PP' for pp in pp_of_interest_list})

    if re_extract_data:
        all_edges = extract_from_plantconnectome(out_dir, roi_dict)
        all_edges.to_csv(out_dir / "totaldf.tsv", sep="\t", index=False)
    else:
        in_file = out_dir / 'totaldf.tsv'
        assert in_file.exists(), ("Have you extracted info from "
                                  "plantconnectome before? If not, "
                                  "do so with the --re-extract-data flag")
        all_edges = pd.read_csv(in_file, sep='\t')

    pruned_edges = keep_only_edges_with_relevant_nodes(
        all_edges, genes_of_interest_list + pp_of_interest_list)

    keywords = ["THERMO", "TEMPERATURE", "DROUGHT", "WATER", "OSMOTIC",
                "SALT"]
    pruned_edges = annotate_from_pmid(pruned_edges, keywords=keywords)
    pruned_edges.to_csv(out_dir / f"pruned_df.tsv", sep="\t", index=False)


if __name__ == '__main__':
    app()
