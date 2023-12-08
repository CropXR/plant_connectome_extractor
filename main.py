"""
Scripts to infer, expand, extract interaction networks from plantconnectome.

Author: Tijmen van Butselaar, Ben Noordijk
"""
import logging
from pathlib import Path
from typing import Annotated

import typer
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

from helpers import (parse_input_list_in_file,
                     keep_only_edges_with_relevant_nodes,
                     extract_from_plantconnectome,
                     annotate_from_pmid)

app = typer.Typer()

def do_one_iter_network_expansion(
        input_edgelist_path: Annotated[
            Path,
            typer.Option(help='Path to edge list tsv of network '
                              'that needs to be expanded')
        ],
        connectome_result_folder: Annotated[
            Path,
            typer.Option(help='Path to folder in which (existing) '
                              'plantconnectome requests are located')
        ],
        keywords_of_interest: Annotated[
            Path,
            typer.Option(
                help='Path to file that lists keywords you are interested in,'
                     ' separated by space. Used for annotation of edges')
        ],
        out_path: Annotated[
            Path,
            typer.Option(help='Path to output directory')
        ],):

    original_edges_df = pd.read_csv(input_edgelist_path, sep='\t')
    original_edges_df['added_as_intermediate'] = False
    original_graph = nx.from_pandas_edgelist(original_edges_df,
                                             source='Source',
                                             target='Target',
                                             edge_attr=True,
                                             create_using=nx.DiGraph)

    connectome_result_folder.mkdir(exist_ok=True)

    # Request plantconnectome data that is not yet present
    nodes_with_no_info = []
    for node in original_graph:
        if len(list(connectome_result_folder.glob(f'{node}.tsv'))) == 0:
            nodes_with_no_info.append(node)

    extract_from_plantconnectome(connectome_result_folder,
                                 {node: "GM"
                                  for node in nodes_with_no_info}
                                 )
    connection_df_list = []

    for node in original_graph:
        try:
            one_df = pd.read_csv(connectome_result_folder / f'{node}.tsv', sep='\t')
            connection_df_list.append(one_df)
        except FileNotFoundError:
            logging.warning(f"{node}.tsv wasn't found")

    new_connections = pd.concat(connection_df_list)

    potential_expansion_network = nx.from_pandas_edgelist(
        new_connections,
        source='Source',
        target='Target',
        edge_attr=True,
        create_using=nx.DiGraph
    )
    nx.set_node_attributes(potential_expansion_network, True, 'extra')

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
    edge_color_map = ['red' if has_been_added else 'black'
                      for source, target, has_been_added
                      in original_graph.edges(data='added_as_intermediate')]
    nx.draw_circular(original_graph, with_labels=True, edge_color=edge_color_map)
    plt.show()

    out_df = nx.to_pandas_edgelist(original_graph)
    out_df.rename(columns={'source': 'Source', 'target': 'Target'}, inplace=True)

    keywords = parse_input_list_in_file(keywords_of_interest)
    out_df = annotate_from_pmid(out_df, keywords)
    out_df.to_csv(out_path, sep='\t', index=False)

    return original_graph


@app.command()
def expand_network_of_interest(
        input_edgelist_path: Annotated[
            Path,
            typer.Option('--input-edgelist', '-i',help='Path to edge list tsv of network '
                              'that needs to be expanded')
        ],
        connectome_result_folder: Annotated[
            Path,
            typer.Option('--connectome-result-folder', '-f',help='Path to folder in which (existing) '
                              'plantconnectome requests are located')
        ],
        keywords_of_interest: Annotated[
            Path,
            typer.Option('--keywords', '-k',
                help='Path to file that lists keywords you are interested in,'
                     ' separated by space. Used for annotation of edges')
        ],
        out_path: Annotated[
            Path,
            typer.Option('--out-path', '-o',
                         help='Path to output directory. Files will be'
                              ' created in this directory automatically.')
        ],
        nr_iters: Annotated[
            int,
            typer.Option('--nr-iters', '-n',
                         help='Number of iterations, i.e. how many steps '
                              'between nodes will be considered for network'
                              ' inference')
        ]):
    """Given an original network with genes of interest, this function will
    add all nodes that are intermediate (i.e. situated between two nodes in
    the original network).
    """
    out_file_basis = input_edgelist_path.stem
    if not out_file_basis.endswith('_expand_iter'):
        out_file_basis = out_file_basis + '_expand_iter'
    input_file = input_edgelist_path
    out_path.mkdir(exist_ok=True)
    for i in range(1, nr_iters + 1):
        output_file = out_path / f'{out_file_basis}{i}.tsv'
        do_one_iter_network_expansion(input_file,
                                      connectome_result_folder,
                                      keywords_of_interest,
                                      output_file)
        input_file = output_file


@app.command()
def create_proto_network(
        out_dir: Annotated[
            Path,
            typer.Option('--out-dir', '-o',
                         help='Directory in which output files will be placed')
        ],
        genes_of_interest: Annotated[
            Path,
            typer.Option('--genes', '-g',
                         help='Path to file that lists genes you are interested in,'
                                ' separated by space')
        ],
        pp_of_interest:  Annotated[
            Path,
            typer.Option('--protein-phenotype', '-p',
                         help='Path to file that lists proteins/phenotypes you are interested in,'
                     ' separated by space')
        ],
        keywords_of_interest: Annotated[
            Path,
            typer.Option('--keywords', '-k',
                         help='Path to file that lists keywords you are interested in,'
                              ' separated by space. Used to annotate edges in cytoscape ')
        ] = None,
        re_extract_data: Annotated[
            bool,
            typer.Option("--re-extract-data", '-r',
                help="If provided, re-extract all data from plantconnectome")
        ] = False):
    """Given a list of genes/molecules and process/phenotypes, extract all their connections
    from plant connectome using the API.

    This is info that can be used to create a new network, or expand an
    existing network."""
    # Gene/Molecule
    genes_of_interest_list = parse_input_list_in_file(genes_of_interest)
    # Dictionary that maps genes and process/phenotype of interest to a type
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
    if keywords_of_interest:
        keywords = parse_input_list_in_file(keywords_of_interest)
        pruned_edges = annotate_from_pmid(pruned_edges, keywords=keywords)
    pruned_edges.to_csv(out_dir / f"pruned_df.tsv", sep="\t", index=False)


if __name__ == '__main__':
    app()
