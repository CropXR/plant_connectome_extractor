import logging
from datetime import datetime
from pathlib import Path
from typing import List, Dict

import numpy as np
import pandas as pd
import requests
from Bio import Entrez as ez
from tqdm import tqdm


def parse_input_list_in_file(in_path: Path) -> List[str]:
    """Read file and split items in first line into list. Needed to read the
    input genes and phenotypes.

    :param in_path: path to input text file
    :return: list of strings
    """
    with in_path.open('r') as f:
        my_list = f.read().split()
    upper_case_list = [name.upper() for name in my_list]
    return upper_case_list


def keep_only_edges_with_relevant_nodes(df: pd.DataFrame,
                                        names_of_interest: List[str]):
    """Keep only edges in dataframe that go directly from one node of interest to another node of interest.
    Removes all other nodes and edges.

    :param df: Dataframe to be filtered, should be in the form of edgelist
    :param names_of_interest: List of gene/molecule/process/phenotype
                              names of interest
    :return: Filtered dataframe, that only contains edges that go to and from
             nodes that are in the names_of_interest
    """
    assert all(name.isupper() for name in names_of_interest), "Assert that gene names are uppercase"
    regex_to_match = '|'.join(names_of_interest)
    index = (df['Source'].str.contains(regex_to_match)
             & df['Target'].str.contains(regex_to_match))
    return df[index]


def extract_from_plantconnectome(out_dir: Path, roi_dict: Dict) -> pd.DataFrame:
    """Given a list of molecules and phenotypes, extract all their connections
    from plant connectome using the API.

    This is info that can be used to create a new network, or expand an
    existing network.

    :param out_dir: Out directory in which output should be saved
    :param roi_dict: Dict that contains name of gene/phenotype as key, and its
                     type (e.g. GM for gene/molecule) as value
    """
    out_dir.mkdir(exist_ok=True)
    # date = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    # sub_result_folder = out_dir / 'per_gene_result'
    # sub_result_folder.mkdir()

    totdf = pd.DataFrame()
    logging.info('Making API requests to plantconnectome...')
    # Save connectome output for it and annotate context
    for regul in tqdm(roi_dict.keys()):
        roitype = roi_dict.get(regul)
        res = requests.get(
            "https://connectome.plant.tools/normal/{0}".format(regul))
        try:
            df1 = pd.read_html(res.content)[0]
        except ValueError as e:
            logging.warning(f"\n{regul} normal search returned Error: '{e}'. "
                            f"Skipping that result.")
            df1 = pd.DataFrame()

        # If gene/molecule

        if roitype == "GM":
            search_type = 'alias'
        elif roitype == "PP":
            search_type = 'substring'
        else:
            raise KeyError(f'{roitype} is not a known interest type. '
                           f'Only GM and PP are supported.')

        res = requests.get(f"https://connectome.plant.tools/{search_type}/{regul}")
        try:
            df2 = pd.read_html(res.content)[0]
        except ValueError as e:
            logging.warning(f"\n{regul} {search_type} search returned Error: '{e}'. "
                            f"Skipping that result.\n")
            df2 = pd.DataFrame()

        # Merge two search type dfs
        df = pd.concat([df1, df2])
        df = df.drop_duplicates()

        df["query"] = regul
        df["querytype"] = roitype

        df.to_csv(out_dir / f"{regul}.tsv", sep="\t", index=False)
        totdf = pd.concat([totdf, df])

    totdf = totdf.drop_duplicates(subset=totdf.columns[:4])
    return totdf


def annotate_from_pmid(df: pd.DataFrame, keywords: list[str]) -> pd.DataFrame:
    """Add extra information to dataframe based on its pubmed-IDs

    :param df: dataframe that contains 'Pubmed ID' column
    :return: df with fields that describe the information of the study
    """
    # List PMIDs with Target or Source columns containing any of the keywords
    keyword_regex = '|'.join(keywords)
    pmids = np.concatenate((
        df[df.Target.str.contains(keyword_regex)]['Pubmed ID'].unique(),
        df[df.Source.str.contains(keyword_regex)]['Pubmed ID'].unique()
    ))
    # Annotate that information on df
    for index, row in df.iterrows():
        if row['Pubmed ID'] in pmids:
            df.loc[index, "HD_context"] = 1
        else:
            df.loc[index, "HD_context"] = 0
    df.loc[:, "HD_context"] = df["HD_context"].astype(int)
    # Fetch information from Pubmed for each PMID
    pmid = list(df['Pubmed ID'].unique())
    ez.email = "tijmen.vanbutselaar@wur.nl"
    handle = ez.efetch(db="pubmed", id=','.join(map(str, pmid)),
                       rettype="xml", retmode="text")
    records = ez.read(handle)
    abstracts = [pubmed_article['MedlineCitation']['Article']['Abstract'][
                     'AbstractText'][0]
                 if 'Abstract' in pubmed_article['MedlineCitation'][
        'Article'].keys()
                 else pubmed_article['MedlineCitation']['Article'][
        'ArticleTitle']
                 for pubmed_article in records['PubmedArticle']]
    rawoutput = []
    for i in records["PubmedArticle"]:
        try:
            rawoutput.append("{0}-{1}".format(
                str(i['MedlineCitation']['MeshHeadingList']),
                str(i['MedlineCitation']['KeywordList'])))
        except:
            rawoutput.append(str(i['MedlineCitation']['KeywordList']))
    abstract_dict = dict(zip(pmid, abstracts))
    raw_dict = dict(zip(pmid, rawoutput))
    df['abstract'] = df['Pubmed ID'].map(abstract_dict)
    df['raw'] = df['Pubmed ID'].map(raw_dict)
    for index, row in df.iterrows():
        # Look up if keyword mentioned in abstract of PMID
        abstract = row['abstract']
        if any(keywd in abstract.upper() for keywd in keywords):
            df.loc[index, "keywd"] = 1
        else:
            df.loc[index, "keywd"] = 0

    # Look up if Arabidopsis mentioned in abstract or keywords

        if ("Arabidopsis" in row['raw']) or (
                "Arabidopsis" in row['abstract']):
            df.loc[index, "Ath"] = 1
        else:
            df.loc[index, "Ath"] = 0

    df['keywd'] = df['keywd'].astype(int)
    df['Ath'] = df['Ath'].astype(int)
    df = df.drop(columns="raw")
    # Look up if root/shoot terminology used in abstract
    for index, row in df.iterrows():
        abstract = row['abstract']
        if "ROOT" in abstract.upper():
            df.loc[index, "root"] = 1
        else:
            df.loc[index, "root"] = 0
        if any(keywd in abstract.upper() for keywd in
               ["SHOOT", "LEAF", "HYPOCOTYL"]):
            df.loc[index, "shoot"] = 1
        else:
            df.loc[index, "shoot"] = 0
    df['root'] = df['root'].astype(int)
    df['shoot'] = df['shoot'].astype(int)
    return df
