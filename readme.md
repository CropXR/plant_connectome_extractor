# What this does
Generate initial networks from groups of genes, or expand already known networks. Or hell; even combine the two of them

# How to install
```
git clone https://github.com/bnoordijk/plant_connectome_extractor.git
```

```
mamba env create --file environment.yml 
```

# How to run
```
python main.py --help
```
```
╭─ Commands ───────────────────────────────────────────────────────────────────╮
│ create-proto-network        Given a list of genes/molecules and/or           │
│                             process/phenotypes, extract all their            │
│                             connections from plant connectome using the API. │
│ expand-network-of-interest  Given an original network with genes of          │
│                             interest, this function will add all nodes that  │
│                             are intermediate (i.e. situated between two      │
│                             nodes in the original network).                  │
╰──────────────────────────────────────────────────────────────────────────────╯
```

## Generate initial network from list of genes
### Example:

```
python main.py create-proto-network create-proto-network
--out-dir
some_directory
--genes-of-interest
input_files/subset_genes.txt
--pp-of-interest
input_files/subset_phenotypes.txt
--keywords-of-interest
input_files/keywords.txt
```


### For more info:
```
python main.py create-proto-network --help
```

## Expand network of genes
### Example:
```
python main.py expand-network-of-interest
--input-edgelist-path
test_dir/small_df.tsv
--connectome-result-folder
test_dir/per_gene_result
--keywords-of-interest
input_files/keywords.txt
--out-path
test_dir/
--nr-iters
3
```

### For more info:
```
python main.py expand-network-of-interest --help
```