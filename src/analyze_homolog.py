import requests
import os

# RCSB PDB Search API URL
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v1/query?json="

# Download URL for PDB files
DOWNLOAD_URL = "https://files.rcsb.org/download"

# Define your search criteria: organism and protein name
organism_name = "Escherichia coli"
protein_name = "lacZ"

# Form the search query
query = {
    "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.taxonomy_lineage.name",
                    "operator": "exact_match",
                    "value": organism_name
                }
            },
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                    "operator": "exact_match",
                    "value": protein_name
                }
            }
        ]
    },
    "return_type": "entry"
}

# Send the search query to the API
response = requests.post(SEARCH_URL, json=query)
if response.status_code != 200:
    raise Exception(f"Failed to retrieve data: {response.text}")

# Parse the response to get PDB IDs
pdb_ids = [hit["identifier"] for hit in response.json()["result_set"]]

# Download the PDB files
for pdb_id in pdb_ids:
    pdb_response = requests.get(f"{DOWNLOAD_URL}/{pdb_id}.pdb")
    if pdb_response.status_code == 200:
        os.makedirs('pdb_files', exist_ok=True)
        # Save the PDB file
        with open(f'pdb_files/{pdb_id}.pdb', 'wb') as file:
            file.write(pdb_response.content)
        print(f'PDB file {pdb_id}.pdb has been downloaded.')
    else:
        print(f'Failed to download PDB file {pdb_id}.pdb.')
print('Download completed')
