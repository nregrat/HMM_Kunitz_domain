{
 "cells": [
  {
   "cell_type": "markdown",
   "id": "3cf3253b",
   "metadata": {},
   "source": [
    "# Data Collection"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "e296863b",
   "metadata": {},
   "source": [
    "Kunitz domain proteins from the PDB were retrieved, filtered and analyzed. Filtering was carried oout by sequence length, resolution, Pfam belonging and duplicated removal."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "b1b2f181",
   "metadata": {},
   "source": [
    "## Data Retrieval and Initial Processing"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "597c367a",
   "metadata": {},
   "source": [
    "Sequence accession numbers was retrieve from an advance search in PDB filtering for proteins annotated as Kunitz domains by Pfam (PF00014), with a resolution lower a 3.5 to ensure high quality experiment and a sequence length between 40 and 80 residues long, to avoid fragmented sequences.\n",
    "\n",
    "```Bash\n",
    "Query Summary:\n",
    "Attributes\n",
    "            Identifier = \"PF00014\"\n",
    "            AND\n",
    "            Annotation Type = \"Pfam\"\n",
    "      AND\n",
    "      Data Collection Resolution < 3.5\n",
    "      AND\n",
    "      Polymer Entity Sequence Length < 80\n",
    "      AND\n",
    "      Polymer Entity Sequence Length > 40\n",
    "```\n",
    "The ids were downloaded and save in the file `Files\\ids_pdb.csv`"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "b3b5bce1",
   "metadata": {},
   "source": [
    "Protein structure and other features were retrieved from the PDB via its GraphQL API. The obtain JSON was process to create a Pandas DataFrame containing relevant attributes for further analysis. The script filters for sequences between 40 and 80 residues and includes only those with Pfam annotations."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "374c7fd6",
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "import requests      # Import requests for HTTP requests to the PDB API"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "id": "fa1237c9",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>PDB_ID</th>\n",
       "      <th>Sequence</th>\n",
       "      <th>Length</th>\n",
       "      <th>Resolution</th>\n",
       "      <th>Annotation_type</th>\n",
       "      <th>Chain</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1AAL</td>\n",
       "      <td>RPDFCLEPPYTGPCRLRIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.85</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>F</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>1AAP</td>\n",
       "      <td>RPDFCLEPPYTGPCLARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.35</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>I</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>1B0C</td>\n",
       "      <td>VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFD...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.80</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>B</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>1BHC</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.22</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>A</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>4</th>\n",
       "      <td>1BPI</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVGGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.80</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>A</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "</div>"
      ],
      "text/plain": [
       "  PDB_ID                                           Sequence  Length  \\\n",
       "0   1AAL  RPDFCLEPPYTGPCRLRIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "1   1AAP  RPDFCLEPPYTGPCLARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "2   1B0C  VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFD...      58   \n",
       "3   1BHC  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "4   1BPI  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVGGGCRAKRNNFK...      58   \n",
       "\n",
       "   Resolution Annotation_type Chain  \n",
       "0        1.85            Pfam     F  \n",
       "1        1.35            Pfam     I  \n",
       "2        1.80            Pfam     B  \n",
       "3        1.22            Pfam     A  \n",
       "4        1.80            Pfam     A  "
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    " #Prepare PDB IDs for the GraphQL query\n",
    "ids_str = ''  # Initialize an empty string to build the list of IDs for the query\n",
    "ids = pd.Series(pd.read_csv('Files/ids_pdb.csv', header=1)['Entry ID'])  # Read IDs into a pandas Series\n",
    "for i in range(ids.shape[0]):\n",
    "    ids_str += '\"' + ids[i] + '\",'  # Append each ID in quotes, separated by commas\n",
    "ids_str = ids_str[:-1]  # Remove the trailing comma\n",
    "\n",
    "# Define the GraphQL query to fetch data for the specified PDB entries\n",
    "query = '''\n",
    "{\n",
    "  entries(entry_ids: [%s])\n",
    "  {\n",
    "    rcsb_id\n",
    "    rcsb_entry_container_identifiers {\n",
    "      entry_id\n",
    "    }\n",
    "    rcsb_entry_info {\n",
    "      diffrn_resolution_high {\n",
    "        value\n",
    "      }\n",
    "    }\n",
    "    polymer_entities {\n",
    "      entity_poly {\n",
    "        pdbx_seq_one_letter_code_can\n",
    "        rcsb_sample_sequence_length\n",
    "      }\n",
    "      polymer_entity_instances {\n",
    "        rcsb_polymer_entity_instance_container_identifiers {\n",
    "          auth_asym_id\n",
    "        }\n",
    "      }\n",
    "      rcsb_polymer_entity_annotation {\n",
    "        type\n",
    "      }\n",
    "      rcsb_polymer_entity_container_identifiers {\n",
    "        reference_sequence_identifiers {\n",
    "          database_accession\n",
    "        }\n",
    "      }\n",
    "    }\n",
    "  }\n",
    "} ''' %ids_str #%s will be replaced with the ids_str\n",
    "\n",
    "# Send the GET request to the PDB GraphQL endpoint\n",
    "p = requests.get('https://data.rcsb.org/graphql?query=%s' % requests.utils.requote_uri(query))\n",
    "\n",
    "# Parse the JSON response\n",
    "j = p.json()\n",
    "pdb = j['data']['entries']  # Extract the list of entries\n",
    "\n",
    "# Initialize an empty list to store dictionaries of extracted data\n",
    "df_pdb = []\n",
    "\n",
    "# Loop through each PDB entry\n",
    "for i in range(ids.shape[0]):\n",
    "    resolution = float(pdb[i]['rcsb_entry_info']['diffrn_resolution_high']['value'])  # Extract resolution as float\n",
    "    for j_entity in pdb[i]['polymer_entities']:  # Loop through polymer entities in the entry\n",
    "        len_seq = int(j_entity['entity_poly']['rcsb_sample_sequence_length'])  # Get sequence length\n",
    "        if len_seq > 80 or len_seq < 40:  # Filter sequences outside 40-80 residues\n",
    "            continue\n",
    "        seq = j_entity['entity_poly']['pdbx_seq_one_letter_code_can']  # Extract sequence\n",
    "        annot_type = j_entity['rcsb_polymer_entity_annotation'][0]['type']  # Extract annotation type\n",
    "        chain = j_entity['polymer_entity_instances'][0]['rcsb_polymer_entity_instance_container_identifiers']['auth_asym_id']  # Extract chain ID\n",
    "\n",
    "        # Append a list of attributes to the df list\n",
    "        df_pdb.append([ids[i], seq, len_seq, resolution, annot_type, chain])\n",
    "\n",
    "# Create a Pandas DataFrame from the dictionary\n",
    "df_pdb = pd.DataFrame(df_pdb, columns=['PDB_ID', 'Sequence', 'Length', 'Resolution', 'Annotation_type', 'Chain'])\n",
    "display(df_pdb.head())"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "id": "e32e5939",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "The obtained PDB entries were: 140\n"
     ]
    }
   ],
   "source": [
    "print( \"The obtained PDB entries were: \" + str(df_pdb.shape[0]))"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "7c488d99",
   "metadata": {},
   "source": [
    "The DataFrame (`df_pdb`) containing key protein attributes: *PDB ID*, *sequence*, *length*, *resolution*, *annotation type*, and *chain identifier*. This DataFrame serves as the foundational dataset for subsequent filtering and analysis steps within the workflow."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "657bdd58",
   "metadata": {},
   "source": [
    "## Removing Duplicate Entries and Filtering by Annotation Type\n",
    "\n",
    "The duplicate entries were remove and only those associated with 'Pfam' annotations were retained. Duplicates can arise from various sources, such as multiple chains within the same PDB entry or redundant data, and their removal is crucial for maintaining data integrity. By specifically filtering for 'Pfam' annotations, we ensure that the subsequent analysis is based exclusively on relevant Kunitz domains, thereby enhancing the specificity and accuracy of our study."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "id": "8e6d01db",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Found 6 PDB_IDs with duplicate entries:\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "0    1D0D\n",
       "1    1EAW\n",
       "2    1ZJD\n",
       "3    3LDJ\n",
       "4    3LDM\n",
       "5    5XX7\n",
       "dtype: str"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# Identify PDB_IDs that appear more than once in the DataFrame\n",
    "duplicated_pdb_entries = df_pdb[df_pdb.duplicated(subset=['PDB_ID'], keep=False)]\n",
    "\n",
    "# Extract the unique PDB_IDs that are duplicated\n",
    "duplicated_pdb_ids = duplicated_pdb_entries['PDB_ID'].unique()\n",
    "\n",
    "print(f\"Found {len(duplicated_pdb_ids)} PDB_IDs with duplicate entries:\")\n",
    "display(pd.Series(duplicated_pdb_ids))"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "6d8e4b12",
   "metadata": {},
   "source": [
    "As observed, some entries were redundantly annotated by both Gene Ontology (GO) and Pfam. For consistency and relevance to Kunitz domains, only those entries specifically annotated by Pfam were retained."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "id": "d946804f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/html": [
       "<div>\n",
       "<style scoped>\n",
       "    .dataframe tbody tr th:only-of-type {\n",
       "        vertical-align: middle;\n",
       "    }\n",
       "\n",
       "    .dataframe tbody tr th {\n",
       "        vertical-align: top;\n",
       "    }\n",
       "\n",
       "    .dataframe thead th {\n",
       "        text-align: right;\n",
       "    }\n",
       "</style>\n",
       "<table border=\"1\" class=\"dataframe\">\n",
       "  <thead>\n",
       "    <tr style=\"text-align: right;\">\n",
       "      <th></th>\n",
       "      <th>PDB_ID</th>\n",
       "      <th>Sequence</th>\n",
       "      <th>Length</th>\n",
       "      <th>Resolution</th>\n",
       "      <th>Annotation_type</th>\n",
       "      <th>Chain</th>\n",
       "    </tr>\n",
       "  </thead>\n",
       "  <tbody>\n",
       "    <tr>\n",
       "      <th>0</th>\n",
       "      <td>1AAL</td>\n",
       "      <td>RPDFCLEPPYTGPCRLRIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.85</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>F</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>1</th>\n",
       "      <td>1AAP</td>\n",
       "      <td>RPDFCLEPPYTGPCLARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.35</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>I</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>2</th>\n",
       "      <td>1B0C</td>\n",
       "      <td>VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFD...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.80</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>B</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>3</th>\n",
       "      <td>1BHC</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.22</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>A</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>4</th>\n",
       "      <td>1BPI</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVGGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.80</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>A</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>...</th>\n",
       "      <td>...</td>\n",
       "      <td>...</td>\n",
       "      <td>...</td>\n",
       "      <td>...</td>\n",
       "      <td>...</td>\n",
       "      <td>...</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>130</th>\n",
       "      <td>7QIR</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.09</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>A</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>131</th>\n",
       "      <td>7QIS</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>2.63</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>J</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>132</th>\n",
       "      <td>7QIT</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>2.80</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>E</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>133</th>\n",
       "      <td>8PTI</td>\n",
       "      <td>VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFD...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.50</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>B</td>\n",
       "    </tr>\n",
       "    <tr>\n",
       "      <th>134</th>\n",
       "      <td>9PTI</td>\n",
       "      <td>RPDFCLEPPYTGPCKARIIRYFYNAKAGLVQTFVYGGCRAKRNNFK...</td>\n",
       "      <td>58</td>\n",
       "      <td>1.60</td>\n",
       "      <td>Pfam</td>\n",
       "      <td>B</td>\n",
       "    </tr>\n",
       "  </tbody>\n",
       "</table>\n",
       "<p>135 rows × 6 columns</p>\n",
       "</div>"
      ],
      "text/plain": [
       "    PDB_ID                                           Sequence  Length  \\\n",
       "0     1AAL  RPDFCLEPPYTGPCRLRIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "1     1AAP  RPDFCLEPPYTGPCLARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "2     1B0C  VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFD...      58   \n",
       "3     1BHC  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "4     1BPI  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVGGGCRAKRNNFK...      58   \n",
       "..     ...                                                ...     ...   \n",
       "130   7QIR  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "131   7QIS  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "132   7QIT  RPDFCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFK...      58   \n",
       "133   8PTI  VREVCSEQAETGPCRAMISRWYFDVTEGKCAPFFYGGCGGNRNNFD...      58   \n",
       "134   9PTI  RPDFCLEPPYTGPCKARIIRYFYNAKAGLVQTFVYGGCRAKRNNFK...      58   \n",
       "\n",
       "     Resolution Annotation_type Chain  \n",
       "0          1.85            Pfam     F  \n",
       "1          1.35            Pfam     I  \n",
       "2          1.80            Pfam     B  \n",
       "3          1.22            Pfam     A  \n",
       "4          1.80            Pfam     A  \n",
       "..          ...             ...   ...  \n",
       "130        1.09            Pfam     A  \n",
       "131        2.63            Pfam     J  \n",
       "132        2.80            Pfam     E  \n",
       "133        1.50            Pfam     B  \n",
       "134        1.60            Pfam     B  \n",
       "\n",
       "[135 rows x 6 columns]"
      ]
     },
     "execution_count": 8,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "# Filter the DataFrame to keep only entries with 'Pfam' annotation type.\n",
    "df_pdb = df_pdb.query('Annotation_type == \"Pfam\"').reset_index(drop=True)\n",
    "df_pdb"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "id": "d48fa3c1",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "The obtained PDB entries were: 135\n"
     ]
    }
   ],
   "source": [
    "print( \"The obtained PDB entries were: \" + str(df_pdb.shape[0]))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "id": "6c9a9a05",
   "metadata": {},
   "outputs": [],
   "source": [
    "# Save the dataframe as a TSV file in the Files folder\n",
    "df_pdb.to_csv(\"Files/df_pdb.tsv\", sep=\"\\t\", index=False)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.13.14"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
