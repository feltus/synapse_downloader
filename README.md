# synapse_downloader

# Background

This script will download files from synapse.org.  This specific example will download MCF10A RNAseq datasets.

# Preparing for data download
#Get a list of Synapse IDs from https://www.synapse.org/Synapse:syn18486863/tables/ and place in a file.
#Example MCF10A RNAseq dataset list is included as 'mcf10a-lincs_rnaseq-datasets.txt'

#Set  up a token at synapse.org and place in ~/.synapseConfig
```
[authentication]
authtoken = <TOKEN>
```

#Install Synapse Data Transfer client.  You might want to install in a conda or venv environment.
```
pip install synapseclient
```
# Download the datasets 

```
python DL_SYNAPSE.py -f mcf10a-lincs_rnaseq-datasets.txt
```
