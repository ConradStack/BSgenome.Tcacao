
# load biomaRt package 
require(biomaRt)

## Global variables
biomart.host = "phytozome.jgi.doe.gov"  # The phytozome biomart url

# List all "marts" that are available on phytozome:
listMarts(host=biomart.host)

# List all datasets available:
listDatasets(useMart("phytozome_mart", host=biomart.host,path="/biomart/martservice"))

# Select the main one (named 'phytozome')
phytozome = useMart("phytozome_mart","phytozome",host=biomart.host,path="/biomart/martservice")

# List available attributes and filters for this dataset
(attro = listAttributes(phytozome))
(filters=listFilters(phytozome))
columns(phytozome)
keytypes(phytozome)
cbind(sort(keytypes(phytozome)))

# get organism ids and names
oids=keys(phytozome, keytype="organism_id")

