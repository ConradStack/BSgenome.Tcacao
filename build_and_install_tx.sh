
#!/bin/bash

# Build source package:
R CMD BUILD TxDb.Tcacao.Geneious.pub3i

# Copy to OneDrive for uploading
ls -l1t *.tar.gz | head -1 | xargs -I {} cp {} ~/OneDrive/code/R_packages/

# (optional) remove and reinstall 
R CMD REMOVE mmutils
R CMD INSTALL mmutils

