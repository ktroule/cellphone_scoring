# Tentative code to score cellphonedb interactions
Code to rank cellphonedb interactions. \
This protcol intendes to score and rank interactions based on a score which ranges between 0 and 100. \
This score is relative to the cells present in the dataset, adding or excluding cell types will make the score to change at some degree. \


### Protocol steps.
The protocol consist of 5 (yet to be optimized). \
**Step 1**: filter genes expressed in less than `min_perc_cell` of cells in a given cluster. \
**Step 2**: calculate the gene's mean expression per cluster. \
**Step 3**: scale the gene's mean expression across clusters. \
**Step 4**: filter genes expressed in less than `min_perc_cell of cells in a given cluster. \
**Step 5**: calculate the ligand-receptor score (0-100) and cry. \
