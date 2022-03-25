def filter_genes_cluster(matrix, metadata, min_perc_cell, cell_column_name):
    """
    This function takes as input a normalized count matrix and for each gene it calculates the 
    percentage of cells expressing it (anything but 0). Then sets to 0 the expression of a given
    gene for al cells of a specific cell type if this gene is expressed in less than min_per_cell.
    
    Parameters
    ----------
    matrix: normalized gene expression matrix (genes x barcodes).
    metadata: index contains the barcode id and a single column named 'cell_type' indicating the group/cell type at which the barcode belongs to.
    min_perce_cell : percentage of cells requiring to express a given gene.
    """
    matrix = matrix.copy()
    
    # -- Cell types present in Metadata
    cell_type_data = set(metadata[cell_column_name])

    # -- Iterate over each cell type
    # --
    for cell_type in cell_type_data:

        # -- Obtain the barcode of the cells annotated under the current cell_type
        idx = metadata[cell_column_name] == cell_type
        cell_barcode = list(metadata.index[idx])

        # -- Calculate percentage of cells expressing the gene (expression !=0)
        gene_expr_perc = (matrix[cell_barcode] != 0).sum(axis = 1) / len(cell_barcode)

        # -- List of genes lowly expressed
        gene_lowly_expr = list(gene_expr_perc.index[gene_expr_perc < min_perc_cell])

        # -- Set to zero those genes that are expressed in a given cell type below the 
        # -- user defined min_perc_cell
        matrix.loc[gene_lowly_expr, cell_barcode] = 0
        
        
    # -- Return filtered matrix
    # --
    return(matrix)

def mean_expression_cluster(matrix, metadata, cell_column_name):
    """
    This functions calculates the mean expression of each gene per group/cell type.
    """
    import pandas as pd 

    matrix = matrix.copy()
    out_dict = {}

    # -- Cell types present in Metadata
    cell_type_data = set(metadata[cell_column_name])

    # -- Iterate over each cell type
    # --
    for cell_type in cell_type_data:

        # -- Obtain the barcode of the cells annotated under the current cell_type
        idx = metadata[cell_column_name] == cell_type
        cell_barcode = list(metadata.index[idx])

        # -- Calculate mean expression per cluster
        out_dict[cell_type] = matrix[cell_barcode].mean(axis = 1)

    # -- Convert the dictionary to a dataframe
    matrix_mean_expr = pd.DataFrame.from_dict(out_dict)
    
    
    # -- Return mean expression matrix
    # --
    return(matrix_mean_expr)



def scale_expression(matrix, upper_range):
    """
    Scale the expression of the gene across all the cell types in the dataset.
    """
    
    from sklearn.preprocessing import MinMaxScaler
    import pandas as pd

    # -- Transpose matrix to apply scaling per row (i.e. scale across cell types)
    scaler = MinMaxScaler(feature_range = (0, upper_range)).fit(matrix.T)
    matrix_scaled = scaler.transform(matrix.T).T

    matrix_scaled = pd.DataFrame(matrix_scaled,
                                 index = matrix.index,
                                 columns = matrix.columns)

    # -- return scaled matrix
    return(matrix_scaled)


def heteromer_geometric_expression(matrix, cellphone_genes, cellphone_complex):
    
    from scipy.stats.mstats import gmean
    import numpy as np
    import pandas as pd

    matrix = matrix.copy()
    cellphone_complex = cellphone_complex.copy()
    cellphone_genes = cellphone_genes.copy()
    
    # -- Define geometric function
    # -- 
    def geomean(x):
        sub_values = list(x)
        sub_prod = np.prod(sub_values)
        geom = np.power(sub_prod, 1/len(sub_values))
        return(geom)


    # -- Subset the mean expression matrix to keep only the genes
    # -- employed in the CellPhoneDB
    print(matrix.shape)
    idx = [gene in list(cellphone_genes['gene_name']) for gene in matrix.index]
    matrix = matrix.loc[idx]
    print(matrix.shape)

    # -- Keep only complex_name and uniprot ID
    cellphone_complex = cellphone_complex.iloc[:,[0,1,2,3]]

    # -- Create a dictionary with uniprot to gene name correspondence
    # -- Map uniprot ID to gene id and then convert it to a dictionary
    uniprot_dict = dict(zip(cellphone_genes['uniprot'], cellphone_genes['gene_name']))

    cellphone_complex = cellphone_complex.stack().map(uniprot_dict).unstack()

    cellphone_complex_dict = dict(
        zip(cellphone_complex.index,
            cellphone_complex.values.tolist()
           ))

    # -- Iterate over the complexes to calculate the geometric mean of them
    # -- If any member of the subunit its not present in the mean matrix
    # -- the geometric mean expression is not calculated
    complex_geom_mean = dict()
    for complex_id in cellphone_complex_dict:

        subunits_list = cellphone_complex_dict[complex_id]

        # -- Remove nan values from heteromers
        subunits_list = [i for i in subunits_list if str(i) != 'nan']

        # -- Test if all the members of the subunit are present in the matrix
        # -- if true then calculate the geometic mean of the complex
        check_subunit = all([sub in matrix.index for sub in subunits_list])

        if(check_subunit):
            complex_geom_mean[complex_id] = matrix.loc[subunits_list,].apply(geomean, axis = 0)


    # -- Convert geometric mean dictionary to dataframe
    complex_geom_mean_df = pd.DataFrame.from_dict(complex_geom_mean,
                                                  orient = 'index')

    # -- Detect and remove genes that have the same name that a complex
    # -- Cases such as OSMR, LIFR, IL27, if I dont to so, two rows assigned
    # -- to the same genes will appear
    idx = [i not in complex_geom_mean_df.index for i in matrix.index]
    matrix = matrix.loc[idx]
    
    final_df = pd.concat([matrix, complex_geom_mean_df])

    # -- Return final df
    # -- 
    return(final_df)


#####################################################################################################
#####################################################################################################

def heteromer_minimum_expression(matrix, cellphone_genes, cellphone_complex):
    
    from scipy.stats.mstats import gmean
    import numpy as np
    import pandas as pd

    matrix = matrix.copy()
    cellphone_complex = cellphone_complex.copy()
    cellphone_genes = cellphone_genes.copy()
    
    # -- Define geometric function
    # -- 
    def min_subunit(x):
        return(min(list(x)))


    # -- Subset the mean expression matrix to keep only the genes
    # -- employed in the CellPhoneDB
    idx = [gene in list(cellphone_genes['gene_name']) for gene in matrix.index]
    matrix = matrix.loc[idx]


    # -- Keep only complex_name and uniprot ID
    cellphone_complex = cellphone_complex.iloc[:,[0,1,2,3]]

    # -- Create a dictionary with uniprot to gene name correspondence
    # -- Map uniprot ID to gene id and then convert it to a dictionary
    uniprot_dict = dict(zip(cellphone_genes['uniprot'], cellphone_genes['gene_name']))

    cellphone_complex = cellphone_complex.stack().map(uniprot_dict).unstack()

    cellphone_complex_dict = dict(
        zip(cellphone_complex.index,
            cellphone_complex.values.tolist()
           ))

    # -- Iterate over the complexes to get the subunit that is the least expressed
    # -- If any member of the subunit its not present in the mean matrix
    # -- the min value is not taken
    complex_geom_mean = dict()
    for complex_id in cellphone_complex_dict:

        subunits_list = cellphone_complex_dict[complex_id]

        # -- Remove nan values from heteromers
        subunits_list = [i for i in subunits_list if str(i) != 'nan']

        # -- Test if all the members of the subunit are present in the matrix
        # -- if true then keep the value of the least expressed gene
        check_subunit = all([sub in matrix.index for sub in subunits_list])

        if(check_subunit):
            complex_geom_mean[complex_id] = matrix.loc[subunits_list,].apply(min_subunit, axis = 0)


    # -- Convert geometric mean dictionary to dataframe
    complex_geom_mean_df = pd.DataFrame.from_dict(complex_geom_mean,
                                                  orient = 'index')

    # -- Detect and remove genes that have the same name that a complex
    # -- Cases such as OSMR, LIFR, IL27, if I dont to so, two rows assigned
    # -- to the same genes will appear
    idx = [i not in complex_geom_mean_df.index for i in matrix.index]
    matrix = matrix.loc[idx]
    
    final_df = pd.concat([matrix, complex_geom_mean_df])

    # -- Return final df
    # -- 
    return(final_df)

def score_product(matrix, cellphone_genes, cellphone_interactions):


    from itertools import combinations_with_replacement
    import numpy as np
    import pandas as pd

    matrix = matrix.copy()
    cellphone_genes = cellphone_genes.copy()
    cellphone_inter = cellphone_interactions.copy()

    # -- Convert the uniprot ids from the interaction file to gene
    # -- symbol and then replace them.
    uniprot_dict = dict(zip(cellphone_genes['uniprot'],
                            cellphone_genes['gene_name']
                           ))

    cellphone_inter = cellphone_inter.loc[:,['partner_a', 'partner_b']].stack().replace(uniprot_dict).unstack()
    cellphone_inter['id_cp_interaction'] = cellphone_inter.index

    # -- Create lists of the interacting LR (and the reverse) to keep in the lr_outer_long only those entries 
    # -- that are described in the cpdb interaction file
    cpdb_list_a = list(cellphone_inter['partner_a']+'|'+cellphone_inter['partner_b'])
    cpdb_list_b = list(cellphone_inter['partner_b']+'|'+cellphone_inter['partner_a'])
    cpdb_list_all = cpdb_list_a + cpdb_list_b


    # -- Calculate all cell_type combinations
    # -- Initialize score dictionary
    combinations_cell_types = list(combinations_with_replacement(matrix.columns, 2))
    score_dict = {}

    # -- Iterate over the combination of pairs of cell types
    counter = 0
    for cell_A, cell_B in combinations_cell_types:

        # -- outer product
        # -- rows: first element outer function
        # -- cols: second element outer function
        lr_outer = pd.DataFrame(np.outer(list(matrix[cell_A]),
                                         list(matrix[cell_B])),
                                index = matrix.index,
                                columns = matrix.index)
        # ----------
        # ----------
        # ----------

        upper_idx = np.triu(np.ones(lr_outer.shape)).astype(bool)
        upper_tri = lr_outer.where(upper_idx)

        lower_idx = np.tril(np.ones(lr_outer.shape)).astype(bool)
        lower_tri = lr_outer.where(lower_idx)

        # -- Upper triangle (mtx_up) represents communication between cell B to cell A
        # -- Lower triangle (mtx_df) represents communication between cell A to cell B
        mtx_up = upper_tri.stack().reset_index()
        mtx_dn = lower_tri.stack().reset_index()

        mtx_up.columns = [cell_A, cell_B, 'Score']
        mtx_dn.columns = [cell_A, cell_B, 'Score']

        lr_outer_long = pd.concat([mtx_up, mtx_dn])
        # ----------
        # ----------
        # ----------


        # -- This operation should be performed just once as all the matrices in should contain
        # -- the same genes in rows and columns.
        if(counter == 0):

            lr_list = list(lr_outer_long.iloc[:,0]+'|'+lr_outer_long.iloc[:,1])
            idx_interactions = [i in cpdb_list_all for i in lr_list]
            counter += 1



        lr_outer_long_filt = lr_outer_long.loc[idx_interactions].copy()


        # -- Add interaction ID
        lr_outer_sorted = ['-'.join(sorted([a, b])) for a, b in zip(lr_outer_long_filt.iloc[:,0], lr_outer_long_filt.iloc[:,1])]
        cpdb_sorted = ['-'.join(sorted([a, b])) for a,b in zip(cellphone_inter['partner_a'], cellphone_inter['partner_b'])]

        dict_cpdb_id = dict(zip(cpdb_sorted, cellphone_inter['id_cp_interaction']))
        lr_outer_long_filt['id_cp_interaction'] = [dict_cpdb_id[i] for i in lr_outer_sorted]


        # -- assign cell names
        # -- first column is the sender
        # -- second column the receiver
        # -- third column the score
        lr_outer_long_filt.columns = [cell_A, cell_B, 'Score', 'id_cp_interaction']
        lr_outer_long_filt = lr_outer_long_filt.reset_index(drop = True)

        # -- Save result in dictionary
        dict_name = '|'.join(sorted([cell_A, cell_B]))
        score_dict[dict_name] = lr_outer_long_filt

    return(score_dict)
