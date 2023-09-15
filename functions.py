import numpy as np
from scipy.stats import pearsonr
from msepm import helpers

def pearson_correlation(meth_matrix: np.array, phenotype: np.array) -> np.array:
    """calculate pearson correlation coefficient between rows of input matrix and phenotype"""
    # calculate mean for each row and phenotype mean
    matrix_means = np.mean(meth_matrix, axis=1)
    phenotype_mean = np.mean(phenotype)
    

    # subtract means from observed values
    transformed_matrix = meth_matrix - matrix_means.reshape([-1,1])
    transformed_phenotype = phenotype - phenotype_mean

    # calculate covariance
    covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
    variance_meth = np.sqrt(np.sum(transformed_matrix ** 2, axis=1))
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))

    return covariance / (variance_meth * variance_phenotype)

#PARTIAL CORRELATION ORDER 1 WITH 2D ARRAYS
def partial_correlation_1(meth_matrix: np.array, phenotype: np.array, phenotypeCo1: np.array) -> np.array:
    pcorr_12 = pearson_correlation(meth_matrix, phenotype)
    pcorr_13 = pearson_correlation(meth_matrix, phenotypeCo1)
    pcorr_23 = pearsonr(phenotype, phenotypeCo1)
    pcorr_23 = pcorr_23[0]
#     print(pcorr_12,pcorr_13,pcorr_23)
    numerator = (pcorr_12 - (pcorr_13*pcorr_23))
    denominator = np.sqrt((1-(pcorr_13**2)*(1-(pcorr_23)**2)))
    pcorr_12_3 = numerator / denominator

    return pcorr_12_3

#PARTIAL CORRELATION ORDER 2 WITH 2D ARRAYS
def partial_correlation_2(meth_matrix: np.array, phenotype: np.array, phenotypeCo1: np.array, phenotypeCo2: np.array) -> np.array:
    pcorr_12_3 = partial_correlation_1(meth_matrix, phenotype, phenotypeCo1)
    pcorr_14_3 = partial_correlation_1(meth_matrix, phenotypeCo2, phenotypeCo1)
    pcorr_24_3 = partial_correlation_single(phenotype, phenotypeCo2, phenotypeCo1)
    
    numerator = (pcorr_12_3 - (pcorr_14_3*pcorr_24_3))
    denominator = np.sqrt((1-(pcorr_14_3**2)*(1-(pcorr_24_3)**2)))
    pcorr_12_34 = numerator / denominator
    return pcorr_12_34

#PARTIAL CORRELATION ORDER 1 FOR 1-D ARRAYS
def partial_correlation_single(trait, traitCo1, traitCo2):
    pcorr_12 = pearsonr(trait, traitCo1)
    pcorr_13 = pearsonr(trait, traitCo2)
    pcorr_23 = pearsonr(traitCo1, traitCo2)
    
    pcorr_12 = pcorr_12[0]
    pcorr_13 = pcorr_13[0]
    pcorr_23 = pcorr_23[0]

    numerator = (pcorr_12 - (pcorr_13*pcorr_23))
    denominator = np.sqrt((1-(pcorr_13**2)*(1-(pcorr_23)**2)))
    pcorr_12_3 = numerator / denominator

    return pcorr_12_3

def matrix_of_coeff(meth, pheno):
    methylation_values = np.array(meth)
    pcc_coefficients = pearson_correlation(methylation_values, pheno)
    return pcc_coefficients

def filter_sites(group, corr_pheno, cutoff):
    methylation_values = np.array(group.T)
    #counts = np.array(df_counts.loc[group.index].select_dtypes(exclude = ['object']).T)

    #print('all samples and sites', np.shape(methylation_values))
    
    # Filter Out Bad Samples (many NaNs or no age information)
#     bad_samples = np.where(sum(np.isnan(methylation_values))>5000)[0]
#     methylation_values = np.delete(methylation_values, bad_samples, axis = 1)
#     corr_pheno = np.delete(corr_pheno, bad_samples)
#     #counts = np.delete(counts, bad_samples, axis = 1)
#     print('removed bad samples', np.shape(methylation_values))
    
#     # Filter out sites with NaNs
#     methylation_values = methylation_values[~np.isnan(methylation_values).any(axis=1)]
#     print('removed sites with NaNs', np.shape(methylation_values))
    
    # Filter Sites with Low Coverage
#     counts = counts[~np.isnan(counts).any(axis=1)]
#     total_counts = np.sum(counts, axis = 1)
#     training_sites = np.where(total_counts>np.median(total_counts))[0]
#     methylation_values = methylation_values[training_sites, :]
#     print('removed sites with low coverage', np.shape(methylation_values))

    abs_pcc_coefficients = abs(pearson_correlation(methylation_values, corr_pheno))
    print(abs_pcc_coefficients)
    # histogram of correlations
    #plt.hist(abs_pcc_coefficients)
    #plt.show()
    
    # return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(abs_pcc_coefficients > np.percentile([abs_pcc_coefficients], cutoff))[0]
#     print('pcc cutoff:', np.percentile([abs_pcc_coefficients], cutoff) )
#     print('training sites:', np.shape(training_sites))
    return methylation_values[training_sites,:]

#http://wzimm.weebly.com/uploads/1/3/5/2/13522665/partial_correlation_intro_1.pdf
def filter_sites_partial(group, corr_pheno, corr_phenoCO, cutoff):
    methylation_values = np.array(group.T)

    abs_pcc_coefficients = abs(partial_correlation_1(methylation_values, corr_pheno, corr_phenoCO))
    
    # return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(abs_pcc_coefficients > np.percentile([abs_pcc_coefficients], cutoff))[0]
    return methylation_values[training_sites,:]

#http://wzimm.weebly.com/uploads/1/3/5/2/13522665/partial_correlation_intro_1.pdf
def filter_sites_partial_order2(group, corr_pheno, corr_phenoCO, corr_phenoCO2, cutoff):
    methylation_values = np.array(group.T)
    #print('all samples and sites', np.shape(methylation_values))
    abs_pcc_coefficients = abs(partial_correlation_2(methylation_values, corr_pheno, corr_phenoCO, corr_phenoCO2))
    
    # return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(abs_pcc_coefficients > np.percentile([abs_pcc_coefficients], cutoff))[0]
    #print('training sites:', len(training_sites))
    return methylation_values[training_sites,:]

def filter_sites_partial_order2_locs(group, corr_pheno, corr_phenoCO, corr_phenoCO2, cutoff, chr_locs):
    methylation_values = np.array(group.T)
    #print('all samples and sites', np.shape(methylation_values))
    abs_pcc_coefficients = abs(partial_correlation_2(methylation_values, corr_pheno, corr_phenoCO, corr_phenoCO2))
    
    # return list of site indices with a high absolute correlation coefficient
    training_sites = np.where(abs_pcc_coefficients > np.percentile([abs_pcc_coefficients], cutoff))[0]
    #print('training sites:', len(training_sites))
    
#     cols = df_meth.columns.values
#     chr_locs = cols[training_sites]
    
    return chr_locs[training_sites]

def filter_sites_multi(group, corr_pheno, cutoffs, traits, locsORvals = "vals"):
    
    methylation_values = np.array(group.select_dtypes(exclude = ['object']).T)
    #counts = np.array(df_counts.loc[group.index].select_dtypes(exclude = ['object']).T)
   
    print('all samples and sites', np.shape(methylation_values))
    
    # Filter Out Bad Samples (many NaNs or no age information)
#     bad_samples = np.where(sum(np.isnan(methylation_values))>5000)[0]
#     methylation_values = np.delete(methylation_values, bad_samples, axis = 1)
#     corr_pheno = np.delete(corr_pheno, bad_samples, axis = 1)
#     counts = np.delete(counts, bad_samples, axis = 1)
#     print('removed bad samples', np.shape(methylation_values))
    
#     # Filter out sites with NaNs
#     methylation_values = methylation_values[~np.isnan(methylation_values).any(axis=1)]
#     print('removed sites with NaNs', np.shape(methylation_values))
    
#     # Filter Sites with Low Coverage
#     counts = counts[~np.isnan(counts).any(axis=1)]
#     total_counts = np.sum(counts, axis = 1)
#     training_sites = np.where(total_counts>np.median(total_counts))[0]
#     methylation_values = methylation_values[training_sites, :]
#     print('removed sites with low coverage', np.shape(methylation_values))
    
    mdpc = abs(helpers.pearson_correlation(corr_pheno, methylation_values))
    
    training_sites = [[] for i in range(corr_pheno.shape[1])]
    
    for i in range(mdpc.shape[1]):
        site = mdpc[:,i]
        sorted_site = np.argsort(site)
        if site[sorted_site[-1]] > cutoffs[sorted_site[-1]]:
            training_sites[sorted_site[-1]].append(i)
                      
    training_sites = sum(training_sites, [])
    
    if locsORvals == 'locs':
        cols = group.columns.values
        chr_locs = cols[training_sites]
        return chr_locs
    else:
        print('training sites:', len(training_sites))
        return methylation_values[training_sites,:]
    
