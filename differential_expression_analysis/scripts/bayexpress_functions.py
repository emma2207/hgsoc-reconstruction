# Please read, follow and use our LICENSE and CITATION.

import scipy.special as sc
import numpy as np


# Optional hyper-parameters (here we use 1 for a flat prior)
u_1 = 1
u_2 = 1


# Calculating log_10 Bayes factors for differential gene expression (BF_21)
def get_BF_21(
    N_1,  # Total number of reads for genes in experiment 1
    n_1,  # Number of reads for a gene in experiment 1
    N_2,  # Total number of reads for genes in experiment 2
    n_2,  # Number of reads for a gene in experiment 2
):

    # betaln returns the natural log of the beta function
    # hence we after convert natural log to log_10
    BF_21 = (
        sc.betaln(u_1 + n_1, u_2 + N_1 - n_1)
        + sc.betaln(u_1 + n_2, u_2 + N_2 - n_2)
        - sc.betaln(u_1 + n_1 + n_2, u_2 + N_1 - n_1 + N_2 - n_2)
    ) / np.log(10)

    # if we wanted to return nans for genes with 0 in all replicates and conditions we can use this
    # BF_21 = np.where((n_1+n_2==0), np.nan, BF_21)

    return BF_21


# ratio of expression
# calculating log_2 fold change
def get_FC(
    N_1,  # Total number of reads for genes in experiment 1
    n_1,  # Number of reads for a gene in experiment 1
    N_2,  # Total number of reads for genes in experiment 2
    n_2,  # Number of reads for a gene in experiment 2
):

    rate_1 = (u_1 + n_1) / (u_2 + N_1 - n_1)
    rate_2 = (u_1 + n_2) / (u_2 + N_2 - n_2)

    FC = np.log2(rate_2 / rate_1)

    # if we wanted to return nans for genes with 0 in all replicates and conditions we can use this
    # FC = np.where((n_1+n_2==0), np.nan, FC)

    return FC


# calculating q (following Laplace's rule of succession)
# assuming u_1 = u_2 = 1
def get_q(
    n, N  # Number of reads mapping to a gene  # Total number of reads in an experiment
):

    return (n + 1) / (N + 2)


# calculating log_10 Bayes factors (BF_k1) for consistency
def get_BF_k1(data):

    evidence2 = np.full(len(data), 0)

    for col in data.columns:
        n_j = data[col]
        N_j = n_j.sum()
        evidence2 += sc.betaln(u_1 + n_j, u_2 + N_j - n_j)

    n_i = data.sum(axis=1, numeric_only=True)
    N = n_i.sum()

    evidence1 = sc.betaln(u_1 + n_i, u_2 + N - n_i)

    # convert natural log to log_10
    return (evidence2 - evidence1) / np.log(10)
