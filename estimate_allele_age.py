#!/usr/bin/python
import numpy as np
import pandas as pd
from scipy.stats import norm
import matplotlib.pyplot as plt
from tqdm import tqdm
import sys
import argparse


def allele_frequency_dist_with_uncertainty(p, n):
    """
    Compute confidence interval for allele frequency
    :param p: float, measured allele frequency
    :param n: int, sample size (individuals)
    :return: initialized normal distribution, lower bound (float), upper bound (float)
    """
    # define std of allele frequency dist
    std = np.sqrt((p * (1 - p)) / (2 * n))
    # initialize pdf
    pdf = norm(loc=p, scale=std)
    # define lower bound --> cannot be lower than 1/2n
    lower_bound = p - 4 * std
    if lower_bound < 1 / (2 * n):
        lower_bound = 1 / (2 * n)
    # define upper bound
    upper_bound = p + 4 * std
    if upper_bound > 1.0:
        upper_bound = 1.0
    # plot PDF
    rvs = pdf.rvs(100000)
    plot_pdf_allele_frequency(rvs, p, lower_bound, upper_bound)
    return pdf, lower_bound, upper_bound


def get_age_distribution(p, sample_size):
    """
    Estimate age distribution in terms of 2N and propagate uncertainty in allele frequency
    :param p: float, allele frequency
    :param sample_size: int, sample size
    :return: array, array: likelihood of age, age
    """
    step_size = 0.0001
    # get pdf of allele frequency and bounds of CI
    pdf_allele_frequency, lower_bound, upper_bound = allele_frequency_dist_with_uncertainty(p, sample_size)
    required_tmax = lambda p: (2 * (-np.log(1 - p) + (2 * sample_size) * np.log(1 - p) - np.log(1))) / \
                              ((2 * sample_size) * (np.log(1 - p) + np.log(1)))
    p = lower_bound
    age_prob = []
    age = []
    # iterate over different values for allele frequency and determine pdf of allele age
    while p < upper_bound:

        t = np.arange(0, required_tmax(p) + 0.00001, 0.000001)
        cum_dist = cumulative_dist_age(p, sample_size, t)
        pdf_t = np.diff(cum_dist)
        pdf_dist = pdf_t * pdf_allele_frequency.pdf(p)
        age_prob.append(pdf_dist)
        p += step_size
    # put pdf of allele ages into array
    age_prob_arr = np.zeros((len(age_prob), max(map(len, age_prob))))
    t = np.arange(0, age_prob_arr.shape[1] * 0.000001, 0.000001)
    for i in range(len(age_prob)):
        age_prob_arr[i, :len(age_prob[i])] += age_prob[i]

    # get probability associated with each age --> sum over all frequency and normalize so that they sum to one
    age_prob_arr = age_prob_arr.sum(axis=0) / age_prob_arr.sum(axis=0).sum()


    return age_prob_arr, t


def cumulative_dist_age(p, n, t):
    """
    Cumulative distribution of allele age as defined by Slatkin and Rannala
    :param p: float, allele frequency
    :param n: int, sample size
    :param t: float (array), age in 2N
    :return: float (array), likelihood
    """
    # cumulative distribution of allele ages according to Slatkin & Rannala
    return (1 - p) ** (-1 + (2 * n) / (1 + (2 * n) * t / 2))


def expected_age(p, n):
    """
    MLE solution for allele age
    :param p: float, allele frequency
    :param n: int, sample size
    :return: float, expected allele age
    """
    # MLE solution for allele age
    return -np.log(1 - p) - 2 / (2 * n)


def plot_pdf_age(t, probs):
    """
    Plot PDF of allele age in terms of 2N
    :param t: array, age
    :param probs: array, likelihoods
    """
    fig, ax = plt.subplots()
    ax.plot(t, probs)
    ax.set_xlabel("time in 2N")
    ax.set_ylabel("Likelihood")
    ax.set_xlim([0, 0.2])
    fig.savefig('pdf_allele_age.png', bbox_inches='tight')
    plt.close()


def plot_pdf_allele_frequency(freqs, p, lower_bound, upper_bound):
    """
    Plot PDF of allele frequencys
    :param freqs: array, rvs samples from allele frequency PDF
    """
    fig, ax = plt.subplots()
    ax.hist(freqs, bins=1000, weights=np.full(freqs.shape[0], 1 / freqs.shape[0]), histtype='step')
    ax.set_xlabel('p')
    ax.set_ylabel('Likelihood')
    ax.set_xlim([0, ax.get_xlim()[1]])
    ax.axvline(p, ls='-', color='k', label="Measured frequency")
    ax.axvline(lower_bound, ls='--', color='k',)
    ax.axvline(upper_bound, ls='--', color='k', label="95% CI")
    fig.savefig('pdf_allele_frequency.png', bbox_inches='tight')
    plt.close()


def plot_pdf_population_size(rvs):
    """
    Plot PDF of effective population size
    :param rvs: array, estimated population sizes
    """
    fig, ax = plt.subplots()
    ax.hist(rvs, bins=10, weights=np.full(rvs.shape[0], 1 / rvs.shape[0]), histtype='step')
    ax.set_xlabel(r'$N_e$')
    ax.set_ylabel('Likelihood')
    fig.savefig('pdf_population_size.png', bbox_inches='tight')
    plt.close()


def main(argv):
    parser = argparse.ArgumentParser(description='Script for estimating the age of the HOXB13 variant found in African men')
    parser.add_argument('-p', '--allele_frequency', help='Allele frequency; default=0.0037', type=float, default=0.0037)
    parser.add_argument('-n', '--sample_size', help='Number of individuals sampled from population; default=1071',
                        type=int, default=1071)
    parser.add_argument('-b', '--bantu_migration', help='Date of Bantu migration in years ago; default=3000',
                        type=int, default=3000)
    parser.add_argument('-u', '--sample_size_uganda', help='Number of individuals in Bantu speaking Ugandan cohort'
                                                           'without observing a single copy; default=677',
                        type=int, default=677)
    parser.add_argument('-g', '--generation_time', help='Generation time in years; default=25', type=int, default=25)
    args = parser.parse_args()
    # read heterozygosity and Ne results
    heterozygosity_df = pd.read_csv("data/heterozygous_counts_plink2.tab", sep='\t', header=0)
    heterozygosity_df.Ne = heterozygosity_df.Ne.astype(int)
    plot_pdf_population_size(heterozygosity_df.Ne.values)
    # Ne values
    Ne = heterozygosity_df.Ne.value_counts().index.values
    # PDF of Ne
    Ne_dist = norm(loc=int(Ne.mean()), scale=int(Ne.std()))

    # get distribution of allele ages P(age) = sum(P(age | frequency) * P(frequency))
    t_dist, t = get_age_distribution(args.allele_frequency, args.sample_size,)
    # plot PDF of ages in 2N
    plot_pdf_age(t, t_dist)
    # run simulations
    nr_simulations = 1000000
    # sample population size
    sampled_N = Ne_dist.rvs(nr_simulations)
    # sample age
    sampled_t = np.random.choice(t, p=t_dist, size=nr_simulations)
    # scale age to two 2N and generation time --> years
    final_age_dist = sampled_t * 2 * sampled_N * args.generation_time
    # get histogram --> probability mass
    hist_unnormalized, bins = np.histogram(final_age_dist, bins=np.arange(0, final_age_dist.max() + 10, 10),
                                           weights=np.full(final_age_dist.shape[0], 1 / final_age_dist.shape[0]))
    # bin centers
    bins = (bins[1:] + bins[:-1]) / 2

    pdf_allele_frequency, _, _ = allele_frequency_dist_with_uncertainty(args.allele_frequency,
                                                                        args.sample_size)
    # since it is uniform I don't need to sample it
    admixture = np.arange(0.5, 0.75, 0.01)
    final_hist = []
    for adm in tqdm(admixture):
        hist = []
        # in total sample 20k replicates for each admixture proportion
        for i in range(10):
            # sample allele frequencies at time of Bantu migration
            p_B = ((bins - args.bantu_migration) / bins)[:, np.newaxis] * pdf_allele_frequency.rvs(2000)
            # probability of not observing copy by chance
            pre_bantu = (1 - p_B * adm) ** (2 * args.sample_size_uganda)
            pre_bantu[bins < args.bantu_migration] = 1
            # combine with allele age distribution
            c_hist = hist_unnormalized[:, np.newaxis] * pre_bantu
            # get the probability per bin
            c_hist = c_hist.sum(axis=1) / c_hist.sum(axis=1).sum()
            hist.append(c_hist)
        # get probability per admixture proportion
        hist = np.array(hist).T
        hist = hist.sum(axis=1) / hist.sum(axis=1).sum()
        final_hist.append(hist)
    # get final probability distribution
    final_hist = np.array(final_hist).T
    final_hist = final_hist.sum(axis=1) / final_hist.sum(axis=1).sum()
    # plot
    fig, ax = plt.subplots()
    ax.plot(bins, final_hist, color='black')
    ax.axvline(args.bantu_migration, ls='--', c='black', label='Time of Bantu migration')
    ax.set_xlabel('Age in years')
    ax.set_ylabel('Probability mass')
    ax.set_xlim([0, 20000])
    ax.legend()
    fig.savefig('age_distribution.png', bbox_inches='tight')
    plt.close()
    print("Not considering Bantu migration")
    print('Median: {:.2f} years'.format(bins[np.where(np.cumsum(hist_unnormalized) >= 0.5)[0][0]]))
    print('Mean: {:.2f} years'.format(np.sum(hist_unnormalized * bins)))
    print('Mode: {:.2f} - {:.2f} years'.format(bins[np.argmax(hist_unnormalized)], bins[np.argmax(hist_unnormalized) + 1]))
    print('95% CI: {:.2f} - {:.2f} years'.format(bins[np.where(np.cumsum(hist_unnormalized) >= 0.025)[0][0]],
                                                 bins[np.where(np.cumsum(hist_unnormalized) >= 0.975)[0][0]]))
    print("Considering Bantu migration:")
    print('Median: {:.2f} years'.format(bins[np.where(np.cumsum(final_hist) >= 0.5)[0][0]]))
    print('Mean: {:.2f} years'.format(np.sum(final_hist * bins)))
    print('Mode: {:.2f} - {:.2f} years'.format(bins[np.argmax(final_hist)], bins[np.argmax(final_hist) + 1]))
    print('95% CI: {:.2f} - {:.2f} years'.format(bins[np.where(np.cumsum(final_hist) >= 0.025)[0][0]],
                                                 bins[np.where(np.cumsum(final_hist) >= 0.975)[0][0]]))


if __name__ == '__main__':
    main(sys.argv[1:])
