import os
import argparse

import numpy as np
import pandas as pd
import allel
from scipy.stats import chi2
from scipy.stats import beta
from itertools import cycle


COLUMNS_TYPES = {
    "chrom": str,
    "pos": int,
    "beta": float,
    "se": float,
    "pval": float,
    "eaf": float,
    "maf": float,
    "imp_info": float
}


class SumstatsFile:

    def __init__(self, file_path:str, map_names):
        # input tsv file
        if file_path.endswith('tsv'):
            self.df = self.open_table_file(file_path, map_names, sep='\t')
        if file_path.endswith('csv'):
            self.df = self.open_table_file(file_path, map_names, sep=',')
        # input vcf file
        elif file_path.endswith('vcf') or file_path.endswith('vcf.gz'):
            self.df = self.open_vcf_file(file_path, map_names)

    def open_vcf_file(self, file_path, map_names):
        # adjust col mapping dict for VCF
        map_names_vcf = {
            'variants/CHROM': "chrom",
            'variants/POS': "pos",
            'INFO': "imp_info",
        }
        for k,v in map_names.items():
            if v not in ['chrom', 'pos', 'imp_info']:
                map_names_vcf[f'calldata/{k}'] = v

        # read vcf
        sumstat_dictionary = allel.read_vcf(file_path, fields=map_names_vcf.keys(), alt_number=1)

        # convert dictionary to dataframe
        df_vcf = pd.DataFrame.from_dict({k:v.flatten() for k,v in sumstat_dictionary.items()})

        # get right info column
        if 'variants/INFO' in df_vcf.columns and not df_vcf['variants/INFO'].isnull().all():
            df_vcf['INFO'] = df_vcf['variants/INFO']
        else:
            df_vcf['INFO'] = np.nan

        df_vcf = df_vcf.rename(columns=map_names_vcf)

        return df_vcf

    def open_table_file(self, file_path, map_names, sep='\t'):
        df = pd.read_table(file_path, sep=sep, dtype=str)
        if 'is_cc' in df.columns:
            df['is_cc'] = df['is_cc'].fillna(0).astype(int)  # Set to False
        df = df.rename(columns=map_names)

        return df


    def format_gwas(self, study_id):
        col_dtypes = { k:v for k,v in COLUMNS_TYPES.items() if k in self.df.columns }
        self.df = self.df.astype(col_dtypes)
       
        self.df['study_id'] = study_id

        if 'eaf' in self.df.columns:
            self.df['maf'] = 1 - self.df['eaf']

        # filter empty se and eaf rows - cannot be used
        self.df = self.df[self.df.se.notna()]
        self.df = self.df[self.df.maf.notna()]

        # ensure maf really is minor
        self.df['inv_maf'] = 1 - self.df['maf']
        self.df['maf'] = self.df[['maf', 'inv_maf']].min(axis=1)

        # only use colums that are needed
        cols_to_use = [ c for c in COLUMNS_TYPES.keys() if c != "eaf" ]
        self.df = self.df[cols_to_use]

    def quality_control(self, out_dir, study_id, n_tot, n_cas):
        self.__generate_frequency_curves__(out_dir, study_id)
        self.__generate_scatter_plots__(out_dir, study_id)
        self.__generate_qq_plot__(out_dir, study_id)
        self.__generate_manhattan_plot__(out_dir, study_id)
        self.__descriptive_statistics__(f"{out_dir}/{study_id}/descriptive_statistics.txt", n_tot, n_cas)
        self.__descriptive_statistics_table__(f"{out_dir}/{study_id}/descriptive_statistics.tsv")

    def __generate_frequency_curves__(self, out_dir, study_id):
        if not os.path.exists(f"{out_dir}/{study_id}"):
            os.makedirs(f"{out_dir}/{study_id}")
        
        data = pd.DataFrame()
        data['INFO'] = self.df['imp_info']
        data['MAF'] = self.df['maf']
        data['Absolute_beta'] = np.abs(self.df['beta'])

        for var in ['INFO', 'MAF', 'Absolute_beta']:
        # for var in ['INFO']:
            if not data[var].isnull().all():
                data[var][~np.isfinite(data[var])] = 0 #remove nan and infs
                count, bins_count = np.histogram(data[var], bins=50)
                pdf = count / sum(count)
                cdf = np.cumsum(pdf)
                df = pd.DataFrame({'x': np.round(bins_count[:-1], 3), 'cdf': cdf}, columns=['x', 'cdf'])
                ax = df.plot.line(y='cdf', figsize=(8,8))
                ax = df.plot.bar(x='x', y='cdf', alpha=0.8, xlabel=f"{var}", ylabel="Frequency", figsize=(8,8), ax=ax)
                ax.set_title(study_id)
                ax.get_figure().savefig(f"{out_dir}/{study_id}/frequencyCurve-{var}.png")

    def __generate_scatter_plots__(self, out_dir, study_id):
        if not os.path.exists(f"{out_dir}/{study_id}"):
            os.makedirs(f"{out_dir}/{study_id}")
        
        # add relevant data to data frame
        data = self.df[['imp_info', 'maf', 'pval']].copy()
        data['log_pval'] = -np.log10(self.df['pval'])
        data['abs_es'] = np.abs(self.df['beta'])
        # bin data based on maf
        if np.min(data['maf']) == 0:
            log_bins = np.insert(np.logspace(np.log10(np.max(data['maf']) * 0.05), np.log10(np.max(data['maf'])), 20), 0, 0)
        else:
            log_bins = np.logspace(np.log10(np.min(data['maf'])), np.log10(np.max(data['maf'])), 21)
        bin_labels = [f'{str(np.round(start, 3))}-{str(np.round(end, 3))}' for start, end in zip(log_bins[:-1], log_bins[1:])]

        data['bins'] = pd.cut(data.maf, log_bins, labels=False)

        # only if INFO metric is present, make QC plots containing INFO metric
        if not self.df['imp_info'].isnull().all():
            ax1 = data.plot.scatter(x='imp_info', y='log_pval', xlabel="INFO", ylabel="-log$_{\mathrm{10}}$ p-value")
            ax1.set_title(study_id)
            ax2_1 = data.plot.scatter(x='maf', y='imp_info', xlabel="MAF", ylabel="INFO")
            ax2_1.set_title(study_id)
            ax2_2 = data.boxplot(column=['imp_info'], by='bins', showfliers=False, figsize=(8,8))
            ax2_2.set_xticklabels(bin_labels, rotation=90)
            ax2_2.set_xlabel('MAF')
            ax2_2.set_ylabel('INFO')
            ax2_2.set_title(study_id)
            ax2_2.get_figure().suptitle('')

            # save figures
            ax1.get_figure().savefig(f"{out_dir}/{study_id}/scatterPlot-info-logp.png")
            ax2_1.get_figure().savefig(f"{out_dir}/{study_id}/scatterPlot-maf-info.png")
            ax2_2.get_figure().savefig(f"{out_dir}/{study_id}/boxPlot-maf-info.png")

        ax3_1 = data.plot.scatter(x='maf', y='log_pval', xlabel="MAF", ylabel="-log$_{\mathrm{10}}$ p-value")
        ax3_1.set_title(study_id)
        ax3_2 = data.boxplot(column=['log_pval'], by='bins', showfliers=False, figsize=(8,8))
        ax3_2.set_xticks(range(1, 21))
        ax3_2.set_xticklabels(bin_labels, rotation=90)
        ax3_2.set_xlabel('MAF')
        ax3_2.set_ylabel('-log$_{\mathrm{10}}$ p-value')
        ax3_2.set_title(study_id)
        ax3_2.get_figure().suptitle('')
        ax4_1 = data.plot.scatter(x='maf', y='abs_es', xlabel="MAF", ylabel="Absolute beta")
        ax4_1.set_title(study_id)
        ax4_2 = data.boxplot(column=['abs_es'], by='bins', showfliers=False, figsize=(8,8))
        ax4_2.set_xticks(range(1, 21))
        ax4_2.set_xticklabels(bin_labels, rotation=90)
        ax4_2.set_xlabel('MAF')
        ax4_2.set_ylabel('Absolute beta')
        ax4_2.set_title(study_id)
        ax4_2.get_figure().suptitle('')
        ax5 = data.plot.scatter(x='abs_es', y='log_pval', xlabel="Absolute beta", ylabel="-log$_{\mathrm{10}}$ p-value")
        ax5.set_title(study_id)

        # save figures
        ax3_1.get_figure().savefig(f"{out_dir}/{study_id}/scatterPlot-maf-p.png")
        ax3_2.get_figure().savefig(f"{out_dir}/{study_id}/boxPlot-maf-p.png")
        ax4_1.get_figure().savefig(f"{out_dir}/{study_id}/scatterPlot-maf-beta.png")
        ax4_2.get_figure().savefig(f"{out_dir}/{study_id}/boxPlot-maf-beta.png")
        ax5.get_figure().savefig(f"{out_dir}/{study_id}/scatterPlot-beta-logp.png")

    def __generate_qq_plot__(self, out_dir, study_id):
        if not os.path.exists(f"{out_dir}/{study_id}"):
            os.makedirs(f"{out_dir}/{study_id}")

        obs_p = -np.log10(self.df['pval']).sort_values(ascending=False)
        exp_p = -np.log10(np.divide(range(len(self.df['pval'])+1, 1, -1), len(self.df['pval'])))
        df = pd.DataFrame({'exp': exp_p, 'obs': obs_p}, columns=['exp', 'obs'])

        # plot data
        ax = df.plot.scatter(x='exp', y='obs', xlabel="Expected (-log$_{\mathrm{10}}$ p-value)", ylabel="Observed (-log$_{\mathrm{10}}$ p-value)")
        
        # plot confidence interval
        q_pos = np.concatenate([np.arange(99.)/len(obs_p), np.logspace(-np.log10(len(obs_p))+2, 0, 100)])
        q_err = np.zeros([len(q_pos),2])
        for i in range(0, len(q_pos)):
            q_err[i, :] = beta.interval(0.95, len(obs_p)*q_pos[i], len(obs_p) - len(obs_p)*q_pos[i])
        q_err[i, q_err[i, :] < 0] = 1e-15
        q_err[0,:] = q_err[1,:]
        q_err[len(q_err)-1,:] = q_err[len(q_err)-2,:] 
        ax.fill_between(-np.log10(q_pos), -np.log10(q_err[:,0]), -np.log10(q_err[:,1]), color='b', alpha=0.2)

        # set limits
        obs_p[~np.isfinite(obs_p)] = 0
        ax.set_ylim(0, np.ceil(max(obs_p)))
        ax.set_xlim(0, np.ceil(max(exp_p)))
        ax.plot([0, 100], [0, 100],'--k')
        ax.set_title(study_id)
        # save figure
        ax.get_figure().savefig(f"{out_dir}/{study_id}/QQ.png")

    def __generate_manhattan_plot__(self, out_dir, study_id):
        if not os.path.exists(f"{out_dir}/{study_id}"):
            os.makedirs(f"{out_dir}/{study_id}")

        data = self.df[['pval', 'chrom', 'pos']].copy()
        data['log_pval'] = -np.log10(self.df['pval'])

        colors = cycle(["#3489ca", "#53BBD5"])
        last_xpos = 0
        x_ticks, x, y, c = [], [], [], []

        data = data[pd.to_numeric(data['chrom'], errors='coerce').notnull()]
        data.chrom = data.chrom.astype(int)
        
        for chrom, group_data in data.groupby(by='chrom', sort=True):
            color = next(colors)
            for pos, pval in zip(group_data['pos'], group_data['log_pval']):
                x.append(last_xpos + pos)
                y.append(pval)
                c.append(color)

            x_ticks.append([chrom, last_xpos + (group_data['pos'].iloc[0] + group_data['pos'].iloc[-1]) / 2])
            last_xpos = x[-1]
        
        df = pd.DataFrame({'x': x, 'y': y}, columns=['x', 'y'])

        ax = df.plot.scatter(
            x='x',
            y='y',
            c=c,
            alpha=0.8,
            edgecolors="none",
            xlabel="Chromosome",
            ylabel="-log$_{\mathrm{10}}$ p-value",
            figsize=(12, 5),
            s=8
            )

        ax.axhline(y=-np.log10(5e-8), color="r")

        ax.set_xticks([x for _, x in x_ticks])
        ax.set_xticklabels([x for x, _ in x_ticks], rotation=45)

        ax.set_xlim(0, x[-1])

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_title(study_id)

        ax.get_figure().savefig(f"{out_dir}/{study_id}/Manhattan.png")

    def __descriptive_statistics__(self, file_name:str, n_total, n_cases):
        parameters = ["beta", "se", "pval", "maf", "imp_info"] if not self.df['imp_info'].isnull().all() else ["beta", "se", "pval", "maf"]
        with open(file_name, 'w') as file:
            for chr in range(23):
                if chr == 0:
                    file.write(f"Chrom: all\n")
                    file.write(f"\tN variants: {len(self.df)}\n")
                    if not None in (n_total, n_cases):
                        file.write(f"\tlambda: {np.round(self.__calculate_lambda__(n_total, n_cases)[0], 2)}\n")
                        file.write(f"\tlambda_1000: {np.round(self.__calculate_lambda__(n_total, n_cases)[1], 2)}\n")
                    else:
                        file.write(f"\tlambda: NA\n")
                        file.write(f"\tlambda_1000: NA\n")
                    for param in parameters:
                        file.write(f"\t{param}:\n")
                        file.write(f"\t\tMin: {np.round(self.df[param].min(), 2)}\n")
                        file.write(f"\t\t25% Quantile: {np.round(self.df[param].quantile(0.25), 2)}\n")
                        file.write(f"\t\tMedian: {np.round(self.df[param].median(), 2)}\n")
                        file.write(f"\t\tMean: {np.round(self.df[param].mean(), 2)}\n")
                        file.write(f"\t\t75% Quantile: {np.round(self.df[param].quantile(0.75), 2)}\n")
                        file.write(f"\t\tMax: {np.round(self.df[param].max(), 2)}\n")
                else:
                    file.write(f"Chrom: {chr}\n")
                    file.write(f"\tN variants: {len(self.df.loc[self.df['chrom'] == str(chr)])}\n")
                    for param in parameters:
                        file.write(f"\t{param}:\n")
                        file.write(f"\t\tMin: {np.round(self.df.loc[self.df['chrom'] == str(chr)][param].min(), 2)}\n")
                        file.write(f"\t\t25% Quantile: {np.round(self.df.loc[self.df['chrom'] == str(chr)][param].quantile(0.25), 2)}\n")
                        file.write(f"\t\tMedian: {np.round(self.df.loc[self.df['chrom'] == str(chr)][param].median(), 2)}\n")
                        file.write(f"\t\tMean: {np.round(self.df.loc[self.df['chrom'] == str(chr)][param].mean(), 2)}\n")
                        file.write(f"\t\t75% Quantile: {np.round(self.df.loc[self.df['chrom'] == str(chr)][param].quantile(0.75), 2)}\n")
                        file.write(f"\t\tMax: {np.round(self.df.loc[self.df['chrom'] == str(chr)][param].max(), 2)}\n")

    def __descriptive_statistics_table__(self, file_name:str):
        parameters = ["beta", "se", "pval", "maf", "imp_info"] if not self.df['imp_info'].isnull().all() else ["beta", "se", "pval", "maf"]
        desc_all = self.df[['chrom']+parameters].describe().T
        desc_all['chrom'] = 'all'
        desc_all = desc_all.reset_index().set_index(['chrom','index'])

        dfg = self.df[['chrom']+parameters].set_index('chrom').groupby('chrom').apply(lambda d: d.describe().T)
        desc_all = pd.concat([desc_all,dfg])
        desc_all['count'] = desc_all['count'].astype(int)
        desc_all.to_csv(file_name, sep='\t', float_format='{:.2f}'.format)

    def __calculate_lambda__(self, n_total, n_cases):
        pm = np.median(self.df['pval'])
        Chi = chi2.ppf(1 - pm, 1)
        inflation_factor = Chi / chi2.ppf(0.5, 1)

        controls = n_total - n_cases

        l1000 = 1 + (inflation_factor-1)*(1/n_cases+1/controls)/(1/1000+1/1000) if controls != 0 else 1 + (inflation_factor-1)*(1/n_cases)/(1/1000)

        return inflation_factor, l1000

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--out_dir', metavar="<directory>",
                        help="Location of qc plots/files", type=str, required=True)
    parser.add_argument('--sumstat', metavar="<tsv/csv/vcf file>",
                        help="Location of the sumstats file", type=str, required=True)
    parser.add_argument('--study_identifier', metavar="identifier <String>", default='CUSTOM_GWAS',
                        help="GWAS study identifier prefix", type=str)
    parser.add_argument('--n-total', metavar="N <Int>", default=None,
                        help="N total", type=int)
    parser.add_argument('--n-cases', metavar="N <Int>", default=None,
                        help="N cases", type=int)
    parser.add_argument('--c-chrom', type=str, default="chrom")
    parser.add_argument('--c-pos', type=str, default="pos")
    parser.add_argument('--c-info', type=str, default="imp_info")
    parser.add_argument('--c-maf', type=str, default="maf")
    parser.add_argument('--c-beta', type=str, default="beta")
    parser.add_argument('--c-se', type=str, default="se")
    parser.add_argument('--c-pval', type=str, default="pval")
    
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()
    study_id = args.study_identifier
    map_names = {
        args.c_chrom: "chrom",
        args.c_pos: "pos",
        args.c_info: "imp_info",
        args.c_maf: "maf",
        args.c_beta: "beta",
        args.c_se: "se",
        args.c_pval: "pval"
    }

    print(f"Working on: {os.path.basename(args.sumstat)} ({study_id})")

    sumstats_file = SumstatsFile(args.sumstat, map_names)

    sumstats_file.format_gwas(study_id)
    if sumstats_file.df.chrom.values[0].startswith('chr'):
        sumstats_file.df['chrom'] = sumstats_file.df['chrom'].str.replace('chr', '')

    sumstats_file.quality_control(args.out_dir, study_id, args.n_total, args.n_cases)
