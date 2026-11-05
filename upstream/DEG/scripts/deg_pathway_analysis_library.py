#!/usr/bin/env python3
# Paths below refer to the upstream processing pipeline that produced the
# deposited intermediates; it is not part of this release's run path. This
# file is included as a record of how the input was produced; it is not
# called by _run_all_panels.sh or ./run.
"""
Unified Differential Gene Expression and Pathway Enrichment Analysis Library
"""

import scanpy as sc
import pandas as pd
import numpy as np
import yaml
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from abc import ABC, abstractmethod
import warnings
warnings.filterwarnings('ignore')

class DEGConfig:

    def __init__(self, config_path: str):
        self.config_path = config_path
        self.config = self._load_config()
        self._validate_config()

    def _load_config(self) -> Dict:
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)

    def _validate_config(self):
        required_sections = ['data', 'analysis']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required config section: {section}")

        if 'adata_path' not in self.config['data']:
            raise ValueError("Missing required config: data.adata_path")

        if 'pathway_enrichment' not in self.config:
            self.config['pathway_enrichment'] = {'enabled': False}

        if 'visualization' not in self.config:
            self.config['visualization'] = {
                'generate_plots': True,
                'plot_format': 'png'
            }

    def get(self, *keys, default=None):
        value = self.config
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
        return value

    def get_comparison_configs(self) -> Dict:
        return {
            '01_treatment_comparison': {
                'column': 'stomach_treatment',
                'group1': 'Pre',
                'group2': 'Post',
                'description': 'Stomach pre vs post treatment comparison'
            },
            '02_stomach_pre_response': {
                'column': 'stomach_pre_grouping',
                'group1': 'No-response',
                'group2': 'Responsed',
                'description': 'Pre-treatment stomach response comparison'
            },
            '03_stomach_post_response': {
                'column': 'stomach_post_grouping',
                'group1': 'No-response',
                'group2': 'Responsed',
                'description': 'Post-treatment stomach response comparison'
            },
            '04_liver_pre_response': {
                'column': 'liver_pre_grouping',
                'group1': 'No-response',
                'group2': 'Responsed',
                'description': 'Pre-treatment liver response comparison'
            },
            '05_pre_stomach_vs_liver': {
                'column': 'pre_stomach_vs_pre_liver',
                'group1': 'Stomach',
                'group2': 'Liver',
                'description': 'Pre-treatment stomach vs liver comparison'
            },
            '06_overall_responders_vs_non_responders': {
                'column': 'combined_stomach_responder_status',
                'group1': 'No-response',
                'group2': 'Responsed',
                'description': 'Overall stomach responders vs non-responders'
            }
        }


class DataLoader:

    def __init__(self, config: DEGConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.adata = None

    def load_h5ad_data(self) -> sc.AnnData:
        adata_path = self.config.get('data', 'adata_path')

        self.adata = sc.read_h5ad(adata_path)

        return self.adata

    def validate_data(self) -> bool:
        if self.adata is None:
            self.logger.error("No data loaded")
            return False

        required_cols = ['major_cell_type', 'minor_cell_state']
        missing_cols = [col for col in required_cols if col not in self.adata.obs.columns]

        if missing_cols:
            self.logger.warning(f"Missing columns: {missing_cols}")

        return True

    def prepare_for_deg(self, preserve_raw: bool = True):
        if preserve_raw and self.adata.raw is None:
            self.adata.raw = self.adata.copy()

        # Apply log1p if not already done
        if 'log1p' not in self.adata.uns:
            sc.pp.log1p(self.adata)

    def filter_by_comparison(self, comparison_config: Dict) -> sc.AnnData:
        group_col = comparison_config['column']
        valid_groups = [comparison_config['group1'], comparison_config['group2']]

        if group_col not in self.adata.obs.columns:
            self.logger.error(f"Column {group_col} not found in data")
            return sc.AnnData()

        mask = self.adata.obs[group_col].isin(valid_groups)
        adata_filtered = self.adata[mask].copy()

        return adata_filtered


class DEGMethodStrategy(ABC):

    @abstractmethod
    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        pass

    @staticmethod
    @abstractmethod
    def check_dependencies() -> bool:
        pass

    @staticmethod
    @abstractmethod
    def get_description() -> str:
        pass






class MASTStrategy(DEGMethodStrategy):

    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        try:
            import rpy2.robjects as robjects
            from rpy2.robjects import pandas2ri, numpy2ri
            from rpy2.robjects.conversion import localconverter
            from rpy2.robjects.packages import importr
            import anndata2ri

            r = robjects.r
            mast = importr('MAST')

            group_col = comparison_config['column']
            reference_group = comparison_config['group1']
            comparison_group = comparison_config['group2']

            groups_to_include = [reference_group, comparison_group]
            mask = adata.obs[group_col].isin(groups_to_include)
            adata_subset = adata[mask].copy()

            if adata_subset.n_obs == 0:
                return pd.DataFrame()

            # Apply 8000 cell cap if needed (for memory efficiency)
            max_cells = kwargs.get('max_cells_mast', 8000)
            if adata_subset.n_obs > max_cells:
                np.random.seed(42)  # For reproducible results
                selected_indices = np.random.choice(adata_subset.n_obs, max_cells, replace=False)
                adata_subset = adata_subset[selected_indices, :].copy()

            min_pct_expressed = kwargs.get('min_pct_expressed', 0.1)
            if hasattr(adata_subset.X, 'toarray'):
                gene_filter = (adata_subset.X.toarray() > 0).mean(axis=0) >= min_pct_expressed
            else:
                gene_filter = (adata_subset.X > 0).mean(axis=0) >= min_pct_expressed
            adata_filtered = adata_subset[:, gene_filter].copy()

            if hasattr(adata_filtered.X, 'toarray'):
                expr_matrix = adata_filtered.X.toarray()
            else:
                expr_matrix = adata_filtered.X.copy()

            cell_metadata = adata_filtered.obs[group_col].values
            cell_names = adata_filtered.obs.index.values
            gene_names = adata_filtered.var.index.values

            sample_id = None
            for col in ['Sample ID', 'sample_id', 'Patient ID', 'patient_id']:
                if col in adata_filtered.obs.columns:
                    sample_id = adata_filtered.obs[col].values
                    break

            if sample_id is None:
                sample_id = np.array([f"sample_{i % 10}" for i in range(len(cell_names))])

            with localconverter(robjects.default_converter + pandas2ri.converter + numpy2ri.converter):
                r.assign('expr_matrix', expr_matrix)
                r.assign('cell_metadata', cell_metadata)
                r.assign('cell_names', cell_names)
                r.assign('gene_names', gene_names)
                r.assign('reference_group', reference_group)
                r.assign('comparison_group', comparison_group)
                r.assign('sample_id', sample_id)

            r('''
            library(MAST)
            library(SingleCellExperiment)

            # Data is already log-transformed from Python
            # Replace NA/NaN/Inf with 0 (unexpressed)
            expr_mat <- expr_matrix
            expr_mat[is.na(expr_mat) | is.nan(expr_mat) | is.infinite(expr_mat)] <- 0

            expr_mat <- t(expr_mat)
            rownames(expr_mat) <- gene_names
            colnames(expr_mat) <- cell_names

            cell_meta <- data.frame(
                condition = cell_metadata,
                sample_id = as.factor(sample_id),
                wellKey = cell_names,
                stringsAsFactors = FALSE
            )
            rownames(cell_meta) <- cell_names

            cdr_values <- colSums(expr_mat > 0)
            cell_meta$cngeneson <- as.numeric(scale(cdr_values))

            gene_meta <- data.frame(
                primerid = gene_names,
                stringsAsFactors = FALSE
            )
            rownames(gene_meta) <- gene_names

            sca <- FromMatrix(
                exprsArray = expr_mat,
                cData = cell_meta,
                fData = gene_meta
            )

            condition_factor <- factor(sca$condition)
            if (reference_group %in% levels(condition_factor)) {
                sca$condition <- relevel(condition_factor, ref = reference_group)
            } else {
                sca$condition <- condition_factor
            }

            zlmCond <- zlm(~ condition + sample_id + cngeneson, sca)

            contrast_name <- paste0("condition", comparison_group)
            summaryCond <- summary(zlmCond, doLRT = contrast_name)
            summaryDt <- summaryCond$datatable

            fcHurdle <- merge(
                summaryDt[contrast==contrast_name & component=="H", .(primerid, `Pr(>Chisq)`)],
                summaryDt[contrast==contrast_name & component=="logFC", .(primerid, coef, ci.hi, ci.lo)],
                by="primerid"
            )

            fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`, method="fdr")

            results <- data.frame(
                gene = fcHurdle$primerid,
                logfoldchanges = fcHurdle$coef,
                pvals = fcHurdle$`Pr(>Chisq)`,
                pvals_adj = fcHurdle$fdr,
                scores = -log10(fcHurdle$`Pr(>Chisq)`),
                stringsAsFactors = FALSE
            )

            results <- results[!is.na(results$pvals) & !is.na(results$logfoldchanges), ]
            results$abs_logfc <- abs(results$logfoldchanges)
            results <- results[order(-results$abs_logfc), ]
            results$abs_logfc <- NULL
            ''')

            with localconverter(robjects.default_converter + pandas2ri.converter + numpy2ri.converter):
                results_r = robjects.conversion.rpy2py(r['results'])

            if not isinstance(results_r, pd.DataFrame):
                results_df = pd.DataFrame(results_r)
            else:
                results_df = results_r

            results_df.set_index('gene', inplace=True)
            results_df = results_df.loc[:, ~results_df.columns.duplicated()]

            return results_df

        except ImportError as e:
            raise ImportError(f"MAST requires rpy2 and R MAST package\nError: {str(e)}")
        except Exception as e:
            raise RuntimeError(f"MAST analysis failed: {str(e)}")

    @staticmethod
    def check_dependencies() -> bool:
        try:
            import rpy2.robjects
            from rpy2.robjects.packages import importr
            import anndata2ri
            mast = importr('MAST')
            return True
        except:
            return False

    @staticmethod
    def get_description() -> str:
        return "MAST hurdle model (requires: rpy2, anndata2ri, R MAST package)"






