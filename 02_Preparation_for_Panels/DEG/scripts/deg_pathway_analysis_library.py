#!/usr/bin/env python3
"""
Unified Differential Gene Expression and Pathway Enrichment Analysis Library

This library provides a complete pipeline for:
1. DEG (Differential Gene Expression) analysis with multiple statistical methods
2. Pathway enrichment analysis (Enrichr + GSEA)
3. Visualization and reporting

Author: Project 4 - Stomach Cancer Single Cell Analysis
Version: 1.0
Date: 2025
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

# ============================================================================
# CONFIGURATION CLASS
# ============================================================================

class DEGConfig:
    """Configuration management for DEG and pathway enrichment analysis"""

    def __init__(self, config_path: str):
        """Initialize configuration from YAML file"""
        self.config_path = config_path
        self.config = self._load_config()
        self._validate_config()

    def _load_config(self) -> Dict:
        """Load configuration from YAML file"""
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)

    def _validate_config(self):
        """Validate required configuration sections"""
        required_sections = ['data', 'analysis']
        for section in required_sections:
            if section not in self.config:
                raise ValueError(f"Missing required config section: {section}")

        # Validate data paths
        if 'adata_path' not in self.config['data']:
            raise ValueError("Missing required config: data.adata_path")

        # Set defaults for optional sections
        if 'pathway_enrichment' not in self.config:
            self.config['pathway_enrichment'] = {'enabled': False}

        if 'visualization' not in self.config:
            self.config['visualization'] = {
                'generate_plots': True,
                'plot_format': 'png'
            }

    def get(self, *keys, default=None):
        """Get nested configuration value"""
        value = self.config
        for key in keys:
            if isinstance(value, dict) and key in value:
                value = value[key]
            else:
                return default
        return value

    def get_comparison_configs(self) -> Dict:
        """Get predefined comparison configurations"""
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


# ============================================================================
# DATA LOADER CLASS
# ============================================================================

class DataLoader:
    """Handles loading and preprocessing of single-cell data"""

    def __init__(self, config: DEGConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.adata = None

    def load_h5ad_data(self) -> sc.AnnData:
        """Load AnnData object from h5ad file"""
        adata_path = self.config.get('data', 'adata_path')
        self.logger.info(f"Loading data from {adata_path}")

        self.adata = sc.read_h5ad(adata_path)
        self.logger.info(f"Loaded {self.adata.n_obs} cells and {self.adata.n_vars} genes")

        return self.adata

    def validate_data(self) -> bool:
        """Validate that data has required annotations"""
        if self.adata is None:
            self.logger.error("No data loaded")
            return False

        # Check for required columns
        required_cols = ['major_cell_type', 'minor_cell_state']
        missing_cols = [col for col in required_cols if col not in self.adata.obs.columns]

        if missing_cols:
            self.logger.warning(f"Missing columns: {missing_cols}")

        return True

    def prepare_for_deg(self, preserve_raw: bool = True):
        """Prepare data for DEG analysis"""
        if preserve_raw and self.adata.raw is None:
            self.logger.info("Preserving raw count data")
            self.adata.raw = self.adata.copy()

        # Apply log1p if not already done
        if 'log1p' not in self.adata.uns:
            self.logger.info("Applying log1p transformation")
            sc.pp.log1p(self.adata)

    def filter_by_comparison(self, comparison_config: Dict) -> sc.AnnData:
        """Filter data for specific comparison groups"""
        group_col = comparison_config['column']
        valid_groups = [comparison_config['group1'], comparison_config['group2']]

        if group_col not in self.adata.obs.columns:
            self.logger.error(f"Column {group_col} not found in data")
            return sc.AnnData()

        mask = self.adata.obs[group_col].isin(valid_groups)
        adata_filtered = self.adata[mask].copy()

        self.logger.info(f"Filtered to {adata_filtered.n_obs} cells for {group_col}")
        return adata_filtered


# ============================================================================
# DEG METHOD STRATEGY CLASSES
# ============================================================================

class DEGMethodStrategy(ABC):
    """Abstract base class for DEG analysis methods"""

    @abstractmethod
    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        """Run DEG analysis and return results DataFrame"""
        pass

    @staticmethod
    @abstractmethod
    def check_dependencies() -> bool:
        """Check if required dependencies are available"""
        pass

    @staticmethod
    @abstractmethod
    def get_description() -> str:
        """Get method description"""
        pass


class WilcoxonStrategy(DEGMethodStrategy):
    """Wilcoxon rank-sum test for DEG analysis"""

    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        """Run Wilcoxon test using Scanpy"""
        try:
            group_col = comparison_config['column']
            reference_group = comparison_config['group1']
            comparison_group = comparison_config['group2']

            # Filter to only include the two groups
            groups_to_include = [reference_group, comparison_group]
            mask = adata.obs[group_col].isin(groups_to_include)
            adata_subset = adata[mask].copy()

            if adata_subset.n_obs == 0:
                return pd.DataFrame()

            # Check if both groups are present
            group_counts = adata_subset.obs[group_col].value_counts()
            if len(group_counts) < 2:
                return pd.DataFrame()

            # Run Wilcoxon test
            sc.tl.rank_genes_groups(
                adata_subset,
                groupby=group_col,
                groups=[comparison_group],
                reference=reference_group,
                method='wilcoxon',
                use_raw=kwargs.get('use_raw', False),
                layer=kwargs.get('layer', None),
                n_genes=adata_subset.n_vars,
                pts=True
            )

            # Extract results
            result_keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
            results_dict = {}

            for key in result_keys:
                if key in adata_subset.uns['rank_genes_groups']:
                    results_dict[key] = adata_subset.uns['rank_genes_groups'][key][comparison_group]

            if 'names' not in results_dict or len(results_dict['names']) == 0:
                return pd.DataFrame()

            # Create results DataFrame
            deg_results = pd.DataFrame(index=results_dict['names'])
            for key in ['logfoldchanges', 'pvals', 'pvals_adj', 'scores']:
                if key in results_dict:
                    deg_results[key] = results_dict[key]

            # Add percentage information
            if 'pts' in adata_subset.uns['rank_genes_groups']:
                pts_data = adata_subset.uns['rank_genes_groups']['pts']
                if comparison_group in pts_data:
                    deg_results['pct_1'] = pts_data[comparison_group]
                if reference_group in pts_data:
                    deg_results['pct_2'] = pts_data[reference_group]

            # Remove NaN rows
            deg_results = deg_results.dropna()

            # Sort by absolute log fold change
            if 'logfoldchanges' in deg_results.columns:
                deg_results['abs_logfc'] = np.abs(deg_results['logfoldchanges'])
                deg_results = deg_results.sort_values('abs_logfc', ascending=False)
                deg_results = deg_results.drop('abs_logfc', axis=1)

            return deg_results

        except Exception as e:
            print(f"Error in Wilcoxon test: {str(e)}")
            return pd.DataFrame()

    @staticmethod
    def check_dependencies() -> bool:
        """Wilcoxon is always available with scanpy"""
        return True

    @staticmethod
    def get_description() -> str:
        return "Wilcoxon rank-sum test (non-parametric, always available)"


class DESeq2Strategy(DEGMethodStrategy):
    """DESeq2-style analysis using PyDESeq2 with pseudobulk"""

    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        """Run DESeq2 analysis with pseudobulk aggregation"""
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats

            group_col = comparison_config['column']
            reference_group = comparison_config['group1']
            comparison_group = comparison_config['group2']

            # Filter to groups
            groups_to_include = [reference_group, comparison_group]
            mask = adata.obs[group_col].isin(groups_to_include)
            adata_subset = adata[mask].copy()

            if adata_subset.n_obs == 0:
                return pd.DataFrame()

            # Determine sample column
            sample_col = None
            for col in ['Sample ID', 'sample_id', 'Patient ID']:
                if col in adata_subset.obs.columns:
                    sample_col = col
                    break

            if sample_col is None:
                # Create artificial samples
                n_cells_per_sample = max(100, adata_subset.n_obs // 20)
                adata_subset.obs['artificial_sample'] = (
                    adata_subset.obs.index.to_series().reset_index(drop=True) // n_cells_per_sample
                ).astype(str)
                sample_col = 'artificial_sample'

            # Get count data
            if hasattr(adata_subset, 'raw') and adata_subset.raw is not None:
                count_data = adata_subset.raw.X
                gene_names = adata_subset.raw.var.index
            else:
                count_data = adata_subset.X
                gene_names = adata_subset.var.index

            # Create pseudobulk
            pseudobulk_data = []
            pseudobulk_metadata = []

            for sample_id in adata_subset.obs[sample_col].unique():
                sample_mask = adata_subset.obs[sample_col] == sample_id
                sample_indices = np.where(sample_mask)[0]

                # Check for mixed conditions
                sample_conditions = adata_subset.obs[group_col].iloc[sample_indices].unique()
                if len(sample_conditions) > 1:
                    continue

                condition = sample_conditions[0]

                # Aggregate counts
                if hasattr(count_data, 'toarray'):
                    pseudobulk_counts = np.array(count_data[sample_indices, :].toarray().sum(axis=0)).flatten()
                else:
                    pseudobulk_counts = np.array(count_data[sample_indices, :].sum(axis=0)).flatten()

                pseudobulk_data.append(pseudobulk_counts)
                pseudobulk_metadata.append({
                    'sample_id': f"{sample_id}_{condition}",
                    'condition': condition
                })

            if len(pseudobulk_data) < 4:
                print(f"Insufficient samples for DESeq2 (need ≥4, found {len(pseudobulk_data)})")
                print("Falling back to Wilcoxon")
                return WilcoxonStrategy().run(adata, comparison_config, **kwargs)

            # Create count matrix
            count_matrix = pd.DataFrame(
                np.array(pseudobulk_data),
                columns=gene_names,
                index=[meta['sample_id'] for meta in pseudobulk_metadata]
            )

            metadata_df = pd.DataFrame(pseudobulk_metadata)
            metadata_df.index = metadata_df['sample_id']

            # Filter low-expression genes
            min_count_threshold = kwargs.get('min_counts', 10)
            gene_filter = (count_matrix > min_count_threshold).sum(axis=0) >= 2
            count_matrix = count_matrix.loc[:, gene_filter]

            # Run DESeq2
            dds = DeseqDataSet(
                adata=None,
                counts=count_matrix,
                metadata=metadata_df,
                design_factors=['condition']
            )
            dds.deseq2()

            stat_res = DeseqStats(
                dds,
                contrast=['condition', comparison_group, reference_group]
            )
            stat_res.summary()

            # Extract and format results
            results_df = stat_res.results_df.copy()
            results_df = results_df.rename(columns={
                'log2FoldChange': 'logfoldchanges',
                'pvalue': 'pvals',
                'padj': 'pvals_adj',
                'stat': 'scores'
            })

            # Sort and clean
            results_df = results_df.dropna(subset=['logfoldchanges', 'pvals'])
            if 'logfoldchanges' in results_df.columns:
                results_df['abs_logfc'] = np.abs(results_df['logfoldchanges'])
                results_df = results_df.sort_values('abs_logfc', ascending=False)
                results_df = results_df.drop('abs_logfc', axis=1)

            return results_df

        except ImportError as e:
            raise ImportError(f"DESeq2 requires PyDESeq2: pip install pydeseq2\nError: {str(e)}")
        except Exception as e:
            raise RuntimeError(f"DESeq2 analysis failed: {str(e)}")

    @staticmethod
    def check_dependencies() -> bool:
        """Check if PyDESeq2 is available"""
        try:
            from pydeseq2.dds import DeseqDataSet
            from pydeseq2.ds import DeseqStats
            return True
        except ImportError:
            return False

    @staticmethod
    def get_description() -> str:
        return "DESeq2-style analysis with pseudobulk (requires: pydeseq2)"


class MASTStrategy(DEGMethodStrategy):
    """MAST hurdle model analysis via R/rpy2"""

    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        """Run MAST analysis"""
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

            # Filter to groups
            groups_to_include = [reference_group, comparison_group]
            mask = adata.obs[group_col].isin(groups_to_include)
            adata_subset = adata[mask].copy()

            if adata_subset.n_obs == 0:
                return pd.DataFrame()

            # Apply 8000 cell cap if needed (for memory efficiency)
            max_cells = kwargs.get('max_cells_mast', 8000)
            if adata_subset.n_obs > max_cells:
                print(f"Applying cell cap: downsampling from {adata_subset.n_obs} to {max_cells} cells")
                np.random.seed(42)  # For reproducible results
                selected_indices = np.random.choice(adata_subset.n_obs, max_cells, replace=False)
                adata_subset = adata_subset[selected_indices, :].copy()
                print(f"After downsampling: {adata_subset.n_obs} cells")

            # Filter genes with low expression
            min_pct_expressed = kwargs.get('min_pct_expressed', 0.1)
            if hasattr(adata_subset.X, 'toarray'):
                gene_filter = (adata_subset.X.toarray() > 0).mean(axis=0) >= min_pct_expressed
            else:
                gene_filter = (adata_subset.X > 0).mean(axis=0) >= min_pct_expressed
            adata_filtered = adata_subset[:, gene_filter].copy()

            # Extract data for R
            if hasattr(adata_filtered.X, 'toarray'):
                expr_matrix = adata_filtered.X.toarray()
            else:
                expr_matrix = adata_filtered.X.copy()

            cell_metadata = adata_filtered.obs[group_col].values
            cell_names = adata_filtered.obs.index.values
            gene_names = adata_filtered.var.index.values

            # Extract sample ID for patient effects
            sample_id = None
            for col in ['Sample ID', 'sample_id', 'Patient ID', 'patient_id']:
                if col in adata_filtered.obs.columns:
                    sample_id = adata_filtered.obs[col].values
                    print(f"Using '{col}' for patient effects in MAST model")
                    break

            if sample_id is None:
                print("Warning: No Sample ID or Patient ID found. Model will not include patient effects.")
                # Create artificial sample IDs as fallback
                sample_id = np.array([f"sample_{i % 10}" for i in range(len(cell_names))])

            # Pass to R using context manager for conversion
            with localconverter(robjects.default_converter + pandas2ri.converter + numpy2ri.converter):
                r.assign('expr_matrix', expr_matrix)
                r.assign('cell_metadata', cell_metadata)
                r.assign('cell_names', cell_names)
                r.assign('gene_names', gene_names)
                r.assign('reference_group', reference_group)
                r.assign('comparison_group', comparison_group)
                r.assign('sample_id', sample_id)

            # Run MAST in R
            r('''
            library(MAST)
            library(SingleCellExperiment)

            # Data is already log-transformed from Python
            # Replace NA/NaN/Inf with 0 (unexpressed)
            expr_mat <- expr_matrix
            expr_mat[is.na(expr_mat) | is.nan(expr_mat) | is.infinite(expr_mat)] <- 0

            # Transpose for MAST
            expr_mat <- t(expr_mat)
            rownames(expr_mat) <- gene_names
            colnames(expr_mat) <- cell_names

            # Create metadata with sample_id for patient effects
            cell_meta <- data.frame(
                condition = cell_metadata,
                sample_id = as.factor(sample_id),
                wellKey = cell_names,
                stringsAsFactors = FALSE
            )
            rownames(cell_meta) <- cell_names

            # Calculate CDR
            cdr_values <- colSums(expr_mat > 0)
            cell_meta$cngeneson <- as.numeric(scale(cdr_values))

            # Create gene metadata
            gene_meta <- data.frame(
                primerid = gene_names,
                stringsAsFactors = FALSE
            )
            rownames(gene_meta) <- gene_names

            # Create MAST object
            sca <- FromMatrix(
                exprsArray = expr_mat,
                cData = cell_meta,
                fData = gene_meta
            )

            # Set reference level
            condition_factor <- factor(sca$condition)
            if (reference_group %in% levels(condition_factor)) {
                sca$condition <- relevel(condition_factor, ref = reference_group)
            } else {
                sca$condition <- condition_factor
            }

            # Fit model with patient effects (sample_id)
            # Use fixed effects for sample_id to account for patient-level variation
            zlmCond <- zlm(~ condition + sample_id + cngeneson, sca)

            # The coefficient name in the model is "conditionComparisonGroup", not just "ComparisonGroup"
            contrast_name <- paste0("condition", comparison_group)
            summaryCond <- summary(zlmCond, doLRT = contrast_name)
            summaryDt <- summaryCond$datatable

            # Extract results
            fcHurdle <- merge(
                summaryDt[contrast==contrast_name & component=="H", .(primerid, `Pr(>Chisq)`)],
                summaryDt[contrast==contrast_name & component=="logFC", .(primerid, coef, ci.hi, ci.lo)],
                by="primerid"
            )

            fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`, method="fdr")

            # Format results
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

            # Get results from R using context manager for conversion
            with localconverter(robjects.default_converter + pandas2ri.converter + numpy2ri.converter):
                results_r = robjects.conversion.rpy2py(r['results'])

            # Convert to pandas DataFrame (results_r might be a recarray)
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
        """Check if R and MAST are available"""
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


class TTestStrategy(DEGMethodStrategy):
    """Student's t-test for DEG analysis"""

    def run(self, adata: sc.AnnData, comparison_config: Dict, **kwargs) -> pd.DataFrame:
        """Run t-test using Scanpy"""
        try:
            group_col = comparison_config['column']
            reference_group = comparison_config['group1']
            comparison_group = comparison_config['group2']

            groups_to_include = [reference_group, comparison_group]
            mask = adata.obs[group_col].isin(groups_to_include)
            adata_subset = adata[mask].copy()

            if adata_subset.n_obs == 0:
                return pd.DataFrame()

            sc.tl.rank_genes_groups(
                adata_subset,
                groupby=group_col,
                groups=[comparison_group],
                reference=reference_group,
                method='t-test',
                use_raw=kwargs.get('use_raw', False),
                layer=kwargs.get('layer', None),
                n_genes=adata_subset.n_vars
            )

            # Extract results (same as Wilcoxon)
            result_keys = ['names', 'logfoldchanges', 'pvals', 'pvals_adj', 'scores']
            results_dict = {}

            for key in result_keys:
                if key in adata_subset.uns['rank_genes_groups']:
                    results_dict[key] = adata_subset.uns['rank_genes_groups'][key][comparison_group]

            if 'names' not in results_dict:
                return pd.DataFrame()

            deg_results = pd.DataFrame(index=results_dict['names'])
            for key in ['logfoldchanges', 'pvals', 'pvals_adj', 'scores']:
                if key in results_dict:
                    deg_results[key] = results_dict[key]

            deg_results = deg_results.dropna()

            if 'logfoldchanges' in deg_results.columns:
                deg_results['abs_logfc'] = np.abs(deg_results['logfoldchanges'])
                deg_results = deg_results.sort_values('abs_logfc', ascending=False)
                deg_results = deg_results.drop('abs_logfc', axis=1)

            return deg_results

        except Exception as e:
            print(f"Error in t-test: {str(e)}")
            return pd.DataFrame()

    @staticmethod
    def check_dependencies() -> bool:
        """T-test is always available with scanpy"""
        return True

    @staticmethod
    def get_description() -> str:
        return "Student's t-test (parametric, always available)"


# ============================================================================
# DEG ANALYZER CLASS
# ============================================================================

class DEGAnalyzer:
    """Main DEG analysis orchestrator using Strategy Pattern"""

    def __init__(self, config: DEGConfig, data_loader: DataLoader, logger: logging.Logger):
        self.config = config
        self.data_loader = data_loader
        self.logger = logger
        self.results = {}

        # Initialize method strategies
        self.methods = {
            'wilcoxon': WilcoxonStrategy(),
            'deseq2': DESeq2Strategy(),
            'mast': MASTStrategy(),
            'ttest': TTestStrategy()
        }

    def get_method_strategy(self, method_name: str) -> DEGMethodStrategy:
        """Get strategy for specified method"""
        if method_name not in self.methods:
            raise ValueError(f"Unknown method: {method_name}. Available: {list(self.methods.keys())}")
        return self.methods[method_name]

    def check_method_available(self, method_name: str) -> bool:
        """Check if method dependencies are available"""
        strategy = self.get_method_strategy(method_name)
        return strategy.check_dependencies()

    def run_single_comparison(self,
                            comparison_name: str,
                            comparison_config: Dict,
                            cell_type_level: str,
                            method: str,
                            output_dir: Path) -> Dict:
        """Run DEG analysis for a single comparison"""

        self.logger.info(f"Running {comparison_name}: {comparison_config['description']}")
        self.logger.info(f"Method: {method}, Cell type level: {cell_type_level}")

        # Check method availability
        if not self.check_method_available(method):
            self.logger.error(f"Method {method} not available - missing dependencies")
            return {}

        # Filter data for this comparison
        adata_comp = self.data_loader.filter_by_comparison(comparison_config)
        if adata_comp.n_obs == 0:
            self.logger.warning(f"No cells for comparison {comparison_name}")
            return {}

        # Determine cell type column
        cell_type_col = 'major_cell_type' if cell_type_level == 'major' else 'minor_cell_state'

        if cell_type_col not in adata_comp.obs.columns:
            self.logger.error(f"Column {cell_type_col} not found")
            return {}

        # Get unique cell types
        cell_types = adata_comp.obs[cell_type_col].unique()
        self.logger.info(f"Found {len(cell_types)} {cell_type_level} cell types")

        # Run DEG for each cell type
        comparison_results = {}
        min_cells = self.config.get('analysis', 'min_cells_per_group', default=10)

        for cell_type in cell_types:
            self.logger.info(f"  Processing {cell_type}...")

            # Filter to current cell type
            cell_mask = adata_comp.obs[cell_type_col] == cell_type
            adata_celltype = adata_comp[cell_mask].copy()

            if adata_celltype.n_obs < min_cells:
                self.logger.info(f"    Skipping {cell_type}: insufficient cells ({adata_celltype.n_obs})")
                continue

            # Check if both groups present
            group_col = comparison_config['column']
            group_counts = adata_celltype.obs[group_col].value_counts()
            if len(group_counts) < 2:
                self.logger.info(f"    Skipping {cell_type}: only one group present")
                continue

            # Run DEG analysis
            try:
                strategy = self.get_method_strategy(method)
                deg_results = strategy.run(adata_celltype, comparison_config)

                if deg_results is not None and len(deg_results) > 0:
                    # Filter significant results
                    sig_results = self._filter_significant_genes(deg_results)

                    comparison_results[cell_type] = {
                        'all_results': deg_results,
                        'significant': sig_results,
                        'n_cells': adata_celltype.n_obs,
                        'group_counts': group_counts.to_dict()
                    }

                    self.logger.info(f"    Found {len(sig_results)} significant DEGs")

            except Exception as e:
                self.logger.error(f"    Error processing {cell_type}: {str(e)}")
                continue

        # Save results
        if comparison_results:
            self._save_comparison_results(
                comparison_results,
                comparison_name,
                cell_type_level,
                method,
                output_dir
            )

        return comparison_results

    def _filter_significant_genes(self, deg_results: pd.DataFrame) -> pd.DataFrame:
        """Filter significant genes based on thresholds"""
        logfc_thresh = self.config.get('analysis', 'logfc_threshold', default=0.25)
        pval_thresh = self.config.get('analysis', 'pval_threshold', default=0.05)

        mask = (np.abs(deg_results['logfoldchanges']) >= logfc_thresh) & \
               (deg_results['pvals_adj'] < pval_thresh)

        return deg_results[mask].copy()

    def _save_comparison_results(self, results: Dict, comparison_name: str,
                                cell_type_level: str, method: str, output_dir: Path):
        """Save DEG results to files"""
        # Create output directories
        level_dir = output_dir / comparison_name / f"{cell_type_level}_cell_types"
        level_dir.mkdir(parents=True, exist_ok=True)

        # Save individual cell type results
        summary_data = []
        for cell_type, cell_results in results.items():
            # Save detailed results
            filename = f"{cell_type.replace(' ', '_').replace('+', 'pos')}_{method}_deg.csv"
            filepath = level_dir / filename
            cell_results['all_results'].to_csv(filepath, index=True)

            # Collect summary statistics
            summary_data.append({
                'cell_type': cell_type,
                'n_cells': cell_results['n_cells'],
                'n_genes_tested': len(cell_results['all_results']),
                'n_significant': len(cell_results['significant']),
                'n_upregulated': len(cell_results['significant'][
                    cell_results['significant']['logfoldchanges'] > 0
                ]),
                'n_downregulated': len(cell_results['significant'][
                    cell_results['significant']['logfoldchanges'] < 0
                ]),
                'group_counts': str(cell_results['group_counts'])
            })

        # Save summary
        if summary_data:
            summary_df = pd.DataFrame(summary_data)
            summary_path = level_dir / f"deg_results_summary_{method}.csv"
            summary_df.to_csv(summary_path, index=False)
            self.logger.info(f"Saved summary to {summary_path}")

    def run_all_comparisons(self, output_dir: Path) -> Dict:
        """Run DEG analysis for all configured comparisons"""
        all_results = {}

        # Get configurations
        all_comparison_configs = self.config.get_comparison_configs()
        selected_comparisons = self.config.get('analysis', 'comparisons',
                                               default=list(all_comparison_configs.keys()))
        methods = self.config.get('analysis', 'methods', default=['wilcoxon'])
        cell_type_level = self.config.get('analysis', 'cell_type_level', default='major')

        for comp_name in selected_comparisons:
            if comp_name not in all_comparison_configs:
                self.logger.warning(f"Unknown comparison: {comp_name}")
                continue

            comp_config = all_comparison_configs[comp_name]

            for method in methods:
                try:
                    results = self.run_single_comparison(
                        comp_name,
                        comp_config,
                        cell_type_level,
                        method,
                        output_dir
                    )

                    key = f"{comp_name}_{cell_type_level}_{method}"
                    all_results[key] = results

                except Exception as e:
                    self.logger.error(f"Error in {comp_name} {method}: {str(e)}")
                    continue

        return all_results


# ============================================================================
# PATHWAY ENRICHMENT ANALYZER CLASS
# ============================================================================

class PathwayEnrichmentAnalyzer:
    """Pathway enrichment analysis using Enrichr and GSEA"""

    def __init__(self, config: DEGConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

        # Check if gseapy is available
        try:
            import gseapy as gp
            self.gp = gp
            self.available = True
        except ImportError:
            self.logger.warning("gseapy not available - pathway enrichment disabled")
            self.available = False

    def load_deg_results(self, filepath: Path) -> pd.DataFrame:
        """Load DEG results from CSV file"""
        try:
            df = pd.read_csv(filepath, index_col=0)
            self.logger.info(f"Loaded {len(df)} genes from {filepath.name}")
            return df
        except Exception as e:
            self.logger.error(f"Error loading {filepath}: {e}")
            return pd.DataFrame()

    def filter_significant_genes(self, df: pd.DataFrame) -> Tuple[List[str], pd.DataFrame]:
        """Filter significant genes for enrichment analysis"""
        pval_threshold = self.config.get('pathway_enrichment', 'pval_threshold', default=0.05)
        logfc_threshold = self.config.get('pathway_enrichment', 'logfc_threshold', default=0.5)

        # Handle different p-value column names (for MAST compatibility)
        pval_col = None
        for col in ['pvals_adj', 'Pr(>Chisq)']:
            if col in df.columns:
                pval_col = col
                break

        if pval_col is None:
            self.logger.warning("No recognized p-value column found")
            return [], pd.DataFrame()

        # Handle different logFC column names (for MAST compatibility)
        logfc_col = None
        for col in ['logfoldchanges', 'coef_logFC']:
            if col in df.columns:
                logfc_col = col
                break

        if logfc_col is None:
            # If no logFC column, just filter by p-value
            significant = df[df[pval_col] < pval_threshold].copy()
        else:
            significant = df[
                (df[pval_col] < pval_threshold) &
                (abs(df[logfc_col]) > logfc_threshold)
            ].copy()

        # Extract gene names - handle both index and 'gene'/'primerid' columns
        if 'primerid' in significant.columns:
            gene_list = significant['primerid'].dropna().tolist()
        elif 'gene' in significant.columns:
            gene_list = significant['gene'].dropna().tolist()
        else:
            gene_list = significant.index.tolist()

        self.logger.info(f"Found {len(gene_list)} significant genes")

        return gene_list, significant

    def separate_genes_by_direction(self, df: pd.DataFrame,
                                   pval_threshold: float = 0.05,
                                   direction_threshold: float = 0.25) -> Tuple[List[str], List[str]]:
        """Separate genes into upregulated and downregulated for directional enrichment"""

        # Handle different p-value column names
        pval_col = None
        for col in ['pvals_adj', 'Pr(>Chisq)']:
            if col in df.columns:
                pval_col = col
                break

        if pval_col is None:
            self.logger.warning("No p-value column found for directional analysis")
            return [], []

        # Handle different logFC column names
        logfc_col = None
        for col in ['logfoldchanges', 'coef_logFC']:
            if col in df.columns:
                logfc_col = col
                break

        if logfc_col is None:
            self.logger.warning("No logFC column found for directional analysis")
            return [], []

        # Filter significant genes
        significant = df[df[pval_col] < pval_threshold].copy()

        # Separate by direction
        upregulated = significant[significant[logfc_col] > direction_threshold].copy()
        downregulated = significant[significant[logfc_col] < -direction_threshold].copy()

        # Extract gene lists
        def extract_genes(df_subset):
            if 'primerid' in df_subset.columns:
                genes = df_subset['primerid'].dropna().tolist()
            elif 'gene' in df_subset.columns:
                genes = df_subset['gene'].dropna().tolist()
            else:
                genes = df_subset.index.tolist()
            return [str(g) for g in genes if g and str(g) != 'nan']

        up_genes = extract_genes(upregulated)
        down_genes = extract_genes(downregulated)

        self.logger.info(f"Upregulated genes (logFC > {direction_threshold}): {len(up_genes)}")
        self.logger.info(f"Downregulated genes (logFC < -{direction_threshold}): {len(down_genes)}")

        return up_genes, down_genes

    def prepare_ranked_gene_list(self, df: pd.DataFrame) -> pd.DataFrame:
        """Prepare ranked gene list for GSEA"""
        df = df.copy()
        df['ranking_metric'] = np.sign(df['logfoldchanges']) * -np.log10(df['pvals'] + 1e-300)
        ranked_genes = df[['ranking_metric']].sort_values('ranking_metric', ascending=False)
        return ranked_genes

    def run_enrichr_analysis(self, gene_list: List[str], direction: str = None) -> Dict:
        """Run Enrichr over-representation analysis

        Args:
            gene_list: List of gene symbols
            direction: Optional direction label ('upregulated', 'downregulated', or None)
        """
        if not self.available:
            return {}

        results = {}
        # Updated default to MSigDB_Hallmark_2020 for more interpretable results
        gene_sets = self.config.get('pathway_enrichment', 'gene_sets',
                                   default=['MSigDB_Hallmark_2020'])
        enrichr_cutoff = self.config.get('pathway_enrichment', 'enrichr_cutoff', default=0.05)
        min_genes = self.config.get('pathway_enrichment', 'min_genes', default=5)

        if len(gene_list) < min_genes:
            direction_label = f" {direction}" if direction else ""
            self.logger.info(f"Too few{direction_label} genes ({len(gene_list)}) for enrichment")
            return {}

        for gene_set in gene_sets:
            try:
                direction_label = f" ({direction})" if direction else ""
                self.logger.info(f"Running Enrichr for {gene_set}{direction_label}")
                enr = self.gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=gene_set,
                    organism='Human',
                    cutoff=enrichr_cutoff
                )

                if not enr.results.empty:
                    results[gene_set] = enr.results
                    self.logger.info(f"Found {len(enr.results)} enriched terms")

            except Exception as e:
                self.logger.error(f"Error in Enrichr {gene_set}: {e}")

        return results

    def run_directional_enrichment(self, deg_df: pd.DataFrame) -> Dict:
        """Run separate enrichment for upregulated and downregulated genes

        Returns:
            Dict with keys 'upregulated', 'downregulated', 'combined'
        """
        pval_threshold = self.config.get('pathway_enrichment', 'pval_threshold', default=0.05)
        direction_threshold = self.config.get('pathway_enrichment', 'direction_threshold', default=0.25)

        # Get directional gene lists
        up_genes, down_genes = self.separate_genes_by_direction(
            deg_df, pval_threshold, direction_threshold
        )

        # Get combined gene list
        all_genes, _ = self.filter_significant_genes(deg_df)

        results = {}

        # Run enrichment for each direction
        if len(up_genes) >= 5:
            self.logger.info("Running enrichment for upregulated genes")
            results['upregulated'] = self.run_enrichr_analysis(up_genes, 'upregulated')
        else:
            self.logger.info("Too few upregulated genes for enrichment")
            results['upregulated'] = {}

        if len(down_genes) >= 5:
            self.logger.info("Running enrichment for downregulated genes")
            results['downregulated'] = self.run_enrichr_analysis(down_genes, 'downregulated')
        else:
            self.logger.info("Too few downregulated genes for enrichment")
            results['downregulated'] = {}

        if len(all_genes) >= 5:
            self.logger.info("Running combined enrichment")
            results['combined'] = self.run_enrichr_analysis(all_genes)
        else:
            self.logger.info("Too few genes for combined enrichment")
            results['combined'] = {}

        return results

    def run_gsea_analysis(self, ranked_genes: pd.DataFrame) -> Dict:
        """Run GSEA analysis"""
        if not self.available:
            return {}

        results = {}
        gene_sets = self.config.get('pathway_enrichment', 'gsea_gene_sets',
                                   default=['GO_Biological_Process_2023', 'KEGG_2021_Human'])
        min_size = self.config.get('pathway_enrichment', 'gsea_min_size', default=15)
        max_size = self.config.get('pathway_enrichment', 'gsea_max_size', default=500)
        permutations = self.config.get('pathway_enrichment', 'gsea_permutations', default=1000)

        for gene_set in gene_sets:
            try:
                self.logger.info(f"Running GSEA for {gene_set}")
                gsea_res = self.gp.gsea(
                    data=ranked_genes,
                    gene_sets=gene_set,
                    processes=4,
                    min_size=min_size,
                    max_size=max_size,
                    permutation_num=permutations
                )

                if not gsea_res.res2d.empty:
                    results[gene_set] = gsea_res.res2d
                    self.logger.info(f"Found {len(gsea_res.res2d)} significant pathways")

            except Exception as e:
                self.logger.error(f"Error in GSEA {gene_set}: {e}")

        return results

    def create_enrichment_plots(self, results: Dict, output_dir: Path, prefix: str):
        """Create visualization plots for enrichment results"""
        if not self.available:
            return

        output_dir.mkdir(parents=True, exist_ok=True)
        plot_top_n = self.config.get('pathway_enrichment', 'plot_top_n', default=20)

        import matplotlib.pyplot as plt

        for gene_set, df in results.items():
            if df.empty:
                continue

            # Get top pathways
            plot_df = df.head(plot_top_n).copy()

            # Determine p-value column
            if 'Adjusted P-value' in plot_df.columns:
                plot_df['-Log10(Adj P-val)'] = -np.log10(plot_df['Adjusted P-value'] + 1e-300)
                y_col = '-Log10(Adj P-val)'
            elif 'FDR q-val' in plot_df.columns:
                plot_df['-Log10(FDR)'] = -np.log10(plot_df['FDR q-val'] + 1e-300)
                y_col = '-Log10(FDR)'
            else:
                continue

            # Create plot
            plt.figure(figsize=(12, 8))
            plot_df['Term_short'] = plot_df['Term'].apply(
                lambda x: x[:60] + '...' if len(x) > 60 else x
            )

            bars = plt.barh(range(len(plot_df)), plot_df[y_col])
            colors = plt.cm.viridis(plot_df[y_col] / plot_df[y_col].max())
            for bar, color in zip(bars, colors):
                bar.set_color(color)

            plt.yticks(range(len(plot_df)), plot_df['Term_short'])
            plt.xlabel(y_col)
            plt.title(f'Top {plot_top_n} Enriched Pathways - {gene_set}')
            plt.gca().invert_yaxis()
            plt.tight_layout()

            plot_file = output_dir / f"{prefix}_{gene_set.replace(' ', '_')}_enrichment.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()

            self.logger.info(f"Saved plot: {plot_file.name}")

    def save_results(self, results: Dict, output_dir: Path, prefix: str):
        """Save enrichment results to CSV files

        Handles both standard results and directional results
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        for key, value in results.items():
            # Check if this is a directional result (nested dict)
            if isinstance(value, dict) and key in ['upregulated', 'downregulated', 'combined']:
                direction_suffix = f"_{key}"
                for gene_set, df in value.items():
                    if df is not None and not df.empty:
                        output_file = output_dir / f"{prefix}_{gene_set.replace(' ', '_')}{direction_suffix}_enrichment.csv"
                        df.to_csv(output_file, index=False)
                        self.logger.info(f"Saved {key} results: {output_file.name}")
            # Standard result (gene_set -> df)
            elif not isinstance(value, dict) and value is not None:
                if hasattr(value, 'empty') and not value.empty:
                    output_file = output_dir / f"{prefix}_{key.replace(' ', '_')}_enrichment.csv"
                    value.to_csv(output_file, index=False)
                    self.logger.info(f"Saved results: {output_file.name}")

    def analyze_comparison(self, deg_dir: Path, output_dir: Path, comparison_name: str):
        """Analyze all DEG results for a comparison"""
        self.logger.info(f"Running pathway enrichment for {comparison_name}")

        # Find all DEG result files
        deg_files = list(deg_dir.glob("*_deg.csv"))

        if not deg_files:
            self.logger.warning(f"No DEG files found in {deg_dir}")
            return

        run_gsea = self.config.get('pathway_enrichment', 'run_gsea', default=True)
        run_directional = self.config.get('pathway_enrichment', 'run_directional', default=True)

        for deg_file in deg_files:
            cell_type = deg_file.stem.split('_')[0]
            self.logger.info(f"Processing {cell_type}")

            # Load DEG results
            deg_df = self.load_deg_results(deg_file)
            if deg_df.empty:
                continue

            enrichment_results = {}

            # Run directional enrichment if enabled
            if run_directional:
                self.logger.info(f"  Running directional enrichment for {cell_type}")
                directional_results = self.run_directional_enrichment(deg_df)
                enrichment_results.update(directional_results)
            else:
                # Standard enrichment (backward compatible)
                gene_list, significant_df = self.filter_significant_genes(deg_df)

                if len(gene_list) < self.config.get('pathway_enrichment', 'min_genes', default=5):
                    self.logger.info(f"  Skipping {cell_type}: too few genes")
                    continue

                enrichment_results = self.run_enrichr_analysis(gene_list)

            # Run GSEA if enabled
            if run_gsea and len(deg_df) > 100:
                ranked_genes = self.prepare_ranked_gene_list(deg_df)
                gsea_results = self.run_gsea_analysis(ranked_genes)
                enrichment_results.update({f"GSEA_{k}": v for k, v in gsea_results.items()})

            # Save results
            if enrichment_results:
                prefix = f"{cell_type}_pathway"
                self.save_results(enrichment_results, output_dir, prefix)
                self.create_enrichment_plots(enrichment_results, output_dir, prefix)


# ============================================================================
# VISUALIZATION GENERATOR CLASS
# ============================================================================

class VisualizationGenerator:
    """Generate visualizations for DEG and pathway results"""

    def __init__(self, config: DEGConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def create_volcano_plot(self, deg_results: pd.DataFrame, output_path: Path,
                          cell_type: str, comparison: str):
        """Create volcano plot for DEG results"""
        try:
            import matplotlib.pyplot as plt

            if deg_results.empty or 'logfoldchanges' not in deg_results.columns:
                return

            logfc_thresh = self.config.get('analysis', 'logfc_threshold', default=0.25)
            pval_thresh = self.config.get('analysis', 'pval_threshold', default=0.05)

            # Prepare data
            df = deg_results.copy()
            df['-log10(p)'] = -np.log10(df['pvals_adj'] + 1e-300)

            # Determine significance
            df['Significant'] = 'No'
            df.loc[(abs(df['logfoldchanges']) >= logfc_thresh) &
                   (df['pvals_adj'] < pval_thresh), 'Significant'] = 'Yes'

            # Create plot
            plt.figure(figsize=(10, 8))
            colors = {'Yes': 'red', 'No': 'gray'}

            for sig_status, color in colors.items():
                subset = df[df['Significant'] == sig_status]
                plt.scatter(subset['logfoldchanges'], subset['-log10(p)'],
                          c=color, label=sig_status, alpha=0.6, s=20)

            plt.axhline(-np.log10(pval_thresh), color='blue', linestyle='--', alpha=0.5)
            plt.axvline(-logfc_thresh, color='blue', linestyle='--', alpha=0.5)
            plt.axvline(logfc_thresh, color='blue', linestyle='--', alpha=0.5)

            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-Log10(Adjusted P-value)')
            plt.title(f'Volcano Plot: {cell_type}\n{comparison}')
            plt.legend()
            plt.tight_layout()

            plt.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close()

            self.logger.info(f"Saved volcano plot: {output_path.name}")

        except Exception as e:
            self.logger.error(f"Error creating volcano plot: {e}")

    def generate_summary_report(self, results_dir: Path, output_path: Path):
        """Generate HTML summary report"""
        try:
            html_content = f"""
            <!DOCTYPE html>
            <html>
            <head>
                <title>DEG Analysis Summary</title>
                <style>
                    body {{ font-family: Arial, sans-serif; margin: 40px; }}
                    h1, h2 {{ color: #333; }}
                    table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
                    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                    th {{ background-color: #f2f2f2; }}
                    .summary {{ background-color: #f9f9f9; padding: 15px; margin: 20px 0; }}
                </style>
            </head>
            <body>
                <h1>DEG and Pathway Enrichment Analysis Summary</h1>
                <div class="summary">
                    <h3>Analysis Overview</h3>
                    <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
                    <p>Results directory: {results_dir}</p>
                </div>
            </body>
            </html>
            """

            with open(output_path, 'w') as f:
                f.write(html_content)

            self.logger.info(f"Generated HTML summary: {output_path}")

        except Exception as e:
            self.logger.error(f"Error generating summary report: {e}")


# ============================================================================
# MAIN PIPELINE CLASS
# ============================================================================

class DEGPathwayPipeline:
    """Main pipeline orchestrator for DEG + Pathway enrichment analysis"""

    def __init__(self, config_path: str):
        self.config = DEGConfig(config_path)
        self.setup_logging()

        # Initialize components
        self.data_loader = DataLoader(self.config, self.logger)
        self.deg_analyzer = DEGAnalyzer(self.config, self.data_loader, self.logger)
        self.pathway_analyzer = PathwayEnrichmentAnalyzer(self.config, self.logger)
        self.viz_generator = VisualizationGenerator(self.config, self.logger)

    def setup_logging(self):
        """Setup logging configuration"""
        output_dir = Path(self.config.get('data', 'output_dir', default='./results'))
        log_dir = output_dir / 'logs'
        log_dir.mkdir(parents=True, exist_ok=True)

        log_file = log_dir / f"deg_pathway_analysis_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.log"

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Logging to {log_file}")

    def run_deg_analysis(self) -> Dict:
        """Run DEG analysis stage"""
        self.logger.info("=" * 60)
        self.logger.info("STAGE 1: DIFFERENTIAL GENE EXPRESSION ANALYSIS")
        self.logger.info("=" * 60)

        # Load and prepare data
        self.logger.info("Loading data...")
        self.data_loader.load_h5ad_data()
        self.data_loader.validate_data()
        self.data_loader.prepare_for_deg()

        # Run DEG analysis
        output_dir = Path(self.config.get('data', 'output_dir')) / '01_DEG_Analysis'
        results = self.deg_analyzer.run_all_comparisons(output_dir)

        self.logger.info(f"DEG analysis completed: {len(results)} result sets")
        return results

    def run_pathway_enrichment(self):
        """Run pathway enrichment stage"""
        if not self.config.get('pathway_enrichment', 'enabled', default=False):
            self.logger.info("Pathway enrichment disabled in config")
            return

        self.logger.info("=" * 60)
        self.logger.info("STAGE 2: PATHWAY ENRICHMENT ANALYSIS")
        self.logger.info("=" * 60)

        if not self.pathway_analyzer.available:
            self.logger.warning("gseapy not available - skipping pathway enrichment")
            return

        # Process each comparison
        deg_base_dir = Path(self.config.get('data', 'output_dir')) / '01_DEG_Analysis'
        pathway_base_dir = Path(self.config.get('data', 'output_dir')) / '02_Pathway_Enrichment'

        all_comparison_configs = self.config.get_comparison_configs()
        selected_comparisons = self.config.get('analysis', 'comparisons',
                                              default=list(all_comparison_configs.keys()))
        cell_type_level = self.config.get('analysis', 'cell_type_level', default='major')

        for comp_name in selected_comparisons:
            deg_dir = deg_base_dir / comp_name / f"{cell_type_level}_cell_types"
            output_dir = pathway_base_dir / comp_name

            if deg_dir.exists():
                self.pathway_analyzer.analyze_comparison(deg_dir, output_dir, comp_name)
            else:
                self.logger.warning(f"DEG directory not found: {deg_dir}")

    def generate_summary_report(self):
        """Generate final summary report"""
        self.logger.info("=" * 60)
        self.logger.info("GENERATING SUMMARY REPORT")
        self.logger.info("=" * 60)

        results_dir = Path(self.config.get('data', 'output_dir'))
        summary_dir = results_dir / '04_Summary_Reports'
        summary_dir.mkdir(parents=True, exist_ok=True)

        summary_file = summary_dir / 'analysis_summary.html'
        self.viz_generator.generate_summary_report(results_dir, summary_file)

    def run_complete_pipeline(self):
        """Run the complete analysis pipeline"""
        self.logger.info("STARTING DEG + PATHWAY ENRICHMENT PIPELINE")
        self.logger.info("=" * 60)

        start_time = pd.Timestamp.now()

        try:
            # Stage 1: DEG Analysis
            deg_results = self.run_deg_analysis()

            # Stage 2: Pathway Enrichment (if enabled)
            self.run_pathway_enrichment()

            # Stage 3: Summary Report
            self.generate_summary_report()

            # Final summary
            end_time = pd.Timestamp.now()
            duration = (end_time - start_time).total_seconds()

            self.logger.info("=" * 60)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY!")
            self.logger.info("=" * 60)
            self.logger.info(f"Total duration: {duration:.1f} seconds")
            self.logger.info(f"Results saved to: {self.config.get('data', 'output_dir')}")

            return {'success': True, 'results': deg_results}

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            return {'success': False, 'error': str(e)}


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

def main():
    """Main function for command-line usage"""
    parser = argparse.ArgumentParser(
        description="Unified DEG and Pathway Enrichment Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  python deg_pathway_analysis_library.py --config config.yaml

  # Run only DEG analysis
  python deg_pathway_analysis_library.py --config config.yaml --stage deg

  # Run only pathway enrichment
  python deg_pathway_analysis_library.py --config config.yaml --stage pathway
        """
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to YAML configuration file'
    )

    parser.add_argument(
        '--stage',
        type=str,
        choices=['all', 'deg', 'pathway', 'summary'],
        default='all',
        help='Analysis stage to run (default: all)'
    )

    args = parser.parse_args()

    # Initialize pipeline
    pipeline = DEGPathwayPipeline(args.config)

    # Run specified stage(s)
    if args.stage == 'all':
        results = pipeline.run_complete_pipeline()
        if results.get('success', False):
            print("\n✓ Complete pipeline finished successfully!")
            exit(0)
        else:
            print(f"\n✗ Pipeline failed: {results.get('error', 'Unknown error')}")
            exit(1)

    elif args.stage == 'deg':
        pipeline.run_deg_analysis()
        print("\n✓ DEG analysis completed!")

    elif args.stage == 'pathway':
        pipeline.run_pathway_enrichment()
        print("\n✓ Pathway enrichment completed!")

    elif args.stage == 'summary':
        pipeline.generate_summary_report()
        print("\n✓ Summary report generated!")


if __name__ == "__main__":
    main()
