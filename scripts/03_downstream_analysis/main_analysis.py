"""
下游分析主程序
包括数据加载、质控、预处理、聚类和注释
"""

import argparse
import os
import yaml
from pathlib import Path
from loguru import logger
import numpy as np
import pandas as pd
import scanpy as sc


def setup_logging(log_dir):
    """设置日志"""
    os.makedirs(log_dir, exist_ok=True)
    logger.add(
        os.path.join(log_dir, "downstream_analysis_{time}.log"),
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level="INFO"
    )


def load_config(config_path):
    """加载配置文件"""
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config


def load_data(data_path):
    """
    加载单细胞数据
    
    Parameters
    ----------
    data_path : str
        数据文件路径 (支持.h5ad, .h5, .csv等)
        
    Returns
    -------
    adata : AnnData
        AnnData对象
    """
    logger.info(f"加载数据: {data_path}")
    
    if data_path.endswith('.h5ad'):
        adata = sc.read_h5ad(data_path)
    elif data_path.endswith('.h5'):
        adata = sc.read_h5ad(data_path)
    elif data_path.endswith('.csv'):
        adata = sc.read_csv(data_path, first_column_names=True)
    else:
        raise ValueError(f"不支持的文件格式: {data_path}")
    
    logger.info(f"数据形状: {adata.n_obs} 细胞 × {adata.n_vars} 基因")
    
    return adata


def quality_control(adata, config):
    """
    质量控制
    
    Parameters
    ----------
    adata : AnnData
        AnnData对象
    config : dict
        配置参数
        
    Returns
    -------
    adata : AnnData
        过滤后的AnnData对象
    """
    logger.info("执行质量控制")
    
    qc_config = config['qc']
    
    # 计算QC指标
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    
    logger.info(f"原始: {adata.n_obs} 细胞")
    
    # 基于基因数过滤
    adata = adata[
        (adata.obs['n_genes_by_counts'] >= qc_config['min_genes']) &
        (adata.obs['n_genes_by_counts'] <= qc_config['max_genes'])
    ]
    logger.info(f"基因过滤后: {adata.n_obs} 细胞")
    
    # 基于线粒体基因过滤
    adata = adata[
        adata.obs['pct_counts_mt'] <= qc_config['mt_percent_threshold']
    ]
    logger.info(f"线粒体过滤后: {adata.n_obs} 细胞")
    
    return adata


def preprocessing(adata, config):
    """
    数据预处理（正规化、高变基因、降维）
    
    Parameters
    ----------
    adata : AnnData
        AnnData对象
    config : dict
        配置参数
        
    Returns
    -------
    adata : AnnData
        预处理后的AnnData对象
    """
    logger.info("执行预处理")
    
    # 正规化
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # 高变基因选择
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=config['preprocessing']['n_top_genes']
    )
    
    # PCA
    sc.tl.pca(
        adata,
        n_comps=config['dimensionality_reduction']['pca']['n_components']
    )
    
    # UMAP
    sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
    sc.tl.umap(adata)
    
    logger.info("预处理完成")
    
    return adata


def clustering(adata, config):
    """
    聚类分析
    
    Parameters
    ----------
    adata : AnnData
        AnnData对象
    config : dict
        配置参数
        
    Returns
    -------
    adata : AnnData
        包含聚类结果的AnnData对象
    """
    logger.info("执行聚类")
    
    # Leiden聚类
    sc.tl.leiden(
        adata,
        resolution=config['clustering']['resolution']
    )
    
    logger.info(f"聚类完成，得到 {adata.obs['leiden'].nunique()} 个cluster")
    
    return adata


def main():
    parser = argparse.ArgumentParser(description="单细胞下游分析")
    parser.add_argument(
        "--config",
        required=True,
        help="配置文件路径"
    )
    parser.add_argument(
        "--input",
        help="输入数据文件路径"
    )
    
    args = parser.parse_args()
    
    # 加载配置
    config = load_config(args.config)
    
    # 设置日志
    setup_logging(config['paths']['logs_dir'])
    
    logger.info("=" * 50)
    logger.info("单细胞下游分析流程启动")
    logger.info("=" * 50)
    
    # TODO: 添加完整的分析流程
    logger.info("分析完成！")


if __name__ == "__main__":
    main()
