import os
import subprocess
import gzip
from modules.log import logger
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import plotly.express as px
import umap
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.mixture import BayesianGaussianMixture
from scipy.cluster import hierarchy
from sklearn.metrics.pairwise import euclidean_distances
import numpy as np
from scipy import linalg
from matplotlib.patches import Ellipse

import sys
import hdbscan

from scipy.cluster.hierarchy import fcluster


class Mosdepth():
    def __init__(self, docker_conf):
        self.docker_path = "/usr/bin/docker"
        self.bed = docker_conf.bed
        self.bed_dir = os.path.dirname(self.bed)
        self.bed_filename = os.path.basename(self.bed)
        self.mosdepth_image = docker_conf.mosdepth["image"]
        self.mosdepth_version = docker_conf.mosdepth["version"]
        self.thresholds = [1, 10, 20, 30, 50, 70, 100, 150, 200, 300, 500, 700, 1000, 2000, 5000]


    def run_mosdepth(self, bam_file, force=False):

        bam_filename = os.path.basename(bam_file)
        sample_name = bam_filename.replace(".bam", "")
        bam_dir = os.path.dirname(bam_file)

        mosdepth_dirname = "Mosdepth"
        mosdepth_dir =  os.path.join(os.path.dirname(bam_file), mosdepth_dirname)
        
        self.bed_regions_filename = f"{sample_name}.regions.bed.gz"
        self.bed_regions_path = os.path.join(mosdepth_dir, sample_name, self.bed_regions_filename)
        self.bed_with_thresholds_filename = f"{sample_name}.thresholds.bed.gz"
        self.bed_with_thresholds_path = os.path.join(mosdepth_dir, sample_name, self.bed_with_thresholds_filename)
        if force == False:
            # check if output already available
            if os.path.exists(self.bed_regions_path):
                logger.info(f"Mosdepth output already available: {self.bed_regions_path}")
                return (self.bed_regions_path, sample_name)

        if not os.path.isdir(mosdepth_dir):
            os.mkdir(mosdepth_dir)
        
        # creating sample directory
        sample_dir = os.path.join(mosdepth_dir, sample_name)
        if not os.path.isdir(sample_dir):
            os.mkdir(sample_dir)

        thresholds_str = ""
        for threshold in self.thresholds:
            if thresholds_str:
                thresholds_str += f",{threshold}"
            else:
                thresholds_str = str(threshold)
        mosdepth_cmd = [
            self.docker_path, "run",
            "-v", f"{self.bed_dir}:/bed_dir",
            "-v", f"{bam_dir}:/bam_dir",
            f"{self.mosdepth_image}:{self.mosdepth_version}",
            "mosdepth", 
            # "--thresholds", thresholds_str,
            "--fast-mode",
            "-n", # don't do per base depth (improve speed)
            "--by", f"/bed_dir/{self.bed_filename}",
            f"/bam_dir/{mosdepth_dirname}/{sample_name}/{sample_name}",
            f"/bam_dir/{bam_filename}"
        ]


        command_str = " ".join(mosdepth_cmd)
        logger.info(
            f"Running Mosdepth on sample {sample_name}:\n{command_str}"
        )
        subprocess.run(mosdepth_cmd, encoding="utf-8", capture_output=True)
        
        
        return(self.bed_regions_path, sample_name)
    
    def parse_bed_with_thresholds(self):
        """
        Parsing file: RB36565.thresholds.bed.gz
        This function is not used as we will be using file RB36565.regions.bed.gz"""
        if not hasattr(self, "bed_with_thresholds_path"):
            raise ValueError(
                f"Run mosdepth before parsing its output file"
            )
        
        bed_filename = os.path.basename(self.bed_with_thresholds_path)
        sample_name = bed_filename.split(".")[0]
        exon_coverage = dict()

        with gzip.open(self.bed_with_thresholds_path, "rb") as f:
            for line in f:
                line = line.strip().decode()
                if line.startswith("#"):
                    # removing the #
                    line = line[1:]
                    header = line.split("\t")
                values = line.split("\t")
                chr = values[0]
                start = values[1]
                end = values[2]
                info = values[3]
                info = dict(item.split(": ") for item in info.str.split(", "))
                bed_thresholds = values[4:]
                exon = info["exon_id"]
                if len(self.thresholds) != len(values[4:]):
                    raise ValueError(
                        f"Number of thresholds defined {self.thresholds} don't match the number of thresholds in file {self.bed_with_thresholds_path}"
                    )
                exon_length = end - start
                max_threshold = 0
                for threshold, value in zip(self.thresholds, bed_thresholds):
                    if value == exon_length:
                        if threshold > max_threshold:
                            max_threshold = threshold
                
                if exon in exon_coverage:
                    raise ValueError(
                        f"Exon {exon} is repeated in {self.bed_with_thresholds_path}"
                    )
                exon_coverage[exon] = max_threshold
        exon_coverage["sample"] = sample_name
        
        return(exon_coverage)

    def parse_mosdepth_regions_bed(self, run_id):
        self.bed_regions_path

        if not hasattr(self, "bed_regions_path"):
            raise ValueError(
                f"Run mosdepth before parsing its output file"
            )
        
        bed_filename = os.path.basename(self.bed_regions_path)
        sample_name = bed_filename.split(".")[0]
        self.mean_exon_coverage = dict()
        gene_coverage = dict()
        self.mean_gene_coverage = dict()
        sample_length = 0
        sample_coverage = 0
        with gzip.open(self.bed_regions_path, "rb") as f:
            for line in f:
                line = line.strip().decode()
                if line.startswith("#"):
                    # removing the #
                    line = line[1:]
                    header = line.split("\t")
                values = line.split("\t")
                chr = values[0]
                start = values[1]
                end = values[2]
                exon_length = int(end) -1 - int(start)
                info = values[3]
                # Remove outer curly braces
                info = info.strip("{}")

                # Split by comma and space to separate key-value pairs
                pairs = info.split(", ")

                # Split each pair by colon to separate key and value, and strip spaces
                info_dict = {pair.split(": ")[0].strip("'"): pair.split(": ")[1].strip("'") for pair in pairs}


                exon_mean_coverage = values[4]
                exon = info_dict["exon_id"]
                gene_id = info_dict["gene_id"]
                
                if exon in self.mean_exon_coverage:
                    raise ValueError(
                        f"Exon {exon} is repeated in {self.bed_with_thresholds_path}"
                    )
                self.mean_exon_coverage[exon] = float(exon_mean_coverage)
                sample_length += exon_length
                sample_coverage += float(exon_length) * float(exon_mean_coverage)
        self.sample_mean_coverage = sample_coverage / sample_length
        self.mean_exon_coverage["mean_coverage"] = self.sample_mean_coverage
        self.mean_exon_coverage["sample"] = sample_name
        self.mean_exon_coverage["run_id"] = run_id
        # print(self.mean_exon_coverage)
        return(self.mean_exon_coverage, self.sample_mean_coverage)


class Mosdepth_df():
    def __init__(self, ref_conf):
        self.exon_coverage_df = pd.DataFrame()
        self.exons_mean_coverage_dicts = list()
        self.normalized_exons_mean_coverage_dicts = list()
        self.main_dir = ref_conf.main_dir
        self.plot_dir = os.path.join(self.main_dir, "Plots")
        if not os.path.exists(self.plot_dir):
            os.mkdir(self.plot_dir)

    def add_mean_coverage_dict(self, mean_coverage_dict):
        self.exons_mean_coverage_dicts.append(mean_coverage_dict)
    
    def add_normalized_mean_coverage_dict(self, mean_coverage_dict, sample_mean_coverage):
        
        not_normalized_cols = ["sample", "run_id"]
        for key in mean_coverage_dict:
            if key in not_normalized_cols:
                continue
            mean_coverage_dict[key] /= float(sample_mean_coverage)

        
        self.normalized_exons_mean_coverage_dicts.append(mean_coverage_dict)

    def get_df_from_exons_coverage(self):
        self.exon_coverage_df = pd.DataFrame(self.exons_mean_coverage_dicts)

        return(self.exon_coverage_df)
    
    def get_df_from_normalized_exons_coverage(self):
        self.normalized_exon_coverage_df = pd.DataFrame(self.normalized_exons_mean_coverage_dicts)

        return(self.normalized_exon_coverage_df)
    
    def apply_pca(self, normalized=True):
        self.pca_plot_dir = os.path.join(self.plot_dir, "PCA")
        if not os.path.exists(self.pca_plot_dir):
            os.mkdir(self.pca_plot_dir)
        if normalized:
            df_copy = self.normalized_exon_coverage_df.copy()
            pca_plot_path = os.path.join(self.pca_plot_dir, "normalized_pca.html")
        else:
            pca_plot_path = os.path.join(self.pca_plot_dir, "non_normalized_PCA.html")
            df_copy = self.exon_coverage_df.copy()

        run_id = df_copy.pop("run_id").tolist()
        sample_id = df_copy.pop("sample").tolist()
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_copy)

        self.pca = PCA(n_components=2)
        pca_result = self.pca.fit_transform(scaled_data)
        self.pca_df = pd.DataFrame(data=pca_result, columns=["PC1", "PC2"])
        
        self.pca_df["run_id"] = run_id
        self.pca_df["sample_id"] = sample_id
        logger.info(f"Creating PCA plot stored in: {pca_plot_path}")
        # Visualize the PCA result
        
        custom_palette = px.colors.qualitative.Vivid

        # Create an interactive scatter plot with color based on run_id and hover information showing sample_id
        fig = px.scatter(self.pca_df, x='PC1', y='PC2', color='run_id', hover_data=['sample_id'], title='PCA Result', color_discrete_sequence=custom_palette)

        # Update layout to show grid
        fig.update_layout(
            xaxis=dict(title='Principal Component 1'),
            yaxis=dict(title='Principal Component 2'),
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="Arial"),
            hovermode='closest',
        )
        fig.write_html(pca_plot_path)
        # Show the plot
        # fig.show()
        # Print explained variance ratio
        logger.info(f"Explained Variance Ratio: {self.pca.explained_variance_ratio_}")


    def get_n_pca_closest_samples(self, sample_id, n=10)-> pd.DataFrame:
        """
        Given a sample_id, it obtain the n points with minimum distance to sample id point
        represented in the PCA dimensional space (extracts similar samples given an id). It
        also returns the same point as the similar points obtained.
        """
        distances = euclidean_distances(self.pca_df.loc[self.pca_df['sample_id'] == sample_id, ["PC1", "PC2"]], self.pca_df[["PC1", "PC2"]])[0]

        # Get the indices of the closest points
        closest_indices = distances.argsort()[:n]

        # Filter the dataframe to get the closest points
        closest_points = self.pca_df.iloc[closest_indices]

        logger.info(f"Closest {n} points to sample {sample_id}:")
        logger.info(closest_points)
        return(closest_points)

    def apply_umap(self, normalized=True):
        self.umap_plot_dir = os.path.join(self.plot_dir, "UMAP")
        if not os.path.exists(self.umap_plot_dir):
            os.mkdir(self.umap_plot_dir)
        if normalized:
            df_copy = self.normalized_exon_coverage_df.copy()
            umap_plot_path = os.path.join(self.umap_plot_dir, "normalized_umap.html")

        else:
            df_copy = self.exon_coverage_df.copy()
            umap_plot_path = os.path.join(self.umap_plot_dir, "non_normalized_UMAP.html")

        run_id = df_copy.pop("run_id").tolist()
        sample_id = df_copy.pop("sample").tolist()
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(df_copy)

        umap_model = umap.UMAP(n_components=2, n_neighbors=15, min_dist=0.1)
        umap_result = umap_model.fit_transform(scaled_data)
        self.umap_df = pd.DataFrame(data=umap_result, columns=["PC1", "PC2"])
        self.umap_df["run_id"] = run_id
        self.umap_df["sample_id"] = sample_id
        logger.info(f"Creating UMAP plot stored in: {umap_plot_path}")
        # Visualize the PCA result
        
        custom_palette = px.colors.qualitative.Dark24

        # Create an interactive scatter plot with color based on run_id and hover information showing sample_id
        fig = px.scatter(self.umap_df, x='PC1', y='PC2', color='run_id', hover_data=['sample_id'], title='UMAP Result', color_discrete_sequence=custom_palette)

        # Update layout to show grid
        fig.update_layout(
            xaxis=dict(title='Principal Component 1'),
            yaxis=dict(title='Principal Component 2'),
            hoverlabel=dict(bgcolor="white", font_size=12, font_family="Arial"),
            hovermode='closest',
        )
        fig.write_html(umap_plot_path)

        # Show the plot
        # fig.show()


    def get_heatmap(self, heatmap_path):
        df_copy = self.exon_coverage_df.copy()
        # Convert all columns to numeric (excluding sample_id and run_id)
        numeric_columns = df_copy.columns.difference(['sample', 'run_id'])
        df_copy[numeric_columns] = df_copy[numeric_columns].apply(pd.to_numeric)
        run_id = df_copy.pop("run_id").tolist()
        # Set sample_id as index
        df_copy.set_index('sample', inplace=True)
        logger.info(f"Creating heatmap of exons coverage depth")
        # Compute hierarchical clustering for both rows (samples) and columns (exons)
        row_cluster = hierarchy.linkage(df_copy, method='average', metric='euclidean')
        col_cluster = hierarchy.linkage(df_copy.transpose(), method='average', metric='euclidean')

        # Create a clustered heatmap
        plt.figure(figsize=(24, 18))  # Adjust the figure size as needed
        sns.heatmap(df_copy.iloc[hierarchy.dendrogram(row_cluster, no_plot=True)['leaves'], 
                        hierarchy.dendrogram(col_cluster, no_plot=True)['leaves']], 
                    cmap='viridis', annot=True, fmt=".2f")

        plt.title('Clustered Heatmap of Exon Coverage Across Samples')
        plt.xlabel('Exon IDs')
        plt.ylabel('Sample IDs')
        plt.savefig(heatmap_path)
