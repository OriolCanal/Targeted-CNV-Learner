import numpy as np
import pandas as pd
import os
from sklearn.model_selection import GridSearchCV
from sklearn.mixture import BayesianGaussianMixture
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import linalg
from matplotlib.patches import Ellipse
from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
from sklearn.metrics import silhouette_score
import hdbscan
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors




class Clustering_Mosdepth():
    reduced_dimensional_method_options = ["UMAP", "PCA", "NORMALIZED", "RAW"]
    def __init__(self, Mosdepth_df, reduced_dimensional_method, ref_conf):
        self.Mosdepth_df = Mosdepth_df
        self.reduced_dimensional_method = reduced_dimensional_method
        self.check_reduced_dimensional_method()
        self.get_data()
        self.plot_dir = os.path.join(ref_conf.main_dir, "Plots", self.reduced_dimensional_method.upper())
        if not os.path.exists(self.plot_dir):
            print("Creating plot dir")
            os.mkdir(self.plot_dir)

    def check_reduced_dimensional_method(self):
        if not self.reduced_dimensional_method.upper() in self.reduced_dimensional_method_options:
            raise ValueError(
                f"Reduced_dimensional method should be in: {self.reduced_dimensional_method_options}"
            )


    def get_data(self):
        if self.reduced_dimensional_method == "PCA":
            PC1_values = self.Mosdepth_df.pca_df['PC1'].values.reshape(-1, 1)
            PC2_values = self.Mosdepth_df.pca_df['PC2'].values.reshape(-1, 1)

            # Concatenate PC1 and PC2 arrays to create a similar structure as X in the example
            self.data = np.concatenate([PC1_values, PC2_values], axis=1)
            self.df = self.Mosdepth_df.pca_df
            return (self.data)
        
        if self.reduced_dimensional_method == "UMAP":
            PC1_values = self.Mosdepth_df.umap_df["PC1"].values.reshape(-1, 1)
            PC2_values = self.Mosdepth_df.umap_df["PC2"].values.reshape(-1, 1)
            # Concatenate PC1 and PC2 arrays to create a similar structure as X in the example
            self.data = np.concatenate([PC1_values, PC2_values], axis=1)
            self.df = self.Mosdepth_df.umap_df

            return (self.data)
    
    def gmm_bic_score(self, estimator, X):
        """Callable to pass to GridSearchCV that will use the BIC score."""
        # Calculate the log-likelihood
        log_likelihood = estimator.score(X)
        # Get the number of samples and number of parameters
        n_samples, n_features = X.shape
        num_parameters = estimator.n_components * (n_features + 1) + estimator.n_components * (estimator.n_components - 1) // 2
        # Compute BIC
        bic_score = np.log(n_samples) * num_parameters - 2 * log_likelihood
        return -bic_score
    
    def apply_bayesian_gaussian_mixture(self):
        gmm_dir = os.path.join(self.plot_dir, "gmm")
        if not os.path.exists(gmm_dir):
            os.mkdir(gmm_dir)

        min_components = 2
        max_components = 10
        n_components_range = range(min_components, max_components+1)
        covariance_types = ["spherical", "tied", "full"]

        param_grid = {
            "n_components": n_components_range,
            "covariance_type": covariance_types
        }

        grid_search = GridSearchCV(BayesianGaussianMixture(), param_grid=param_grid, scoring=self.gmm_bic_score)

        # Use the data attribute directly
        X = self.data
        grid_search.fit(X)

        df = pd.DataFrame(grid_search.cv_results_)[
            ["param_n_components", "param_covariance_type", "mean_test_score"]
        ]
        df["mean_test_score"] = -df["mean_test_score"]
        df = df.rename(
            columns={
                "param_n_components": "Number of components",
                "param_covariance_type": "Type of covariance",
                "mean_test_score": "BIC score",
            }
        )

        bic_scores_path = os.path.join(gmm_dir, "bic_pca_scores.png")

        sns.catplot(
            data=df,
            kind="bar",
            x="Number of components",
            y="BIC score",
            hue="Type of covariance",
        )
        plt.savefig(bic_scores_path)

        color_iter = sns.color_palette("tab10", 2)[::-1]
        Y_ = grid_search.predict(X)
        self.hierarchical_index_path = os.path.join(gmm_dir, "cluster_indices.npy")
        np.save(self.hierarchical_index_path, Y_)
        fig, ax = plt.subplots()

        for i, (mean, cov, color) in enumerate(
            zip(
                grid_search.best_estimator_.means_,
                grid_search.best_estimator_.covariances_,
                color_iter,
            )
        ):
            if not cov.shape:
                continue
            print(cov.shape)
            num_rows, num_cols = cov.shape
            # checking that the matrix is sqaure
            if num_rows != num_cols:
                continue
            v, w = linalg.eigh(cov)
            if not np.any(Y_ == i):
                continue
            plt.scatter(X[Y_ == i, 0], X[Y_ == i, 1], 0.8, color=color)

            angle = np.arctan2(w[0][1], w[0][0])
            angle = 180.0 * angle / np.pi  # convert to degrees
            v = 2.0 * np.sqrt(2.0) * np.sqrt(v)
            ellipse = Ellipse(mean, v[0], v[1], angle=180.0 + angle, color=color)
            ellipse.set_clip_box(fig.bbox)
            ellipse.set_alpha(0.5)
            ax.add_artist(ellipse)

        plt.title(
            f"Selected GMM: {grid_search.best_params_['covariance_type']} model, "
            f"{grid_search.best_params_['n_components']} components"
        )
        plt.axis("equal")
        cluster_bic_path = os.path.join(self.plot_dir, "bic_cluster.png")

        plt.savefig(cluster_bic_path)

    def do_hierarchical_clustering(self):

        hierarchical_clustering_dir = os.path.join(self.plot_dir, "hierarchical_clustering")   
        if not os.path.exists(hierarchical_clustering_dir):
            os.mkdir(hierarchical_clustering_dir)

        # Perform hierarchical clustering
        Z = linkage(self.data, method='ward')  # You can use different linkage methods

        dendogram_path = os.path.join(hierarchical_clustering_dir, "dendogram.png")
        # Visualize the dendrogram
        plt.figure(figsize=(12, 8))
        dendrogram(Z)
        plt.xlabel('Sample Index')
        plt.ylabel('Distance')
        plt.title('Hierarchical Clustering Dendrogram')
        plt.show()
        plt.savefig(dendogram_path)

        # Optionally, you can also cut the dendrogram to obtain clusters
        # For example, you can use fcluster to cut the dendrogram at a specific height or number of clusters:
        max_d = 350  # Max distance threshold to cut the dendrogram
        clusters = fcluster(Z, max_d, criterion='distance')
        self.hierarchical_clustering_tag = f"{self.reduced_dimensional_method}_hierarchical_cluster_label"
        # Assign cluster labels to the samples
        self.df[self.hierarchical_clustering_tag] = clusters

        
        hierarchical_clustering_path = os.path.join(hierarchical_clustering_dir, "hierarchical_cluster.png")
        # Visualize the clusters on the PCA plot
        plt.figure(figsize=(10, 6))
        for cluster_label in set(clusters):
            cluster_data = self.df[self.df[self.hierarchical_clustering_tag] == cluster_label]
            plt.scatter(cluster_data['PC1'], cluster_data['PC2'], label=f'Cluster {cluster_label}')

        plt.xlabel('Principal Component 1')
        plt.ylabel('Principal Component 2')
        plt.title('PCA with Cluster Labels')
        plt.legend()
        plt.grid(True)
        plt.show()
        plt.savefig(hierarchical_clustering_path)

    
    def apply_hdbscan(self, min_cluster_size, min_samples):


        hdbscan_dir = os.path.join(self.plot_dir, "HDBSCAN")
        if not os.path.exists(hdbscan_dir):
            os.mkdir(hdbscan_dir)

        hdbscan_path = os.path.join(hdbscan_dir, f"HDBSCAN_{min_cluster_size}.clustersize_{min_samples}.minsamples.png")

        # Plot clusters using best parameters
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
        cluster_labels = clusterer.fit_predict(self.data)
        print(cluster_labels)
        print(type(cluster_labels))
        self.hdbscan_tag = f"{self.reduced_dimensional_method}_HDBSCAN_labels"
        self.df[self.hdbscan_tag] = cluster_labels
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=self.data[:, 0], y=self.data[:, 1], hue=cluster_labels, palette='viridis')
        plt.title('HDBSCAN Clustering with Best Parameters')
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.legend(title='Cluster')
        plt.savefig(hdbscan_path)

    def apply_dbscan(self, epsilon, min_pts, k_graph=False):


        dbscan_dir = os.path.join(self.plot_dir, "DBSCAN")
        if not os.path.exists(dbscan_dir):
            os.mkdir(dbscan_dir)

        dbscan_path = os.path.join(dbscan_dir, f"HDBSCAN_{min_pts}.minsamples_{epsilon}.epsilon.png")

        # Plot clusters using best parameters
        clusterer = DBSCAN(eps=epsilon, min_samples=min_pts)
        cluster_labels = clusterer.fit_predict(self.data)
        print(f"dbscan labels: {cluster_labels}")
        print(type(cluster_labels))
        self.dbscan_tag = f"{self.reduced_dimensional_method}_DBSCAN_labels"
        self.df[self.dbscan_tag] = cluster_labels
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=self.data[:, 0], y=self.data[:, 1], hue=cluster_labels, palette='viridis')
        plt.title('DBSCAN Clustering')
        plt.xlabel('PC1')
        plt.ylabel('PC2')
        plt.legend(title='Cluster')
        plt.savefig(dbscan_path)


        if not k_graph:
            return(cluster_labels)        
        
        # Step to create k-distance graph
        k_distance_path = os.path.join(dbscan_dir, f"k_distance_{min_pts}.png")
        k = min_pts
        neighbors = NearestNeighbors(n_neighbors=k)
        neighbors_fit = neighbors.fit(self.data)
        distances, indices = neighbors_fit.kneighbors(self.data)
        
        # Sort the distances to the k-th nearest neighbor
        k_distances = np.sort(distances[:, k-1], axis=0)
        
        # Plot the k-distance graph
        plt.figure(figsize=(8, 6))
        plt.plot(k_distances)
        plt.title('k-Distance Graph')
        plt.xlabel('Points sorted by distance')
        plt.ylabel(f'{k}-th Nearest Neighbor Distance')
        plt.savefig(k_distance_path)

        return(cluster_labels)

    def set_hdbscan_cluster_samples(self, cohort_samples):
        cluster_samples = dict()
        for index, row in self.df.iterrows():
            cluster_id = row[self.hdbscan_tag]
            id = row["sample_id"]

            sample = cohort_samples[id]

            sample.cluster = cluster_id

            if cluster_id not in cluster_samples:
                cluster_samples[cluster_id] = [sample]
            else:
                cluster_samples[cluster_id].append(sample)
        
        return cluster_samples
    
    def set_dbscan_cluster_samples(self, cohort_samples):

        cluster_samples = dict()
        for index, row in self.df.iterrows():
            cluster_id = row[self.dbscan_tag]
            id = row["sample_id"]

            sample = cohort_samples[id]

            sample.cluster = cluster_id

            if cluster_id not in cluster_samples:
                cluster_samples[cluster_id] = [sample]
            else:
                cluster_samples[cluster_id].append(sample)
        
        return cluster_samples
    



