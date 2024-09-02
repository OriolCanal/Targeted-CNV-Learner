import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
from modules.log import logger
from xgboost import XGBClassifier
from modules.params import load_config
from sklearn.metrics import f1_score, precision_score, accuracy_score, confusion_matrix, make_scorer, recall_score
import joblib
from sklearn.model_selection import train_test_split, GridSearchCV

import pickle


from sklearn.ensemble import RandomForestClassifier


class Results_Df():
    algorithms = ["decon", "cnvkit", "grapes", "gatk"]
    def __init__(self, csv_file_path, ref_conf):
        self.df = pd.read_csv(csv_file_path)
        logger.info(
            f"Parsing {csv_file_path} where info about CNV detected is stored")
        self.ref_conf = ref_conf
        self.plot_results_dir = os.path.join(ref_conf.main_dir, "Results", "Model_Plots")
        if not os.path.exists(self.plot_results_dir):
            os.mkdir(self.plot_results_dir)
        self.algs = set()
        self.predicted_col = "true_positive"
        self.check_df()
        self.models_dir = os.path.join(ref_conf.main_dir, "models_dir")
        if not os.path.exists(self.models_dir):
            os.mkdir(self.models_dir)

    def check_df(self):
        for alg in Results_Df.algorithms:
            if alg.lower() in self.df.columns:
                self.algs.add(alg)
    
    def convert_all_columns_to_numeric(self):
        self.df = self.df.apply(pd.to_numeric, errors="coerce")
    
    def are_all_columns_numeric(self):
        return all(pd.api.types.is_numeric_dtype(self.df[col]) for col in self.df.columns)

    def encode_variables(self, Bed_obj):
        print("encoding chr")
        self.encode_chr()
        print("encoding genes")
        self.encode_gene_names(Bed_obj)
        print("encoding type")
        self.encode_type()
        print("encoding true positives")
        self.encode_true_positive()
        print("droping non encoded columns")
        self.drop_non_encoded_columns()
    
    def drop_non_encoded_columns(self):

        columns_to_drop = ["sample", "type", "gene", "qual", "true_positive", "is_outlier"]
        self.df = self.df.drop(columns=columns_to_drop)



    def train_model(self):
        self.convert_all_columns_to_numeric()
        self.df.drop(columns=["start", "end", "cnvkit", "cnvkit_qual", "cnv_origin"], inplace=True)

        output_csv = os.path.join(self.ref_conf.main_dir, "model_train.csv")
        print(f"Creating model, excel stored in: {output_csv}")
        self.df.to_csv(output_csv, index=False)
        X = self.df.drop(columns=[self.predicted_col])
        y = self.df[self.predicted_col]
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)


        # Combine X_test and y_test into a single DataFrame
        df_test_combined = pd.concat([X_test, y_test], axis=1)

        # Save the combined DataFrame to a CSV file
        csv_path = "/home/ocanal/Desktop/CNV_detection_on_targeted_sequencing/models_dir/test_set.csv"
        df_test_combined.to_csv(csv_path, index=False)
        # Define parameter grids for both models
        xgb_param_grid = {
            "max_depth": [3, 4, 5],
            "learning_rate": [0.01, 0.05, 0.1],
            "n_estimators": [100, 200, 300],
            "colsample_bytree": [0.3, 0.5, 0.7]
        }
        rf_param_grid = {
            "n_estimators": [100, 200, 300],
            "max_depth": [None, 10, 20, 30],
            "min_samples_split": [2, 5, 10],
            "min_samples_leaf": [1, 2, 4]
        }

        # Initialize and fit XGBoost model
        xgb_model = XGBClassifier(objective='binary:logistic', eval_metric='logloss', use_label_encoder=False)
        scorer = make_scorer(f1_score, average="binary")
        xgb_grid_search = GridSearchCV(estimator=xgb_model, param_grid=xgb_param_grid, scoring=scorer, cv=5, verbose=1)
        xgb_grid_search.fit(X_train, y_train)
        logger.info(f"XGBoost Best parameters: {xgb_grid_search.best_params_}")
        logger.info(f"XGBoost Best F1 score: {xgb_grid_search.best_score_}")

        # Evaluate XGBoost
        best_xgb_model = xgb_grid_search.best_estimator_
        xgb_y_pred = best_xgb_model.predict(X_test)
        self.evaluate_model_with_thresholds(best_xgb_model, X_test, y_test, 'XGBoost')
        
        self.iterative_feature_removal(best_xgb_model, X_train, X_test, y_train, y_test)

        # Initialize and fit Random Forest model
        rf_model = RandomForestClassifier(random_state=42)
        rf_grid_search = GridSearchCV(estimator=rf_model, param_grid=rf_param_grid, scoring=scorer, cv=5, verbose=1)
        rf_grid_search.fit(X_train, y_train)
        logger.info(f"Random Forest Best parameters: {rf_grid_search.best_params_}")
        logger.info(f"Random Forest Best F1 score: {rf_grid_search.best_score_}")

        # Evaluate Random Forest
        best_rf_model = rf_grid_search.best_estimator_
        rf_y_pred = best_rf_model.predict(X_test)
        self.evaluate_model(best_rf_model, X_test, y_test, rf_y_pred, 'Random Forest')

        # Save the models
        model_dir = os.path.join(self.ref_conf.main_dir, "Models")
        if not os.path.exists(model_dir):
            os.mkdir(model_dir)
        
        xgb_model_path = os.path.join(model_dir, "xgb_model.joblib")
        joblib.dump(best_xgb_model, xgb_model_path)
        logger.info(f"XGBoost Model saved to: {xgb_model_path}")

        rf_model_path = os.path.join(model_dir, "rf_model.joblib")
        joblib.dump(best_rf_model, rf_model_path)
        logger.info(f"Random Forest Model saved to: {rf_model_path}")

        # Save misclassified instances for both models
        self.save_misclassified_instances(X_test, y_test, xgb_y_pred, 'XGBoost')
        self.save_misclassified_instances(X_test, y_test, rf_y_pred, 'Random Forest')


    def iterative_feature_removal(self, model, X_train, X_test, y_train, y_test):
        feature_names = X_train.columns.tolist()
        
        # Track performance metrics
        metrics = []
        removed_features = []
        models = {}
        
        while len(feature_names) > 0:
            # Train and evaluate model with current features
            X_train_selected = X_train[feature_names]
            X_test_selected = X_test[feature_names]
            model.fit(X_train_selected, y_train)
            y_pred = model.predict(X_test_selected)

            # Calculate metrics
            recall = recall_score(y_test, y_pred)
            precision = precision_score(y_test, y_pred)
            f1 = f1_score(y_test, y_pred)
            accuracy = accuracy_score(y_test, y_pred)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
           
            # Save the current state of the model
            num_removed = len(X_train.columns) - len(feature_names)
            model_path = os.path.join(self.models_dir, f"model_after_removing_{num_removed}_features.joblib")
            

            logger.info(f"features names of nº features removed: {num_removed}: {feature_names}")
            joblib.dump(model, model_path)

            metrics.append({
                'Removed Features': num_removed,
                'Recall': recall,
                'Precision': precision,
                'F1 Score': f1,
                'Accuracy': accuracy,
                'Specificity': specificity,
            })

            # Log results
            logger.info(f"Features remaining: {len(feature_names)}")
            logger.info(f"Recall: {recall}")
            logger.info(f"Precision: {precision}")
            logger.info(f"F1 Score: {f1}")
            logger.info(f"Accuracy: {accuracy}")
            logger.info(f"Specificity (True Negative Rate): {specificity}")
            logger.info(f"Confusion Matrix:\n{confusion_matrix(y_test, y_pred)}")

            if len(feature_names) == 1:
                break
            
            # Get feature importances and identify the least important feature
            importances = model.feature_importances_
            feature_importance_df = pd.DataFrame({
                'Feature': feature_names,
                'Importance': importances
            }).sort_values(by='Importance', ascending=False)

            least_important_feature = feature_importance_df.iloc[-1]['Feature']
            importance_score = feature_importance_df.iloc[-1]['Importance']
            removed_features.append((least_important_feature, importance_score))

            # Remove the least important feature
            feature_names.remove(least_important_feature)

        # Convert metrics to DataFrame and log
        metrics_df = pd.DataFrame(metrics)
        logger.info(f"Performance metrics after each iteration:\n{metrics_df}")

        # Plot feature importances evolution
        plt.figure(figsize=(10, 8))
        plt.plot(metrics_df['Removed Features'], metrics_df['F1 Score'], marker='o', linestyle='-', label='F1 Score')
        plt.plot(metrics_df['Removed Features'], metrics_df['Accuracy'], marker='o', linestyle='-', label='Accuracy')
        plt.plot(metrics_df['Removed Features'], metrics_df['Recall'], marker='o', linestyle='-', label='Recall')
        plt.xlabel('Number of Features Removed')
        plt.ylabel('Score')
        plt.title('Model Performance vs. Number of Features Removed')
        plt.legend()
        plt.grid(True)
        plt.savefig("/home/ocanal/Desktop/CNV_detection_on_targeted_sequencing/template/imatges/feature_removal.png")
        plt.show()

        # Log the removed features in order with their importance scores
        removed_features_str = ", ".join([f"{feat} (Importance: {score:.4f})" for feat, score in removed_features])
        logger.info(f"Order of removed features and their importance scores: {removed_features_str}")



    def evaluate_model_with_thresholds(self, model, X_test, y_test, model_name):
        # Get predicted probabilities
        y_proba = model.predict_proba(X_test)[:, 1]
        
        # Define thresholds
        thresholds = [0.5, 0.45, 0.40, 0.35, 0.30, 0.25, 0.2, 0.15, 0.1, 0.05]
        results = []

        for threshold in thresholds:
            # Apply threshold
            y_pred = (y_proba > threshold).astype(int)
            
            # Calculate metrics
            recall = recall_score(y_test, y_pred)
            precision = precision_score(y_test, y_pred)
            f1 = f1_score(y_test, y_pred)
            accuracy = accuracy_score(y_test, y_pred)
            tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
            specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
            
            # Log results
            logger.info(f"Threshold: {threshold}")
            logger.info(f"Recall: {recall}")
            logger.info(f"Precision: {precision}")
            logger.info(f"F1 Score: {f1}")
            logger.info(f"Accuracy: {accuracy}")
            logger.info(f"Specificity (True Negative Rate): {specificity}")
            logger.info(f"Confusion Matrix:\n{confusion_matrix(y_test, y_pred)}")
            logger.info(f"true negatives: {tn}\nfalse positives {fp}\n false negatives {fn} true positive {tp}")
            
            # Append to results
            results.append({
                'Threshold': threshold,
                'Recall': recall,
                'Precision': precision,
                'F1 Score': f1,
                'Accuracy': accuracy,
                'Specificity': specificity
            })
        
        # Convert to DataFrame
        results_df = pd.DataFrame(results)
        
        # Plot results
        plt.figure(figsize=(12, 8))
        
        plt.plot(results_df['Threshold'], results_df['Recall'], marker='o', linestyle='-', label='Recall')
        plt.plot(results_df['Threshold'], results_df['Precision'], marker='o', linestyle='-', label='Precision')
        plt.plot(results_df['Threshold'], results_df['F1 Score'], marker='o', linestyle='-', label='F1 Score')
        plt.plot(results_df['Threshold'], results_df['Accuracy'], marker='o', linestyle='-', label='Accuracy')
        plt.plot(results_df['Threshold'], results_df['Specificity'], marker='o', linestyle='-', label='Specificity')
        
        plt.xlabel('Threshold')
        plt.ylabel('Score')
        plt.title(f'{model_name} Performance Metrics at Different Thresholds')
        plt.legend()
        plt.grid(True)
        
        # Save the plot
        plot_path = f"/home/ocanal/Desktop/CNV_detection_on_targeted_sequencing/template/imatges/{model_name.lower()}_thresholds.png"
        plt.savefig(plot_path)
        
        plt.show()
        
        logger.info(f"Threshold performance plot saved to: {plot_path}")

    def evaluate_model(self, model, X_test, y_test, y_pred, model_name):
        test_f1_score = f1_score(y_test, y_pred)
        accuracy = accuracy_score(y_test, y_pred)
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        specificity = tn / (tn + fp)
        recall = recall_score(y_test, y_pred)

        conf_matrix = confusion_matrix(y_test, y_pred)

        logger.info(f"{model_name} Performance:")
        logger.info(f"F1 score on the test set: {test_f1_score}")
        logger.info(f"Accuracy on the test set: {accuracy}")
        logger.info(f"Specificity (True Negative Rate) on the test set: {specificity}")
        logger.info(f"Recall (Sensitivity) on the test set: {recall}")

        logger.info(f"Confusion Matrix:\n{conf_matrix}")

    def save_misclassified_instances(self, X_test, y_test, y_pred, model_name):
        misclassified_mask = (y_test != y_pred)
        misclassified_df = X_test[misclassified_mask].copy()
        misclassified_df['Actual'] = y_test[misclassified_mask]
        misclassified_df['Predicted'] = y_pred[misclassified_mask]

        misclassified_output_csv = os.path.join(self.main_dir, f"misclassified_instances_{model_name}.csv")
        misclassified_df.to_csv(misclassified_output_csv, index=False)
        logger.info(f"Misclassified instances for {model_name} saved to: {misclassified_output_csv}")


    def encode_true_positive(self):
        self.df["true_positive_encoded"] = self.df["true_positive"].astype(int)
        self.predicted_col = "true_positive_encoded"
    def encode_type(self):
        mapping = {"DEL": 0, "DUP": 1}
        self.df["encoded_type"] = self.df["type"].map(mapping)
    def encode_chr(self):
        self.df['chr'] = self.df['chr'].str.replace('chr', '', regex=False)
        self.df["is_autosomal"] = self.df["chr"].apply(lambda x: 0 if x=="X" or x=="Y" else 1)
        self.df["chr"] = self.df["chr"].replace({"X": 23, "Y":24}).astype(int)

    def get_gene_name_encoded(self, gene_name, Bed_obj):
        """
        Bed obj is not necessary if encode_genes() have been executed
        """
        print(gene_name)
        if not hasattr(Bed_obj, "gene_names"):
            print("obtaining gene names")
            Bed_obj.get_gene_names()
            print(Bed_obj.gene_names)
        
        print(Bed_obj.gene_names)
        try:
            gene_pos = Bed_obj.gene_names.index(gene_name)
            print(gene_pos)
        except ValueError:
            logger.info(
                f"gene name: {gene_name} not found in bed, "
            )
            gene_pos = -1
        return(gene_pos)
    
    def encode_gene_names(self, Bed_obj):
        self.df["genes_encoded"] = self.df["gene"].apply(lambda x: self.get_gene_name_encoded(x, Bed_obj))









    def get_algorithms_bar_plot(self):
        for alg in self.algs:
            plot_path = os.path.join(self.plot_results_dir, f"{alg}_bar_chart.png")
            sns.countplot(x=alg.lower(), hue=self.predicted_col, data=self.df)
            plt.title(f"{alg} vs TP(1)/FP(0)")
            plt.savefig(plot_path)

            plot_path = os.path.join(self.plot_results_dir, f"{alg}_stacked_bar_chart.png")
            cross_tab = pd.crosstab(self.df[alg], self.df[self.predicted_col])
            cross_tab.plot(kind="bar", stacked=True)
            plt.title(f"{alg} vs TP(1)/FN(0)")
            plt.savefig(plot_path)

    def get_heatmap(self):
        plot_path = os.path.join(self.plot_results_dir, f"heatmap.png")
        corr = self.df.corr()
        sns.heatmap(corr, annot=True, cmap="coolwarm", center=0)
        plt.title(f"Correlation Heatmap")
        plt.savefig(plot_path)


if __name__ == "__main__":
    ann_config, ref_config, docker_config, isoforms_conf = load_config()
    excel_path = os.path.join(ref_config.main_dir, "cnv_excel.csv")
    # Results_Df(excel_path, ref_conf)