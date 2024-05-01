from __future__ import annotations
from typing import Optional
import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import DecisionTreeClassifier
from multiprocessing import Pool


class RandomForestClassifierCustom(BaseEstimator):
    """
    A custom implementation of Random Forest classifier.
    Parameters:
    n_estimators (int, default=10): The number of trees in the forest.
    max_depth (int or None, default=None): The maximum depth of the tree. If None, then nodes are expanded
    until all leaves are pure or until all leaves contain less than min_samples_split samples.
    max_features (int or None, default=None): The number of features to consider when looking for the best split.
    random_state (int or None, default=None): Controls the randomness of the bootstrapping and the individual trees.
    Attributes:
    classes_ (list): The classes labels.
    trees (list): The list of DecisionTreeClassifier models.
    feat_ids_by_tree (list): The list of feature indices used by each tree.
    """

    def __init__(
            self, n_estimators: int = 10, max_depth: Optional[int] = None,
            max_features: Optional[int] = None, random_state: Optional[int] = None
    ) -> None:
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def _fit_single_tree(self, args) -> tuple:
        """
        Fit a single decision tree.
        Args:
        args (tuple): A tuple containing the following elements:
            X (numpy.ndarray): The input features.
            y (numpy.ndarray): The target labels.
            max_depth (int): The maximum depth of the tree.
            max_features (int): The maximum number of features to consider when splitting.
            random_state (int): The random seed for reproducibility.
            i (int): The index of the tree.
        Returns:
        tuple: A tuple containing the trained decision tree and the indices of the selected features.
        """
        X, y, max_depth, max_features, random_state, i = args
        np.random.seed(random_state + i)
        n_samples, n_features = X.shape

        np.random.seed(self.random_state + i)

        feature_idx = np.random.choice(range(n_features),
                                       size=n_features if self.max_features is None else self.max_features,
                                       replace=False)
        bootstrap_idx = np.random.choice(range(n_samples), size=n_samples, replace=True)
        X_bootstrap = X[bootstrap_idx][:, feature_idx]
        y_bootstrap = y[bootstrap_idx]

        tree = DecisionTreeClassifier(
            max_depth=self.max_depth,
            max_features=self.max_features,
            random_state=self.random_state + i
        )
        tree.fit(X_bootstrap, y_bootstrap)
        self.trees.append(tree)

        return tree, feature_idx


    def fit(self, X: np.ndarray, y: np.ndarray, n_jobs: int = 1) -> RandomForestClassifierCustom:
        """
        Fit the Random Forest classifier to the training data.
        Parameters:
        X (array-like of shape (n_samples, n_features)): The input samples.
        y (array-like of shape (n_samples,)): The target values.
        n_jobs (int, default=1): The number of processes to use for parallel training.
        If n_jobs=-1, all available processes will be used.
        Returns:
        self (RandomForestClassifierCustom): The fitted estimator.
        """
        self.trees = []
        self.feat_ids_by_tree = []
        self.classes_ = sorted(np.unique(y))

        with Pool(n_jobs) as pool:
            results = pool.map(self._fit_single_tree,
                                   [(X, y, self.max_depth, self.max_features, self.random_state, i) for i in
                                    range(self.n_estimators)])

        self.trees, self.feat_ids_by_tree = zip(*results)

        return self

    def _predict_proba_single_tree(self, args) -> np.ndarray:
        """
        Predict class probabilities for X using a single decision tree.
        Parameters:
        args (tuple): A tuple containing the following elements:
            tree (DecisionTreeClassifier): The decision tree model.
            feature_idx (numpy.ndarray): The feature indices used by the tree.
            X (numpy.ndarray): The input samples.
        Returns:
        probas (numpy.ndarray): The class probabilities of the input samples predicted by the tree.
        """
        tree, feature_idx, X = args
        X_subset = X[:, feature_idx]
        single_tree_probas = tree.predict_proba(X_subset)
        return single_tree_probas

    def predict_proba(self, X: np.ndarray, n_jobs: int = 1) -> np.ndarray:
        """
        Predict class probabilities for X.
        Parameters:
        X (array-like of shape (n_samples, n_features)): The input samples.
        n_jobs (int, default=1): The number of processes to use for parallel prediction.
        If n_jobs=-1, all available processes will be used.
        Returns:
        probas (array-like of shape (n_samples, n_classes)): The class probabilities of the input samples.
        """
        with Pool(n_jobs) as pool:
            probas = pool.map(self._predict_proba_single_tree,
                             [(tree, feature_idx, X) for tree, feature_idx in zip(self.trees, self.feat_ids_by_tree)]
                              )

        avg_probas = np.sum(probas, axis=0) / self.n_estimators
        return avg_probas

    def predict(self, X: np.ndarray) -> np.ndarray:
        """
        Predict class labels for X.
        Parameters:
        X (array-like of shape (n_samples, n_features)): The input samples.
        Returns:
        predictions (array-like of shape (n_samples,)): The predicted class labels.
        """
        probas = self.predict_proba(X)
        predictions = np.argmax(probas, axis=1)

        return predictions
