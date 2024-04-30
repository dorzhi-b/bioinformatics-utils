import numpy as np
from concurrent.futures import ThreadPoolExecutor
from sklearn.tree import (DecisionTreeClassifier)


class RandomForestClassifierCustom:
    def __init__(
        self, n_estimators=10, max_depth=None, max_features=None, random_state=None
    ):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state

        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs=1):
        self.classes_ = np.unique(y)
        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            futures = []
            for i in range(self.n_estimators):
                np.random.seed(self.random_state + i)

                feat_ids = np.random.choice(range(X.shape[1]), self.max_features, replace=False)
                self.feat_ids_by_tree.append(feat_ids)

                indices = np.random.choice(range(X.shape[0]), X.shape[0], replace=True)
                X_bootstrap = X[indices]
                y_bootstrap = y[indices]

                tree = DecisionTreeClassifier(
                    max_depth=self.max_depth, random_state=self.random_state
                )
                futures.append(executor.submit(tree.fit, X_bootstrap[:, feat_ids], y_bootstrap))
                self.trees.append(tree)

            for future in futures:
                future.result()

        return self

    def predict_proba(self, X, n_jobs=1):
        probas = np.zeros((X.shape[0], len(self.classes_)))

        with ThreadPoolExecutor(max_workers=n_jobs) as executor:
            futures = []
            for i, tree in enumerate(self.trees):
                feat_ids = self.feat_ids_by_tree[i]
                futures.append(executor.submit(tree.predict_proba, X[:, feat_ids]))

            for i, future in enumerate(futures):
                probas += future.result()

        probas /= len(self.trees)

        return probas

    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs)
        predictions = np.argmax(probas, axis=1)

        return predictions
    