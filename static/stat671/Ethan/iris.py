'''Iris
@author: Ethan Lew
'''

from sklearn.datasets import load_iris
import pandas as pd
import numpy as np

from classifier import PartitionData,  SimpleClassifier, risk

def load_iris_data(s0, s1, ratio):
    # Get Iris Data
    iris_sk = load_iris()
    species_lp = {'I. setosa': 0, 'I. versicolor': 1, 'I. virginica': 2}
    iris_df = {'sepal_length': iris_sk['data'][:, 0], 'sepal_width': iris_sk['data'][:, 1],
               'petal_length': iris_sk['data'][:, 2],
               'petal_width': iris_sk['data'][:, 3], 'species': iris_sk['target']}
    iris_df = pd.DataFrame(data=iris_df)
    iris_df = iris_df.loc[(iris_df['species'] == species_lp[s0]) | (iris_df['species'] == species_lp[s1])]

    # Load as R^4 observations and map labels to binary values
    X = np.array(iris_df)[:, :-1]
    Y = np.array(iris_df)[:, -1]
    unique = set(np.array(Y, dtype=np.int))
    Y[Y == min(unique)] = -1
    Y[Y == max(unique)] = 1

    # Partition the data
    iris_data = PartitionData()
    iris_data.add_data(X, Y)
    iris_data.partition(ratio)

    # Return data object and map to species name
    indices_lp = {species_lp[k] : k for k in species_lp}
    return iris_data, {1: indices_lp[max(unique)], -1: indices_lp[min(unique)]}


def classify_species(s0, s1, ratio):
    sc = SimpleClassifier()
    iris, names = load_iris_data(s0, s1, ratio)
    sc.add_data(*iris.training)
    sc.train()
    valr = np.array([sc.classify(x) for x in iris.training[0]])
    val = np.array([sc.classify(x) for x in iris.validation[0]])
    return risk(valr, iris.training[1]), risk(val, iris.validation[1])

def problem_1():
    risk1, error1 = classify_species('I. setosa', 'I. versicolor', 0.8)
    risk2, error2 = classify_species('I. virginica', 'I. versicolor', 0.8)
    risk3, error3 = classify_species('I. virginica', 'I. setosa', 0.8)
    print('Empirical Risk of I. Setosa, I. Versicolor: ', risk1)
    print('Empirical Risk of I. Virginica, I. Versicolor: ', risk2)
    print('Empirical Risk of I. Virginica, I. Setosa: ', risk3)
    print('Classification Error of I. Setosa and I. Versicolor: ', error1)
    print('Classification Error of I. Virginica and I. Versicolor: ', error2)
    print('Classification Error of I. Virginica and I. Setosa: ', error3)

def problem_2():
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from sklearn import datasets
    from sklearn.decomposition import PCA
    from matplotlib import rc
    rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
    rc('text', usetex=True)

    iris = datasets.load_iris()
    X = iris.data[:, :2]  # we only take the first two features.
    y = iris.target

    x_min, x_max = X[:, 0].min() - .5, X[:, 0].max() + .5
    y_min, y_max = X[:, 1].min() - .5, X[:, 1].max() + .5

    plt.figure(2, figsize=(8 / 1.5, 6 / 1.5))
    plt.clf()

    # Plot the training points
    plt.scatter(X[y == 0, 0], X[y == 0, 1], c='r', cmap=plt.cm.Set1,
                edgecolor='k')
    plt.scatter(X[y == 1, 0], X[y == 1, 1], c='g', cmap=plt.cm.Set1,
                edgecolor='k')
    plt.scatter(X[y == 2, 0], X[y == 2, 1], c='b', cmap=plt.cm.Set1,
                edgecolor='k')
    plt.xlabel('Sepal Length (cm)')
    plt.ylabel('Sepal Width (cm)')

    plt.legend(["I. Setosa", "I. Versicolor", "I. Virginica"])

    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.xticks(())
    plt.yticks(())

    fig = plt.figure(1, figsize=(8, 6))
    ax = Axes3D(fig, elev=-150, azim=110)
    X_reduced = PCA(n_components=3).fit_transform(iris.data)
    ax.scatter(X_reduced[:, 0], X_reduced[:, 1], X_reduced[:, 2], c=y,
               cmap=plt.cm.Set1, edgecolor='k', s=40)
    ax.set_title("First three PCA directions")
    ax.set_xlabel("1st eigenvector")
    ax.w_xaxis.set_ticklabels([])
    ax.set_ylabel("2nd eigenvector")
    ax.w_yaxis.set_ticklabels([])
    ax.set_zlabel("3rd eigenvector")
    ax.w_zaxis.set_ticklabels([])

    plt.show()

if __name__ == "__main__":
    problem_1()
    problem_2()