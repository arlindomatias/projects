import pandas as pd
from sklearn.tree import DecisionTreeClassifier
from dtreeviz import model

# Carrega dados
filename = 'https://github.com/lmassaron/datasets/releases/download/1.0/tennis.feather'
tennis = pd.read_feather(filename)

# Features e target
X = tennis[['outlook', 'temperature', 'humidity', 'wind']]
X = pd.get_dummies(X)  # transforma categóricas em numéricas
y = tennis.play

# Treina Decision Tree
dt = DecisionTreeClassifier()
dt.fit(X, y)

# Cria visualização com a nova API
m = model(
    dt,
    X,
    y,
    target_name='play_tennis',
    feature_names=X.columns,
    class_names=["No", "Yes"])

# Abre gráfico no visualizador externo (browser ou visualizador padrão do sistema)
m.view()
