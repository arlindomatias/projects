import pandas as pd
from urllib.request import Request, urlopen
from bs4 import BeautifulSoup
from datetime import datetime
import dtale
import requests
import numpy as np
import matplotlib.pyplot as plt
import pandas_datareader.data as web
import yfinance as yf

### Criar tabela geral de indicadores com dados do fundamentus

# URL principal
url_total = 'https://fundamentus.com.br/resultado.php'
req = Request(url_total, headers={'User-Agent': 'Mozilla/5.0'})
html_total = urlopen(req)
soup_total = BeautifulSoup(html_total, 'html.parser')

# Extrair tabela principal
tabela = soup_total.find('table', {'id': 'resultado'})
linhas = tabela.find('tbody').find_all('tr')
nomes_colunas = [th.get_text(strip=True) for th in tabela.find('thead').find_all('th')]

# Coletar dados básicos
todos_papeis = []
for linha in linhas:
    colunas = linha.find_all('td')
    dados = [c.get_text(strip=True) for c in colunas]
    todos_papeis.append(dados)

df = pd.DataFrame(todos_papeis, columns=nomes_colunas)

# Converter strings para float
def str_para_float(x):
    if isinstance(x, str):
        x = x.replace('.', '').replace(',', '.').replace('%', '')
        try:
            return float(x)
        except:
            return x
    return x

for col in df.columns[1:]:
    df[col] = df[col].apply(str_para_float)

# --- Novo trecho: criar segundo DataFrame com filtro ---
df = df[(df['P/L'] > 0) & (df['EV/EBITDA'] > 0)]

# Criar colunas de ranking
df['classificação EV/EBITDA'] = df['EV/EBITDA'].rank(method='min', ascending=True).astype(int)
df['classificação ROE'] = df['ROE'].rank(method='min', ascending=False).astype(int)
df['ranking combinado'] = df['classificação EV/EBITDA'] + df['classificação ROE']
df = df.sort_values('ranking combinado')

# Coloca o papel como índice
df.set_index('Papel', inplace=True)

# Abrir como documento html
dt = dtale.show(df)
dt.open_browser()

# Exportar para excel
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
out_xlsx = f"fundamentus_{timestamp}.xlsx"

### Obter dados do Yahoo Finanças
plt.style.use('seaborn-v0_8-darkgrid')
papel = 'SAPR4' # Desempenho de papel específico, para desempenho do IBOVESPA, usar ^BVSP sem SA
histórico_papel = yf.download(papel+'.SA', start = '2000-01-01')

plt.figure(figsize=(22, 8))
ax = plt.gca()  # obtém o eixo atual

# Fechamento
ax.plot(histórico_papel.index, histórico_papel['Close'], label=f'{papel} (Fechamento)', linewidth=1.5)

# médias móveis
ax.plot(histórico_papel.index, histórico_papel['Close'].rolling(21).mean(), label='MM21', linestyle='--', linewidth=1.2)
ax.plot(histórico_papel.index, histórico_papel['Close'].rolling(200).mean(), label='MM200', linestyle='--', linewidth=1.2)


plt.legend()
plt.title(f'Desempenho de {papel}')
plt.xlabel('Data')
plt.ylabel('Preço (R$)')
plt.style.use('seaborn-v0_8-darkgrid')
plt.show()

