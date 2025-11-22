#!/usr/bin/env python3
import os
import subprocess
import csv
from concurrent.futures import ThreadPoolExecutor

# Termo de busca
TERM = "chemobrain"
# Número máximo de artigos
MAX = 20
# Diretório para salvar arquivos
OUTDIR = "./pubmed_abstracts"
os.makedirs(OUTDIR, exist_ok=True)

# Função para executar comando shell
def run_cmd(cmd):
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Erro ao executar: {cmd}\n{result.stderr}")
    return result.stdout.strip()

# 1️⃣ Buscar IDs
print(f"Buscando IDs para o termo: {TERM}")
ids_output = run_cmd(f'esearch -db pubmed -query "{TERM}" -retmax {MAX} | efetch -format uid')
ids_list = ids_output.splitlines()

# CSV output
csv_file = os.path.join(OUTDIR, "pubmed_data.csv")

# Limpa CSV anterior, se existir
with open(csv_file, "w", newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(["ID", "Title", "Authors", "PubDate", "AbstractFile"])

# Função para baixar e processar cada ID
def process_id(ID):
    abstract_file = os.path.join(OUTDIR, f"abstract_{ID}.txt")
    # Baixar abstract
    run_cmd(f'efetch -db pubmed -id {ID} -format abstract > "{abstract_file}"')
    
    # Obter título, autores e data
    summary = run_cmd(f'esummary -_
