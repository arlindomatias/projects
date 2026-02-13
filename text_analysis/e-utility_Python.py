#!/usr/bin/env python3
import os
import subprocess
import csv
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez


# Configurações
Entrez.email = "arlindo.matias@usp.br"  # obrigatório pelo NCBI
TERM = "chemobrain"
MAX = 20
OUTDIR = "./pubmed_abstracts"
os.makedirs(OUTDIR, exist_ok=True)
CSV_FILE = os.path.join(OUTDIR, "pubmed_data.csv")

# Função para executar comando shell
def run_cmd(cmd):
    result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Erro ao executar: {cmd}\n{result.stderr}")
    return result.stdout.strip()

# Limpa CSV anterior, se existir
with open(csv_file, "w", newline='', encoding='utf-8') as f:
    writer = csv.writer(f)
    writer.writerow(["ID", "Title", "Authors", "PubDate", "AbstractFile"])

# Buscar IDs
print(f"Buscando IDs para o termo: {TERM}")
with Entrez.esearch(db="pubmed", term=TERM, retmax=MAX) as handle:
    search_results = Entrez.read(handle)
ids_list = search_results["IdList"]

print(f"Encontrados {len(ids_list)} artigos.")

# Função para processar cada ID
def process_id(pubmed_id):
    # Arquivo do abstract
    abstract_file = os.path.join(OUTDIR, f"abstract_{pubmed_id}.txt")
    
    # Baixar resumo
    with Entrez.efetch(db="pubmed", id=pubmed_id, rettype="abstract", retmode="text") as handle:
        abstract_text = handle.read()
    with open(abstract_file, "w", encoding="utf-8") as f:
        f.write(abstract_text)
    
    # Obter título, autores e data
    with Entrez.esummary(db="pubmed", id=pubmed_id, retmode="xml") as handle:
        summary = Entrez.read(handle)[0]
    title = summary.get("Title", "")
    authors = ", ".join([a["Name"] for a in summary.get("AuthorList", [])])
    pubdate = summary.get("PubDate", "")
    
    # Salvar no CSV
    with open(CSV_FILE, "a", newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([pubmed_id, title, authors, pubdate, abstract_file])
    
    print(f"Processado: {pubmed_id}")


# Rodar com múltiplas threads
with ThreadPoolExecutor(max_workers=5) as executor:
    executor.map(process_id, ids_list)