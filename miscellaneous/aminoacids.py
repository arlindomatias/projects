#!/usr/bin/env python3

import matplotlib
matplotlib.use("TkAgg")

from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

dados_aminoacidos = {
    "Alanina": {"abrev3": "Ala", "abrev1": "A", "smiles": "CC(C(=O)O)N",
        "massa_molecular": 89.09, "solubilidade_agua": 167.2, "logP": -2.85,
        "pKa": {"COOH": 2.34, "NH3+": 9.69}, "pI": 6.00,
        "codons": ["GCU", "GCC", "GCA", "GCG"],
        "classe": "Alifático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Arginina": {"abrev3": "Arg", "abrev1": "R", "smiles": "NC(CCCNC(=N)N)C(=O)O",
        "massa_molecular": 174.20, "solubilidade_agua": 182.6, "logP": -4.20,
        "pKa": {"COOH": 2.17, "NH3+": 9.04, "cadeia lateral": 12.48}, "pI": 10.76,
        "codons": ["CGU", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "classe": "Cátion fixo", "polaridade": "Polar básico", "carga_pH74": "Positiva"},
    "Asparagina": {"abrev3": "Asn", "abrev1": "N", "smiles": "NC(CC(=O)N)C(=O)O",
        "massa_molecular": 132.12, "solubilidade_agua": 25.5, "logP": -3.82,
        "pKa": {"COOH": 2.02, "NH3+": 8.80}, "pI": 5.41,
        "codons": ["AAU", "AAC"],
        "classe": "Amida", "polaridade": "Polar", "carga_pH74": "Neutra"},
    "Aspartato": {"abrev3": "Asp", "abrev1": "D", "smiles": "NC(CC(=O)O)C(=O)O",
        "massa_molecular": 133.10, "solubilidade_agua": 5.0, "logP": -3.89,
        "pKa": {"COOH": 1.88, "NH3+": 9.60, "cadeia lateral": 3.65}, "pI": 2.77,
        "codons": ["GAU", "GAC"],
        "classe": "Ânion", "polaridade": "Base de Brønsted", "carga_pH74": "Negativa"},
    "Cisteína": {"abrev3": "Cys", "abrev1": "C", "smiles": "NC(CS)C(=O)O",
        "massa_molecular": 121.16, "solubilidade_agua": 250.0, "logP": -2.49,
        "pKa": {"COOH": 1.96, "NH3+": 10.28, "cadeia lateral": 8.18}, "pI": 5.07,
        "codons": ["UGU", "UGC"],
        "classe": "Tiol", "polaridade": "Ácido de Brønsted", "carga_pH74": "Neutra"},
    "Fenilalanina": {"abrev3": "Phe", "abrev1": "F", "smiles": "NC(Cc1ccccc1)C(=O)O",
        "massa_molecular": 165.19, "solubilidade_agua": 27.9, "logP": -1.52,
        "pKa": {"COOH": 1.83, "NH3+": 9.13}, "pI": 5.48,
        "codons": ["UUU", "UUC"],
        "classe": "Aromático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Glicina": {"abrev3": "Gly", "abrev1": "G", "smiles": "NCC(=O)O",
        "massa_molecular": 75.07, "solubilidade_agua": 249.9, "logP": -3.21,
        "pKa": {"COOH": 2.34, "NH3+": 9.60}, "pI": 5.97,
        "codons": ["GGU", "GGC", "GGA", "GGG"],
        "classe": "Alifático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Glutamato": {"abrev3": "Glu", "abrev1": "E", "smiles": "NC(CCC(=O)O)C(=O)O",
        "massa_molecular": 147.13, "solubilidade_agua": 8.6, "logP": -3.69,
        "pKa": {"COOH": 2.19, "NH3+": 9.67, "cadeia lateral": 4.25}, "pI": 3.22,
        "codons": ["GAA", "GAG"],
        "classe": "Ânion", "polaridade": "Base de Brønsted", "carga_pH74": "Negativa"},
    "Glutamina": {"abrev3": "Gln", "abrev1": "Q", "smiles": "NC(CCC(=O)N)C(=O)O",
        "massa_molecular": 146.15, "solubilidade_agua": 25.0, "logP": -3.64,
        "pKa": {"COOH": 2.17, "NH3+": 9.13}, "pI": 5.65,
        "codons": ["CAA", "CAG"],
        "classe": "Amida", "polaridade": "Polar", "carga_pH74": "Neutra"},
    "Histidina": {"abrev3": "His", "abrev1": "H", "smiles": "NC(Cc1c[nH]cn1)C(=O)O",
        "massa_molecular": 155.16, "solubilidade_agua": 41.9, "logP": -3.52,
        "pKa": {"COOH": 1.82, "NH3+": 9.17, "cadeia lateral": 6.00}, "pI": 7.59,
        "codons": ["CAU", "CAC"],
        "classe": "Catiônico", "polaridade": "Ácido e base de Brønsted",
        "carga_pH74": "Positiva (10%) / Neutra (90%)"},
    "Isoleucina": {"abrev3": "Ile", "abrev1": "I", "smiles": "CCC(C)C(N)C(=O)O",
        "massa_molecular": 131.17, "solubilidade_agua": 34.5, "logP": -1.70,
        "pKa": {"COOH": 2.36, "NH3+": 9.60}, "pI": 6.02,
        "codons": ["AUU", "AUC", "AUA"],
        "classe": "Alifático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Leucina": {"abrev3": "Leu", "abrev1": "L", "smiles": "CC(C)CC(N)C(=O)O",
        "massa_molecular": 131.17, "solubilidade_agua": 21.5, "logP": -1.52,
        "pKa": {"COOH": 2.36, "NH3+": 9.60}, "pI": 5.98,
        "codons": ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
        "classe": "Alifático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Lisina": {"abrev3": "Lys", "abrev1": "K", "smiles": "NCCCCC(N)C(=O)O",
        "massa_molecular": 146.19, "solubilidade_agua": 739.0, "logP": -3.05,
        "pKa": {"COOH": 2.18, "NH3+": 8.95, "cadeia lateral": 10.53}, "pI": 9.74,
        "codons": ["AAA", "AAG"],
        "classe": "Cátion", "polaridade": "Ácido de Brønsted", "carga_pH74": "Positiva"},
    "Metionina": {"abrev3": "Met", "abrev1": "M", "smiles": "CSCCC(N)C(=O)O",
        "massa_molecular": 149.21, "solubilidade_agua": 33.8, "logP": -1.87,
        "pKa": {"COOH": 2.28, "NH3+": 9.21}, "pI": 5.74,
        "codons": ["AUG"],
        "classe": "Tioéter", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Prolina": {"abrev3": "Pro", "abrev1": "P", "smiles": "OC(=O)C1CCCN1",
        "massa_molecular": 115.13, "solubilidade_agua": 1620.0, "logP": -2.54,
        "pKa": {"COOH": 1.99, "NH2+": 10.60}, "pI": 6.30,
        "codons": ["CCU", "CCC", "CCA", "CCG"],
        "classe": "Cíclico", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Serina": {"abrev3": "Ser", "abrev1": "S", "smiles": "NC(CO)C(=O)O",
        "massa_molecular": 105.09, "solubilidade_agua": 250.0, "logP": -3.29,
        "pKa": {"COOH": 2.21, "NH3+": 9.15}, "pI": 5.68,
        "codons": ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
        "classe": "Hidroxílico", "polaridade": "Polar", "carga_pH74": "Neutra"},
    "Tirosina": {"abrev3": "Tyr", "abrev1": "Y", "smiles": "NC(Cc1ccc(O)cc1)C(=O)O",
        "massa_molecular": 181.19, "solubilidade_agua": 0.45, "logP": -2.26,
        "pKa": {"COOH": 2.20, "NH3+": 9.11, "cadeia lateral": 10.07}, "pI": 5.66,
        "codons": ["UAU", "UAC"],
        "classe": "Aromático", "polaridade": "Ácido de Brønsted", "carga_pH74": "Neutra"},
    "Treonina": {"abrev3": "Thr", "abrev1": "T", "smiles": "CC(O)C(N)C(=O)O",
        "massa_molecular": 119.12, "solubilidade_agua": 98.1, "logP": -2.94,
        "pKa": {"COOH": 2.11, "NH3+": 9.62}, "pI": 5.60,
        "codons": ["ACU", "ACC", "ACA", "ACG"],
        "classe": "Hidroxílico", "polaridade": "Polar", "carga_pH74": "Neutra"},
    "Triptofano": {"abrev3": "Trp", "abrev1": "W", "smiles": "NC(Cc1c[nH]c2ccccc12)C(=O)O",
        "massa_molecular": 204.23, "solubilidade_agua": 13.6, "logP": -1.06,
        "pKa": {"COOH": 2.38, "NH3+": 9.39}, "pI": 5.89,
        "codons": ["UGG"],
        "classe": "Aromático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
    "Valina": {"abrev3": "Val", "abrev1": "V", "smiles": "CC(C)C(N)C(=O)O",
        "massa_molecular": 117.15, "solubilidade_agua": 58.1, "logP": -2.26,
        "pKa": {"COOH": 2.32, "NH3+": 9.62}, "pI": 5.96,
        "codons": ["GUU", "GUC", "GUA", "GUG"],
        "classe": "Alifático", "polaridade": "Apolar", "carga_pH74": "Neutra"},
}

def mostrar_estrutura(nome, smiles):
    mol = Chem.MolFromSmiles(smiles)
    img = Draw.MolToImage(mol, size=(400, 400))
    plt.figure()
    plt.imshow(img)
    plt.axis("off")
    plt.title(nome)
    plt.show()

def mostrar_caracteristicas(nome, dados):
    print(f"\n--- {nome} ---")
    print(f"Abreviação (3 letras): {dados['abrev3']}")
    print(f"Abreviação (1 letra): {dados['abrev1']}")
    print(f"Classe: {dados['classe']}")
    print(f"Polaridade química: {dados['polaridade']}")
    print(f"Net charge (pH 7.4): {dados['carga_pH74']}")
    print(f"Massa molecular: {dados['massa_molecular']} g/mol")
    print(f"Solubilidade em água (25°C): {dados['solubilidade_agua']} g/L")
    print(f"LogP: {dados['logP']}")
    print("pKa dos hidrogênios ionizáveis:")
    for grupo, valor in dados["pKa"].items():
        print(f"  {grupo}: {valor}")
    print(f"Ponto isoelétrico (pI): {dados['pI']}")
    print(f"Códons: {', '.join(dados['codons'])}")

lista_ordenada = sorted(dados_aminoacidos.keys())

while True:
    print("\nAminoácidos disponíveis (ordem alfabética):")
    for i, aa in enumerate(lista_ordenada, start=1):
        print(f"{i} - {aa}")
    print("0 - Sair")

    escolha = input("\nDigite o número correspondente ao aminoácido desejado: ").strip()

    if escolha == "0":
        print("Encerrando o programa.")
        break

    if escolha.isdigit() and 1 <= int(escolha) <= len(lista_ordenada):
        nome_escolhido = lista_ordenada[int(escolha) - 1]
        mostrar_caracteristicas(nome_escolhido, dados_aminoacidos[nome_escolhido])
        mostrar_estrutura(nome_escolhido, dados_aminoacidos[nome_escolhido]["smiles"])
    else:
        print("Número inválido. Tente novamente.")