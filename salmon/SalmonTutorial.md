# Pipeline RNA-seq com Salmon
**BioProject:** [PRJNA1254109](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1254109&o=acc_s%3Aa) — GSE295281  
**Referência:** [Documentação Salmon](https://salmon.readthedocs.io/en/latest/)

---

## 1. Configuração do Ambiente

```bash
conda create -n salmon -c bioconda salmon
conda activate salmon

conda install -c bioconda -c conda-forge \
    sra-tools \
    ca-certificates \
    certifi \
    openssl

conda install -c bioconda fastqc
vdb-config --interactive
```

---

## 2. Variáveis Globais

> Definir aqui uma única vez. Todos os comandos avulsos do tutorial usam essas variáveis.

```bash
SRR="SRR33247678"   # Amostra de teste — alterar conforme necessário
THREADS=4           # Número de threads disponíveis na máquina
```

---

## 3. Lista de Amostras

Criar o arquivo `srr_list.txt` com os IDs das amostras:

```
SRR33247676
SRR33247677
SRR33247678
SRR33247679
SRR33247680
SRR33247681
SRR33247682
SRR33247683
SRR33247684
SRR33247685
SRR33247686
```

---

## 4. Organização dos Diretórios

```bash
mkdir -p fastq salmon_quant
```

Estrutura esperada do projeto:

```
.
├── Snakefile
├── srr_list.txt
├── envs/
│   └── salmon.yaml
├── fastq/
├── qc/
├── salmon_quant/
└── ref/
```

---

## 5. Download dos Metadados

```bash
curl -L \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1254109&result=read_run&fields=run_accession,sample_accession,experiment_accession,study_accession,library_layout,library_strategy,instrument_model,scientific_name,sample_title,sample_description&format=tsv" \
  -o Metadados.tsv
```

---

## 6. Transcriptoma de Referência (Mouse GRCm39 / mm39)

### Para o Salmon

```bash
wget -P ref/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.transcripts.fa.gz
gunzip ref/gencode.vM35.transcripts.fa.gz
head ref/gencode.vM35.transcripts.fa
```
---

## 7. Índice Salmon

> Executar **uma vez** por dataset/referência.

```bash
salmon index \
  -t ref/gencode.vM35.transcripts.fa \
  -i salmon_mm39_index \
  -k 31 \
  --gencode
```

---

---

# Fase 1: Piloto (1 amostra)

> Objetivo: validar o pipeline de ponta a ponta antes de escalar. Se qualquer etapa falhar ou der resultado suspeito, **pare aqui**.

---

## 8. Download da Amostra de Teste

```bash
prefetch $SRR
fasterq-dump $SRR --split-files --threads $THREADS
mv ${SRR}_*.fastq fastq/
```

### Verificar o download

```bash
ls -lh fastq/${SRR}_*.fastq
wc -l fastq/${SRR}_1.fastq   # Número de linhas ÷ 4 = número de reads
```

> O número de reads deve ser da ordem de dezenas de milhões. Arquivos vazios ou muito pequenos indicam download corrompido.

### Alternativa: download manual via EBI

```bash
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/476/${SRR}/${SRR}_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/476/${SRR}/${SRR}_2.fastq.gz
gunzip fastq/*.gz
```

> Para listar o conteúdo do diretório FTP:
> ```bash
> wget --spider -r ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/476/${SRR}/
> ```

---

## 9. Controle de Qualidade

```bash
fastqc fastq/*.fastq -o fastq/
```

---

## 10. Quantificação (Piloto)

```bash
salmon quant \
  -i salmon_mm39_index \
  -l A \
  -1 fastq/${SRR}_1.fastq \
  -2 fastq/${SRR}_2.fastq \
  -p $THREADS \
  --validateMappings \
  --gcBias \
  --seqBias \
  -o salmon_quant/${SRR}
```

### Verificar o log

```bash
less salmon_quant/${SRR}/logs/salmon_quant.log
```

Checar especialmente:
- `Detected lib type` — deve ser `IU`, `ISF` ou `ISR` para paired-end
- `Mapping rate` — esperado > 70%; abaixo de 30% é sinal de problema grave
- `Warnings` — qualquer aviso sobre biblioteca ou viés

> ⚠️ Se a taxa de mapeamento for **< 30%**, pare. Causas comuns: genoma de referência errado, FASTQs corrompidos ou biblioteca single-end tratada como paired-end.

---

## 11. Verificação dos Resultados (Piloto)

### Estrutura esperada após o piloto

```
salmon_quant/
└── SRR33247678/
    ├── quant.sf
    ├── aux_info/
    ├── libParams/
    └── logs/
```

### Checar o quant.sf

```bash
head salmon_quant/${SRR}/quant.sf
```

Deve conter as colunas `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads` — sem valores `NaN`.

**Sobre os valores de TPM:**
- TPM normaliza pelo comprimento do transcript e pela profundidade de sequenciamento, facilitando comparação entre amostras
- Transcritos com TPM muito alto (> 10.000) frequentemente são genes housekeeping (*Actb*, *Gapdh*) — esperado
- Transcritos com TPM = 0 em quase todas as amostras são candidatos a filtro antes da DE
- Nunca use TPM diretamente para DE — use `NumReads` via `tximport`

### Sanidade dos TPMs

```bash
# Soma dos TPMs — deve dar ~1.000.000 (por definição da métrica)
awk 'NR>1 {s+=$4} END {print s}' salmon_quant/${SRR}/quant.sf

# Total de reads mapeados — deve ser dezenas de milhões
awk 'NR>1 {s+=$5} END {print s}' salmon_quant/${SRR}/quant.sf
```

---

## 12. Limpeza Pós-Piloto (opcional)

> Se o espaço em disco for limitado, os FASTQs e arquivos `.sra` podem ser removidos após a quantificação. O Salmon não os utiliza nas etapas seguintes.

```bash
rm fastq/${SRR}_*.fastq
rm -rf ${SRR}/   # diretório .sra gerado pelo prefetch
```

Manter obrigatoriamente:
- `salmon_quant/${SRR}/quant.sf`
- `Metadados.tsv`
- `salmon_mm39_index/`

---

## 13. Decisão: Escalar ou Não?

Antes de baixar todas as amostras, avaliar com base no piloto:

| Critério | Esperado | Ação se falhar |
|---|---|---|
| Mapping rate | > 70% | Revisar referência ou parâmetros |
| Soma dos TPMs | ~1.000.000 | Verificar integridade do quant.sf |
| Total de reads | Dezenas de milhões | Verificar download |
| Tempo por SRR | Estimativa para escala | Ajustar paralelismo |
| Espaço em disco | Suficiente para todas | Planejar limpeza por lotes |

---

---

# Fase 2: Produção (todas as amostras)

> Só iniciar após validar o piloto com sucesso.

---

## 14. Download em Loop

```bash
while read SRR; do
  echo "=== Baixando $SRR ==="
  prefetch $SRR
  fasterq-dump $SRR --split-files --threads $THREADS
  mv ${SRR}_*.fastq fastq/
done < srr_list.txt
```

---

## 15. Quantificação em Loop

```bash
for r1 in fastq/*_1.fastq; do
  base=$(basename $r1 _1.fastq)
  r2=fastq/${base}_2.fastq

  salmon quant \
    -i salmon_mm39_index \
    -l A \
    -1 $r1 \
    -2 $r2 \
    -p $THREADS \
    --validateMappings \
    --gcBias \
    --seqBias \
    -o salmon_quant/$base
done
```

---

## 16. Estrutura Final Esperada

```
salmon_quant/
├── SRR33247676/
│   ├── quant.sf
│   └── logs/
├── SRR33247677/
│   ├── quant.sf
│   └── logs/
└── ...
```

---

## 17. Limpeza em Lote (opcional)

```bash
# Remover todos os FASTQs após quantificação completa
rm fastq/*.fastq

# Remover diretórios .sra
rm -rf ~/ncbi/public/sra/SRR332477*/
```

---

## 18. Troubleshooting

| Sintoma | Causa provável | Solução |
|---|---|---|
| Mapping rate < 30% | Referência errada (ex: humano em vez de camundongo) | Verificar e recriar o índice com o genoma correto |
| Mapping rate < 30% | Biblioteca single-end tratada como paired-end | Usar `-l U` ou `-l A` e revisar o design experimental |
| `NaN` no quant.sf | FASTQ corrompido ou truncado | Refazer o download |
| Soma TPM ≠ ~1.000.000 | Linha de header contada no awk | Adicionar `NR>1` no awk (já corrigido acima) |
| `fasterq-dump` lento | Falta de cache local configurado | Rodar `vdb-config --interactive` e definir cache |
| Índice não encontrado | Caminho relativo errado | Usar caminho absoluto com `$(pwd)/salmon_mm39_index` |

---

---

# Fase 3: Análise no R

---

## 19. Próximos Passos

Com todas as amostras quantificadas:

1. Importar com `tximport` (Bioconductor) — nível gene ou transcript
2. Verificar correlação entre amostras (heatmap)
3. Testar PCA — amostras do mesmo grupo devem agrupar
4. Validar estabilidade dos números entre réplicas
5. Só então seguir para análise de **expressão diferencial (DE)** com `DESeq2` ou `edgeR`
