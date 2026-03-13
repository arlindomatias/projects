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

## 2. Lista de Amostras

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

## 3. Download dos Dados

### Testar com uma amostra

> Definir a variável `SRR` com o ID da amostra de teste. Todos os comandos avulsos abaixo usam essa variável.

```bash
SRR="SRR33247676"

prefetch $SRR                                      # Baixa o arquivo .sra
fasterq-dump $SRR --split-files --threads 8        # Converte para FASTQ
```

### Download em loop

```bash
while read SRR; do
  echo "Baixando $SRR"
  prefetch $SRR
  fasterq-dump $SRR --split-files --threads 8
done < srr_list.txt
```

### Alternativa: download manual via EBI

```bash
# SRR deve estar definido (ver acima)
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/476/${SRR}/${SRR}_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/476/${SRR}/${SRR}_2.fastq.gz
gunzip *.gz
```

> Para listar o conteúdo do diretório FTP:
> ```bash
> wget --spider -r ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/476/${SRR}/
> ```

---

## 4. Download dos Metadados

```bash
curl -L \
  "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1254109&result=read_run&fields=run_accession,sample_accession,experiment_accession,study_accession,library_layout,library_strategy,instrument_model,scientific_name,sample_title,sample_description&format=tsv" \
  -o Metadados.tsv
```

---

## 5. Organização dos Diretórios

```bash
mkdir -p fastq salmon_quant
mv SRR*_*.fastq fastq/
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

## 6. Controle de Qualidade

```bash
fastqc fastq/*.fastq -o fastq/
```

---

## 7. Transcriptoma de Referência (Mouse GRCm39 / mm39)

### Para o Salmon

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.transcripts.fa.gz
gunzip gencode.vM35.transcripts.fa.gz
head gencode.vM35.transcripts.fa
```

### Para o R (tximport)

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M35/gencode.vM35.annotation.gtf.gz
gunzip gencode.vM35.annotation.gtf.gz
```

---

## 8. Índice Salmon

> Executar **uma vez** por dataset/referência.

```bash
salmon index \
  -t gencode.vM35.transcripts.fa \
  -i salmon_mm39_index \
  -k 31 \
  --gencode
```

---

## 9. Quantificação com Salmon

### Testar com uma amostra

```bash
# SRR deve estar definido (ver seção 3)
salmon quant \
  -i salmon_mm39_index \
  -l A \
  -1 fastq/${SRR}_1.fastq \
  -2 fastq/${SRR}_2.fastq \
  -p 8 \
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
- `Detected lib type`
- `Mapping rate`
- `Warnings`

> ⚠️ Se a taxa de mapeamento for **< 30%**, pare. Algo está errado.

### Quantificação em loop

```bash
for r1 in fastq/*_1.fastq; do
  base=$(basename $r1 _1.fastq)
  r2=fastq/${base}_2.fastq

  salmon quant \
    -i salmon_mm39_index \
    -l A \
    -1 $r1 \
    -2 $r2 \
    -p 8 \
    --validateMappings \
    --gcBias \
    --seqBias \
    -o salmon_quant/$base
done
```

---

## 10. Verificação dos Resultados

### Estrutura esperada

```
salmon_quant/
├── SRR33247676/
│   ├── quant.sf
│   └── logs/
├── SRR33247677/
│   ├── quant.sf
│   └── logs/
└── SRR33247678/
    ├── quant.sf
    ├── aux_info/
    ├── libParams/
    └── logs/
```

### Checar o quant.sf

```bash
# SRR deve estar definido (ver seção 3)
head salmon_quant/${SRR}/quant.sf
```

Deve conter colunas: `Name`, `Length`, `EffectiveLength`, `TPM`, `NumReads` — sem valores `NaN`.

### Sanidade dos TPMs

```bash
# Soma dos TPMs — deve dar ~1.000.000
awk '{s+=$4} END {print s}' salmon_quant/${SRR}/quant.sf

# Total de reads mapeados — deve ser dezenas de milhões
awk '{s+=$5} END {print s}' salmon_quant/${SRR}/quant.sf
```

---

## 11. Gestão de Espaço em Disco

> Caso o espaço seja limitado, é possível baixar e processar em partes. O Salmon **não precisa dos FASTQs** após a quantificação. Tendo o `quant.sf`, pode-se apagar os FASTQs e arquivos `.sra`.

Manter obrigatoriamente:
- `salmon_quant/SRRxxxx/quant.sf`
- Metadados (`Metadados.tsv` / `sample_table`)
- Índice (`salmon_mm39_index/`)

---

## 12. Decisão Antes de Escalar

Após o piloto com 2–3 amostras, avaliar:
- Tempo por SRR
- Uso de disco
- Taxa de mapeamento
- Variabilidade entre amostras

---

## 13. Próximos Passos no R

Com 2–3 SRRs quantificados corretamente:

1. Importar com `tximport` (Bioconductor)
2. Verificar correlação entre amostras
3. Testar PCA
4. Validar estabilidade dos números entre amostras
5. Decidir nível de análise: **transcript** vs **gene**
6. Só então seguir para análise de **expressão diferencial (DE)**
