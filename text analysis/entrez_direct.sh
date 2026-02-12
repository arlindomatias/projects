# Instalarção
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"


# Adicionar ao PATH
echo "export PATH=/home/codespace/edirect:\${PATH}" >> ${HOME}/.bashrc
export PATH=${HOME}/edirect:${PATH}
source ~/.bashrc

# esearch performs a new Entrez search using terms in indexed fields.
# elink looks up neighbors (within a database) or links (between databases).
# efilter filters or restricts the results of a previous query.
# efetch downloads records or reports in a designated format.
# xtract converts XML into a table of data values.
# einfo obtains information on indexed fields in an Entrez database.
# epost uploads unique identifiers (UIDs) or sequence accession numbers.
# nquire sends a URL request to a web page or CGI service.

# Para ver ajuda:
esearch --help

# Query
esearch -db pubmed -query "chemobrain"

# Usando pipe
esearch -db pubmed -query "chemobrain" | elink -related

# Query com múltiplas linhas usando \
esearch -db pubmed -query "opsin gene conversion" | \
elink -related | \
elink -target protein

# Selecionando o formato com efetch
esearch -db pubmed -query "chemobrain" | \
efetch -format medline

# Filtrar os resultados por tempo
esearch -db pubmed -query "chemobrain" | \ 
efilter -days 60 -datetype PDAT| \ 
efetch -format medline

# Sequencias de biomoléculas
esearch -db protein -query "lycopene cyclase" | \
efetch -format fasta | \
grep '.' # Essa linha serve para remover os espaços em branco


[AFFL] Affiliation [FILT] Filter [MESH] MeSH Terms
[ALL] All Fields [JOUR] Journal [PTYP] Publication Type
[AUTH] Author [LANG] Language [WORD] Text Word
[FAUT] Author - First [MAJR] MeSH Major Topic [TITL] Title
[LAUT] Author - Last [SUBH] MeSH Subheading [TIAB] Title/Abstract
[PDAT] Date - Publication [UID] UID


WordAtATime() {
sed 's/[^a-zA-Z0-9]/ /g' |
tr 'A-Z' 'a-z' |
xargs -n 1
}
alias word-at-a-time='WordAtATime'

SortUniqCountRank() {
sort -f |
uniq -i -c |
perl -pe 's/\s*(\d+)\s(.+)/$1\t$2/' |
sort -t $'\t' -k 1,1nr -k 2f
}
alias sort-uniq-count-rank='SortUniqCountRank'