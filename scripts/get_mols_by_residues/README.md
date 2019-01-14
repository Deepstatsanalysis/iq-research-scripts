# Analise de Relação de Distância entre Molécula e Resíduo

A problemática que essa análise visa resolver é entender como os arquivos exportados pelo programa Gold em .mol2 e uma base de dados de Proteina (.pdb) obtida de uma outra pesquisa se relacionam na intenção de localizar residuos da proteina que reagem com a molécula em uma determinada distancia máxima.

Para executar o script basta configurar um arquivo `.env` em qualquer local da maquina. Esse arquivo precisa ter a seguinte estrutura:

```
INPUT_LIST=/path/to/data/chemscore2.csv
MOL2_FILE=/path/to/data/chemscore_2_solutions.mol2
PDBS_FILE_FOLDER=/path/to/data/pdbs
OUTPUT_LIST=/path/to/data/filter.csv
RESIDUES=RESIDUES,RESIDUES

MAX_DISTANCE=3.0
````

Feito isso basta rodar o seguinte comando: `python script.py path/to/.env`, esse comando gerar um `.csv` no path que foi configurado no `.env`.

### Váriaveis do .env

