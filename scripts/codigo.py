'''
Entrega un script en Python que implemente de forma detallada un análisis
funcional de los genes COX4I2, ND1 y ATP6. El análisis debe utilizar librerías
de Python para evaluar los procesos biológicos asociados a estos genes.
Asegúrate de que tu código esté bien documentado y describa claramente
los métodos y bases de datos utilizadas para obtener información funcional.
'''

#######################################
## 0. Librerias
#######################################
import mygene
import pandas as pd
import gseapy as gp
import requests

#######################################
## 1. Entrada de datos
#######################################

# Se podria hacer una funcion como la siguiente si tuvieramos los genes en un fichero
def importar_genes(fichero):
    lista_genes = []
    with open(fichero, 'r') as f:
        for l in f:
            lista_genes.append(l.strip())
    return lista_genes

genes = ['COX4I2', 'ND1', 'ATP6']
especie = "human"

# Algunos IDs de genes hay que estandarizarlos
# Se podra hacer una funcion automatica?
def normalize_mt_symbols_hs(genes: list[str]) -> list[str]:
    """Normaliza símbolos mitocondriales comunes en humano."""
    mapping = {"ND1": "MT-ND1", "ATP6": "MT-ATP6"}
    return [mapping.get(g.upper(), g) for g in genes]

genes_norm = normalize_mt_symbols_hs(genes)

#######################################
#  2. Mapeo y anotacion de genes #
#######################################
def map_with_mygene(genes, species) -> pd.DataFrame:
    """
    Mapea símbolos a IDs (Entrez/Ensembl) y obtiene un resumen/aliases.
    species: 'human' (9606), 'mouse' (10090), etc. MyGene acepta nombre/taxid.
    """
    mg = mygene.MyGeneInfo()
    fields = "symbol,name,entrezgene,ensembl.gene,summary,alias,taxid"
    df = mg.querymany(
        genes,
        scopes="symbol,alias,ensemblgene,entrezgene",
        fields=fields,
        species=species,
        as_dataframe=True,
        verbose=False,
    )

    # Ya tenemos el mapeo, hay que aplanar los datos
    # Renombramos etiquetas a nombres más consistentes
    df = df.rename(columns={
        "symbol": "symbol",
        "name": "gene_name",
        "entrezgene": "entrez_id",
        "ensembl.gene": "ensembl_id",
        "summary": "summary",
        "alias": "aliases",
        "taxid": "taxid",
    })

    # Aplanar listas/dicts en alias/ensembl
    def flatten(x):
        if isinstance(x, list):
            return ";".join(sorted(set(map(str, x))))
        if isinstance(x, dict) and "gene" in x:
            return x["gene"]
        return x

    df["ensembl_id"] = df["ensembl_id"].map(flatten)
    df["aliases"] = df["aliases"].map(flatten)

    # Orden estético
    if "symbol" in df.columns:
        df = df.sort_values("symbol")

    # Lista de símbolos oficiales únicos para usar en enriquecimiento
    normalized_symbols = df["symbol"].dropna().unique().tolist()

    return df, normalized_symbols

datos, ids = map_with_mygene(genes_norm, especie)