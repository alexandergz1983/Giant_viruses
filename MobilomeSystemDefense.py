##################################################################################################################################

#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# CREATED BY: J-Alexander Garcia-Zea
# jagarciaz2@eafit.edu.co
# python 3.8

##################################################################################################################################

import pandas as pd
import re
import sys
from Bio import SeqIO
import subprocess
import csv
import os
import pandas as pd

##################################################################################################################################

# Prefilter
# Mandatory columns structure: sys_id, sys_beg, sys_end, name_of_profiles_in_sys, specie

paninmune = "test.csv" # change with your paninmune data csv 
df_paninmune = pd.read_csv(paninmune, header=0, on_bad_lines='skip') # change  ,sep=";" "\t" for you date
print(df_paninmune)

print("File uploaded successfully.")
df_paninmune["sys_beg_final"] = df_paninmune["sys_beg"].str.split("_").str[-1]
df_paninmune["sys_end_final"] = df_paninmune["sys_end"].str.split("_").str[-1]
#print(df_paninmune)

def run_diamond(query_protein, db_name, e_value=1e-15, query_cover=None, subject_cover=None):
    """
    Ejecuta DIAMOND con los parámetros especificados, incluyendo la cobertura mínima de la consulta
    y del sujeto si es necesario.
    
    Parámetros:
    - query_protein (str): Archivo multifasta de proteínas de consulta.
    - db_name (str): Nombre de la base de datos DIAMOND.
    - e_value (float): Valor de e-value para filtrar resultados (por defecto 1e-15).
    - query_cover (int): Cobertura mínima de la consulta (%).
    - subject_cover (int): Cobertura mínima del sujeto (%).
    """
    print("Running run_diamond...")  # Mensaje antes de ejecutar subprocess.run
    try:
        with open('output.txt', 'a') as f:  # Abre el archivo de salida en modo append
            # Ejecuta el comando de DIAMOND con los nuevos parámetros
            result = subprocess.run(["diamond", "blastp", 
                                     "-q", query_protein, 
                                     "-d", db_name, 
                                     "-o", "output.txt", 
                                     "--evalue", str(e_value), 
                                     "--outfmt", "6", 
                                     "--max-target-seqs", "1",
                                     f"--query-cover", str(query_cover),  # Cobertura mínima de la consulta optional: 80
                                     f"--subject-cover", str(subject_cover)],  # Cobertura mínima del sujeto optional: 80
                                    stdout=f,  # Redirige stdout al archivo de salida
                                    text=True,
                                    stderr=subprocess.PIPE)

        print("Resultado de subprocess.run obtenido...")  # Mensaje tras ejecutar subprocess.run
        print(f"Salida estándar: {result.stdout}")  # Muestra la salida estándar
        print(f"Salida de error: {result.stderr}")  # Muestra la salida de error (si la hay)

    except Exception as e:
        print(f"Ocurrió un error al ejecutar run_diamond: {e}")

# Base de datos y archivo multifasta de consulta
db_multifasta_file = "mobileOG-db_beatrix-1.6.faa"  # Cambia por la ruta correcta
diamond_db_name = "mobileOG_db"

# Crear la base de datos DIAMOND
subprocess.run(["diamond", "makedb", "--in", db_multifasta_file, "--db", diamond_db_name, "--threads", "8"])

# Ruta del archivo multifasta de proteínas de consulta
input_multifasta_file = "giantvirus.faa"  # Cambia por la ruta correcta

# Leer el archivo CSV con Pandas (si es necesario para otros análisis)
input_csv_file2 = "mobileOG-db-beatrix-1.6.csv"
df_mobilome = pd.read_csv(input_csv_file2)

# Ejecutar la función con valores de cobertura personalizados
run_diamond(input_multifasta_file, diamond_db_name, e_value=1e-15, query_cover=None, subject_cover=None)

##################################################################################################################################

# Defines the sys_id function
def sys_id(sequences = list(SeqIO.parse(input_multifasta_file, 'fasta'))):
    subgrupos = {}
    
    for sys_id in df_paninmune["sys_id"].unique():
        subgrupo = df_paninmune[df_paninmune["sys_id"] == sys_id].copy()
        subgrupos[sys_id] = subgrupo
    
    return subgrupos

# Call the sys_id function with your DataFrame df_panimmune
subgrupos = sys_id(df_paninmune)
#print(df_paninmune)

def neighboring_genes(df_paninmune):
    # Create columns for neighboring genes
    df_paninmune["upstream_gene"] = ""
    df_paninmune["downstream_gene"] = ""

    # Iterating over strings in sys_id
    for sys_id, group in df_paninmune.groupby("sys_id"):
        group = group.sort_values(by=["sys_beg_final"])

        upstream_genes = []
        downstream_genes = []

        # Get numeric positions of sys_beg and sys_end as strings
        sys_beg_values = group["sys_beg"].str.extract(r'_(\d+)', expand=False).astype(str)
        sys_end_values = group["sys_end"].str.extract(r'_(\d+)', expand=False).astype(str)

        # Search for backward neighboring genes (upstream)
        for i, row in group.iterrows():
            upstream_genes = []
            sys_beg_value = int(sys_beg_values[i])

            for j, other_row in group.iterrows():
                if i != j:
                    other_sys_end_value = int(sys_end_values[j])

                    # Verify whether the numeric positions are continuous or discontinuous with jumps of up to 5
                    if 0 < (sys_beg_value - other_sys_end_value) <= 5:
                        upstream_genes.append(other_row["name_of_profiles_in_sys"])

            df_paninmune.at[i, "upstream_gene"] = ", ".join(upstream_genes)

        # Buscar genes vecinos hacia adelante (downstream)
        for i, row in group.iterrows():
            downstream_genes = []
            sys_end_value = int(sys_end_values[i])

            for j, other_row in group.iterrows():
                if i != j:
                    other_sys_beg_value = int(sys_beg_values[j])

                    # Verify whether the numeric positions are continuous or discontinuous with jumps of up to 5
                    if 0 < (other_sys_beg_value - sys_end_value) <= 5:
                        downstream_genes.append(other_row["name_of_profiles_in_sys"])

            df_paninmune.at[i, "downstream_gene"] = ", ".join(downstream_genes)

    return df_paninmune

# Call the function and pass df_panimmune as argument
df_paninmune = neighboring_genes(df_paninmune)
#print(df_paninmune)

##################################################################################################################################

def neighboring_genes_systems(df_paninmune, input_multifasta_file, diamond_db_name, df_mobilome):
    sequences = list(SeqIO.parse(input_multifasta_file, 'fasta'))
    df_paninmune["upstream_neighboring_genes"] = ""
    df_paninmune["upstream_neighboring_pos"] = ""
    df_paninmune["downstream_neighboring_genes"] = ""
    df_paninmune["downstream_neighboring_pos"] = ""

    for index, row in df_paninmune.iterrows():
        sys_beg = row['sys_beg']
        sys_end = row['sys_end']

        # Buscar coincidencias para sys_beg (upstream)
        for i, seq in enumerate(sequences):
            if str(sys_beg) in seq.id:
                upstream_indices = range(max(0, i - 5), i)
                upstream_seqs = [sequences[j] for j in upstream_indices]
                SeqIO.write(upstream_seqs, "temp_upstream_island.faa", "fasta")
                
                run_diamond("temp_upstream_island.faa", diamond_db_name)

                with open('output.txt', 'a+') as f:
                    f.seek(0)
                    blast_result = f.readlines()

                # Copiar los resultados a result_blast.txt
                with open('result_blast.txt', 'a') as result_f:
                    result_f.writelines(blast_result)

                # Obtener la segunda columna de BLAST
                upstream_genes_blast = [line.split()[1] for line in blast_result]

                # Buscar en df_mobilome la coincidencia en "mobileOG fasta Header"
                upstream_genes = []
                for gene in upstream_genes_blast:
                    match = df_mobilome[df_mobilome["mobileOG fasta Header"] == gene]
                    if not match.empty:
                        # Si hay coincidencia, obtener el valor de la columna "mobileOG fasta Header"
                        upstream_genes.append(match["mobileOG fasta Header"].values[0])

                upstream_ids = [seq.id for seq in upstream_seqs if any(seq.id in line for line in blast_result)]

                # Guardar el valor en df_paninmune usando los valores de la columna "mobileOG fasta Header"
                df_paninmune.at[index, "upstream_neighboring_genes"] = ",".join(upstream_genes)
                df_paninmune.at[index, "upstream_neighboring_pos"] = ",".join(upstream_ids)

        # Buscar coincidencias para sys_end (downstream)
        for i, seq in enumerate(sequences):
            if str(sys_end) in seq.id:
                downstream_indices = range(i + 1, min(len(sequences), i + 6))
                downstream_seqs = [sequences[j] for j in downstream_indices]
                SeqIO.write(downstream_seqs, "temp_downstream_island.faa", "fasta")
                
                run_diamond("temp_downstream_island.faa", diamond_db_name)

                with open('output.txt', 'a+') as f:
                    f.seek(0)
                    blast_result = f.readlines()

                # Copiar los resultados a result_blast.txt
                with open('result_blast.txt', 'a') as result_f:
                    result_f.writelines(blast_result)

                # Obtener la segunda columna de BLAST
                downstream_genes_blast = [line.split()[1] for line in blast_result]

                # Buscar en df_mobilome la coincidencia en "mobileOG fasta Header"
                downstream_genes = []
                for gene in downstream_genes_blast:
                    match = df_mobilome[df_mobilome["mobileOG fasta Header"] == gene]
                    if not match.empty:
                        # Si hay coincidencia, obtener el valor de la columna "mobileOG fasta Header"
                        downstream_genes.append(match["mobileOG fasta Header"].values[0])

                downstream_ids = [seq.id for seq in downstream_seqs if any(seq.id in line for line in blast_result)]

                # Guardar el valor en df_paninmune usando los valores de la columna "mobileOG fasta Header"
                df_paninmune.at[index, "downstream_neighboring_genes"] = ",".join(downstream_genes)
                df_paninmune.at[index, "downstream_neighboring_pos"] = ",".join(downstream_ids)

    return df_paninmune

# Llamar a la función con los DataFrames adecuados
df_paninmune_updated = neighboring_genes_systems(df_paninmune, input_multifasta_file, diamond_db_name, df_mobilome)
#print(df_paninmune_updated)

##################################################################################################################################

df_paninmune_updated["init_island_pos"] = ""
df_paninmune_updated["end_island_pos"] = ""
df_paninmune_updated["init_all_pos_around_island"] = ""
df_paninmune_updated["end_all_pos_around_island"] = ""

def neighboring_pos_final(df):
    for index, row in df.iterrows():
        # Procesar upstream_neighboring_pos
        upstream_ids = row["upstream_neighboring_pos"].split(',')
        upstream_numbers = [int(re.search(r'_(\d+)', uid).group(1)) for uid in upstream_ids if re.search(r'_(\d+)', uid)]
        if upstream_numbers:
            df.at[index, "init_island_pos"] = min(upstream_numbers)
            df.at[index, "init_all_pos_around_island"] = ",".join(map(str, upstream_numbers))  # Almacenar todos los números

        # Procesar downstream_neighboring_pos
        downstream_ids = row["downstream_neighboring_pos"].split(',')
        downstream_numbers = [int(re.search(r'_(\d+)', uid).group(1)) for uid in downstream_ids if re.search(r'_(\d+)', uid)]
        if downstream_numbers:
            df.at[index, "end_island_pos"] = max(downstream_numbers)
            df.at[index, "end_all_pos_around_island"] = ",".join(map(str, downstream_numbers))  # Almacenar todos los números
    
    return df

df_paninmune_updated = neighboring_pos_final(df_paninmune_updated)

############################################################################################################################

# Creando las nuevas columnas en el DataFrame
df_paninmune_updated["upstream_genes_names"] = ""
df_paninmune_updated["downstream_genes_names"] = ""

def annotation_genes(df_paninmune, df_mobilome, result_blast_file="result_blast.txt"):
    # Leer el archivo result_blast.txt y tomar la segunda columna
    with open(result_blast_file, 'r') as f:
        blast_results = f.readlines()
    
    # Extraer la segunda columna de BLAST result
    blast_genes = [line.split()[1] for line in blast_results]

    # Crear diccionario para almacenar las anotaciones a partir de df_mobilome
    gene_annotations = {}
    
    # Buscar coincidencias en df_mobilome y llenar el diccionario
    for gene in blast_genes:
        match = df_mobilome[df_mobilome["mobileOG fasta Header"] == gene]
        if not match.empty:
            gene_annotations[gene] = match["Name"].values[0]
        else:
            gene_annotations[gene] = "No Match"

    # Procesar cada fila en df_paninmune
    for index, row in df_paninmune.iterrows():
        # Procesar upstream_neighboring_genes
        upstream_ids = row["upstream_neighboring_genes"].split(',')
        upstream_names = []
        for uid in upstream_ids:
            # Buscar anotaciones para cada upstream ID
            gene_name = gene_annotations.get(uid, "")
            upstream_names.append(gene_name)
        
        # Asignar los nombres encontrados a la columna de upstream
        df_paninmune.at[index, "upstream_genes_names"] = ",".join(upstream_names)

        # Procesar downstream_neighboring_genes
        downstream_ids = row["downstream_neighboring_genes"].split(',')
        downstream_names = []
        for did in downstream_ids:
            # Buscar anotaciones para cada downstream ID
            gene_name = gene_annotations.get(did, "")
            downstream_names.append(gene_name)

        # Asignar los nombres encontrados a la columna de downstream
        df_paninmune.at[index, "downstream_genes_names"] = ",".join(downstream_names)
    
    return df_paninmune

# Llamar la función con los DataFrames necesarios
df_paninmune_updated = annotation_genes(df_paninmune_updated, df_mobilome)

############################################################################################################################

columns_to_add = [
    "upstream_Major_mobileOG_Category",
    "upstream_Minor_mobileOG_Categories",
    "upstream_Taxonomy",
    "downstream_Major_mobileOG_Category",
    "downstream_Minor_mobileOG_Categories",
    "downstream_Taxonomy"
]

for col in columns_to_add:
    df_paninmune_updated[col] = ""

def mobilome_annotation(df_paninmune_updated, df_mobilome):
    for index, row in df_paninmune_updated.iterrows():
        
        # Proceso para upstream_genes_names
        upstream_ids = str(row["upstream_genes_names"]).split(",")
        upstream_major_categories = []
        upstream_minor_categories = []
        upstream_taxonomies = []

        for uid in upstream_ids:
            matching_row = df_mobilome[df_mobilome["Name"] == uid.strip()]
            if not matching_row.empty:
                major_category = matching_row["Major mobileOG Category"].iloc[0]
                minor_categories = str(matching_row["Minor mobileOG Categories"].iloc[0]).split(",") if not pd.isna(matching_row["Minor mobileOG Categories"].iloc[0]) else []
                taxonomies = str(matching_row["Taxonomy"].iloc[0]).split(",") if not pd.isna(matching_row["Taxonomy"].iloc[0]) else []
                
                if major_category:
                    upstream_major_categories.append(major_category)
                upstream_minor_categories.extend(minor_categories)
                upstream_taxonomies.extend(taxonomies)
        
        df_paninmune_updated.at[index, "upstream_Major_mobileOG_Category"] = ",".join(upstream_major_categories)
        df_paninmune_updated.at[index, "upstream_Minor_mobileOG_Categories"] = ",".join(upstream_minor_categories)
        df_paninmune_updated.at[index, "upstream_Taxonomy"] = ",".join(upstream_taxonomies)

        # Proceso para downstream_genes_names
        downstream_ids = str(row["downstream_genes_names"]).split(",")
        downstream_major_categories = []
        downstream_minor_categories = []
        downstream_taxonomies = []

        for did in downstream_ids:
            matching_row = df_mobilome[df_mobilome["Name"] == did.strip()]
            if not matching_row.empty:
                major_category = matching_row["Major mobileOG Category"].iloc[0]
                minor_categories = str(matching_row["Minor mobileOG Categories"].iloc[0]).split(",") if not pd.isna(matching_row["Minor mobileOG Categories"].iloc[0]) else []
                taxonomies = str(matching_row["Taxonomy"].iloc[0]).split(",") if not pd.isna(matching_row["Taxonomy"].iloc[0]) else []
                
                if major_category:
                    downstream_major_categories.append(major_category)
                downstream_minor_categories.extend(minor_categories)
                downstream_taxonomies.extend(taxonomies)
        
        df_paninmune_updated.at[index, "downstream_Major_mobileOG_Category"] = ",".join(downstream_major_categories)
        df_paninmune_updated.at[index, "downstream_Minor_mobileOG_Categories"] = ",".join(downstream_minor_categories)
        df_paninmune_updated.at[index, "downstream_Taxonomy"] = ",".join(downstream_taxonomies)

    return df_paninmune_updated


df_paninmune_updated = mobilome_annotation(df_paninmune_updated, df_mobilome)
df_paninmune_updated.to_csv("df_paninmune_final_rectificador.csv", index=False)
#print(df_paninmune_updated)

############################################################################################################################

