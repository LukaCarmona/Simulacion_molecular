# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:20:51 2024
@author: BICARILU
"""

import streamlit as st
# from st_speckmol import speck_plot
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image
import json
import time
from backend_molecular_180724 import calculate_outputs, write_hamiltonians
# from mock_backend_molecular_py import calculate_outputs, write_hamiltonians
import os
import base64
import subprocess
 
# Using system() method to
# execute shell commands
subprocess.Popen('echo "Geeks 4 Geeks"', shell=True)

# Inicializar estado de sesión para 'mostrar' si no está ya establecido
if 'mostrar' not in st.session_state:
    st.session_state.mostrar = False
    st.session_state.pulsado = False
    st.session_state.selected_molecule = ''
    st.session_state.selected_step = 0.3
    st.session_state.resultado = ''
    st.session_state.selected_option = 'Un Punto'
    st.session_state.selected_range = (0, 0)

# Path al archivo CSS
htmlpath = Path(__file__).parent / 'style.css'

# Cargar el archivo CSS
if htmlpath.exists():
    with open(htmlpath) as f:
        css = f.read()
    st.markdown(f'<style>{css}</style>', unsafe_allow_html=True)
else:
    st.error(f"El archivo CSS en {htmlpath} no se encontró.")

# Subir archivo XYZ (ejemplo)
# Cargar datos desde el JSON
with open("datos.json", "r") as data:
    datos_json = json.load(data)

# Mostrar la lista de átomos actual y permitir la edición
edited_atoms = []
electrons = []
orbits = []

with st.sidebar:
    imagen_path = Path(__file__).parent / 'logotipo-web-alpha.png'
    image = Image.open(imagen_path)
    st.image(image, use_column_width=True)
    
    # Cargar las opciones de moléculas desde el JSON
    moleculas = list(datos_json.keys())
    molecula = st.selectbox("Molécula", moleculas, key='molecule')
        
    datos_molecula = datos_json[molecula]['case_1']
    
    energias_fijas = datos_molecula['Electrones_activos']
    energy = st.selectbox("Electrones activos", energias_fijas, key='energy')
    
    numeros_orbitas = datos_molecula['Orbitales_moleculares']
    orbitas = st.selectbox("Orbitales moleculares", numeros_orbitas, key='orbits')

    col1, col2 = st.columns(2)

    with col1:
        option = st.radio("Seleccione el tipo de distancia", ("Un Rango", "Un Punto"), key='option')
        distancias = datos_molecula['distance']
    if option == "Un Rango":
        new_distancias = []
        min_distancias = min(distancias)
        max_distancias = max(distancias)
        current_value = min_distancias
        step = st.number_input("Seleccione el step para el gráfico", value=0.3, step=0.1)
            
        num_values = int((max_distancias - min_distancias) / step) + 1
        
        for _ in range(num_values):
            new_distancias.append(round(current_value, 2))
            current_value += step
       
        range_values = (min(distancias), max(distancias))
       
        range_values = st.select_slider(
            "Selecciona un rango de distancias",
            options=(new_distancias),
            value=(min_distancias, max_distancias)
        )
        
    elif option == "Un Punto":
        distancia_min = st.number_input('Especifique la distancia en la que quiere calcular', value=min(distancias), step=st.session_state.selected_step)
    
    col1, col2 = st.columns(2)

    with col1:
        if st.button('Aplicar cambios'):
            st.session_state.pulsado = True
            st.session_state.mostrar = True
            st.session_state.selected_electrones = energias_fijas
            st.session_state.selected_orbitas = orbitas
            st.session_state.selected_molecule = molecula
            st.session_state.selected_option = option
            if option == "Un Rango":
                st.session_state.selected_step = step
                st.session_state.selected_range = range_values
                resultado = calculate_outputs(st.session_state.selected_molecule, 1, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0], st.session_state.selected_range[1], st.session_state.selected_step)
                st.session_state.resultado = resultado
                write_hamiltonians(st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, resultado[0], resultado[2])
            else:
                st.session_state.selected_range = (distancia_min, distancia_min)
                resultado = calculate_outputs(st.session_state.selected_molecule, 1, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0])
                st.session_state.resultado = resultado
                write_hamiltonians(st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, [distancia_min], resultado[2])
                  
            st.rerun()
            
          
    if st.session_state.pulsado:
        with col2:
            if st.session_state.selected_option == "Un Rango": 
                file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],st.session_state.selected_step]}_nl1.txt"
                with open(file_path, 'r') as download_file:
                    file_content = download_file.read()
                if os.path.exists(file_path):
                    btn = st.download_button(
                        label="Descargar Hamiltoniano",
                        data=file_content,
                        file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],st.session_state.selected_step]}_nl1.txt",
                        mime='text/plain'
                    )
                else:
                    st.write("No se ha podido crear el archivo")
            else:
                file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[distancia_min]}_nl1.txt"
                with open(file_path, 'r') as download_file:
                    file_content = download_file.read()
                if os.path.exists(file_path):
                    btn = st.download_button(
                        label="Descargar Hamiltoniano",
                        data=file_content,
                        file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[distancia_min]}_nl1.txt",
                        mime='text/plain'
                    )
                else:
                    st.write("No se ha podido crear el archivo")
           

if st.session_state.mostrar:
    if st.session_state.pulsado:
        col1, col2, col3 = st.columns(3)
        with col1: 
            ''
        with col2:
            with open("loading_gif.gif", "rb") as file_:
                contents = file_.read()
            
            # Codifica el contenido en base64
            data_url = base64.b64encode(contents).decode("utf-8")
            
            # Crea un contenedor vacío
            placeholder = st.empty()
            
            # Mostrar el GIF en el contenedor con CSS para centrarlo y agregar margen superior
            placeholder.markdown(
                f"""
                <div style="display: flex; justify-content: center; align-items: center; height: 55vh;">
                    <img src="data:image/gif;base64,{data_url}" alt="loading gif" style="margin-top: 25px;">
                </div>
                """,
                unsafe_allow_html=True,
            )
            
            # Espera 3 segundos
            time.sleep(3)
            
            # Limpia el contenedor (elimina el GIF)
            placeholder.empty()
        with col3:
            ''
        
    st.session_state.pulsado = False
    edited_atoms = []
    electrons = []
    orbits = []

    with open(f"Archivos_xyz/{st.session_state.selected_molecule}.xyz", 'r') as file:
        xyz_content = file.read()
        
    col1, col2, col3 = st.columns([0.2,1.5,0.2])

    # Definimos el contenido de col1 y col3 (vacíos según tu código)
    with col1:
        pass
    
    with col2:
        # Definimos el contenido de col2 con el contenido dinámico de la variable titulo
        titulo = f'<h1 style="color: #ad44ff; padding: 10px;">Grafico y molecula de {st.session_state.selected_molecule}</h1>'
        st.markdown(titulo, unsafe_allow_html=True)
    
    with col3:
        pass
    
    titulo = '<p style="color: #ad44ff; padding: 5px;">''</p>'
    st.markdown(titulo, unsafe_allow_html=True)
    
    with st.container():
        col1, col2 = st.columns(2)
        
        with col1:
            # res = speck_plot(xyz_content) 
            pass       

        with col2:
            if st.session_state.selected_option == "Un Rango" and len(st.session_state.resultado[1]) > 1:
                if st.session_state.selected_range is not None:
                    
                    energias_completas = st.session_state.selected_electrones
                    distancia_inicio, distancia_fin = st.session_state.selected_range
                    
                    energias = st.session_state.resultado[1][0]
                    distancias = st.session_state.resultado[0]
                    hartree_fall = st.session_state.resultado[1][1]
                    exacto = st.session_state.resultado[1][2]
    
                    x_label_list = {"LiH": "Li-H", "SnO": "Sn-O", "H2S": "H-S", "LiSH": "Li-S", "Li2S": "Li-S"}
    
                    selected_molecule = st.session_state.selected_molecule  
                    
                    x_label = None  
                    
                    for key in x_label_list:
                        if key.startswith(selected_molecule):
                            x_label = x_label_list[key] 
                            break 
                            
                    if len(distancias) == len(energias):
                        plt.figure(facecolor='#0E1117')
                        plt.rcParams.update({'font.size': 18, 'axes.facecolor': '#0E1117', 'figure.facecolor': '#0E1117', 'axes.edgecolor': 'white', 'axes.labelcolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white', 'legend.facecolor': '#0E1117', 'legend.edgecolor': 'white'})
                        
                        fig = plt.figure(figsize=(10, 7), facecolor='#0E1117')
                        ax = fig.add_subplot(111, facecolor='#0E1117')
                        
                        ax.plot(distancias, hartree_fall, label = 'Hartree-Fock', linestyle = '--', linewidth=2, color = "white")
                        ax.plot(distancias, energias, label='VQE ideal', marker='o', linewidth=0, color='#32C7AF', markersize=10)
                        ax.plot(distancias, exacto, label = 'Exacto', linestyle='-', color = '#AD44FF')
    
                        ax.set_xlabel('$R_{'+st.session_state.selected_molecule+'}$ [Å]')
                        ax.set_ylabel('Energía [Ha]')
                        ax.legend(labelcolor='#FFFFFF')
                        
                        step = (distancia_fin - distancia_inicio) / 10  
                        # ax.set_xticks([distancia_inicio + i * step for i in range(int((distancia_fin - distancia_inicio) / step) + 1)])
                        
                        st.pyplot(fig)
                    else:
                        st.error("Las listas de distancias y energías no tienen la misma longitud.")
                        
            else:
                energias_completas = st.session_state.selected_electrones
                distancia_inicio, distancia_fin = st.session_state.selected_range
                
                energias = st.session_state.resultado[1][0]
                distancias = st.session_state.resultado[0]

                x_label_list = {"LiH": "Li-H", "SnO": "Sn-O", "H2S": "H-S", "LiSH": "Li-S", "Li2S": "Li-S"}

                selected_molecule = st.session_state.selected_molecule  
                
                x_label = None  
                
                for key in x_label_list:
                    if key.startswith(selected_molecule):
                        x_label = x_label_list[key] 
                        break 
                        
                if len(distancias) == len(energias):
                    plt.figure(facecolor='#0E1117')
                    plt.rcParams.update({'font.size': 18, 'axes.facecolor': '#0E1117', 'figure.facecolor': '#0E1117', 'axes.edgecolor': 'white', 'axes.labelcolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white', 'legend.facecolor': '#0E1117', 'legend.edgecolor': 'white'})
                    
                    fig = plt.figure(figsize=(10, 7), facecolor='#0E1117')
                    ax = fig.add_subplot(111, facecolor='#0E1117')
                    
                    ax.plot(distancias, energias, linestyle='-', color='darkviolet', label='VQEideal')

                    ax.set_xlabel('Iteraciones')
                    ax.set_ylabel('Energía [Ha]')
                    
                    step = (distancia_fin - distancia_inicio) / 10  # Definir el tamaño del paso para los xticks
                    # ax.set_xticks([distancia_inicio + i * step for i in range(int((distancia_fin - distancia_inicio) / step) + 1)])
                    
                    st.pyplot(fig)
                else:
                    st.error("Las listas de distancias y energías no tienen la misma longitud.")

else:
    titulo = '<h1 style="color: #ad44ff; padding: 10px;">Esto es una página web para visualizar moléculas y ver su comportamiento además de poder ajustarlas y ver qué sucedería</h1>'
    st.markdown(titulo, unsafe_allow_html=True)
