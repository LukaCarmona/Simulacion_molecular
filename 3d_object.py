# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:20:51 2024
@autor: BICARILU
"""

from ase import Atoms
import streamlit as st
from ase.io import write
from st_speckmol import speck_plot
import io
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from PIL import Image
import json


# Inicializar estado de sesión para 'mostrar' si no está ya establecido
if 'mostrar' not in st.session_state:
    st.session_state.mostrar = False

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
    # Cargar las opciones de moléculas desde el JSON
    moleculas = list(datos_json.keys())
    molecula = st.selectbox("Molécula", moleculas, key='molecule')

    with open(f"Archivos_xyz/{molecula}.xyz", 'r') as file:
        xyz_content = file.readlines()

    # Extraer el número de átomos de la primera línea
    try:
        num_atoms = int(xyz_content[0].strip())
    except ValueError:
        st.error("Error: La primera línea debe contener el número de átomos.")

    # Procesar los átomos
    atoms = []
    for line in xyz_content[1:num_atoms + 1]:
        parts = line.split()
        symbol = parts[0]
        x, y, z = map(float, parts[1:4])
        atoms.append((symbol, (x, y, z)))
        
    if molecula in datos_json:
        datos_molecula = datos_json[molecula]['case_1']
        energias_fijas = datos_molecula['Electrones_activos']
        energy = st.selectbox("Energia", energias_fijas, key='energy')
        
        numeros_orbitas = list(range(1, len(energias_fijas) + 1))  # Número de energías disponibles
        orbitas = st.selectbox("Órbitas", numeros_orbitas, key='orbits')

        # Distancias disponibles para esa molécula
        distancias = datos_molecula['distance']
        new_distancias = []
        current_value = 0.5
        while current_value <= 5.0:
            # Verificar si current_value está en distancias originales
            if any(abs(d - current_value) < 0.01 for d in distancias):
                new_distancias.append(round(current_value,2))
            current_value += 0.3
            
        range_values = (min(distancias), max(distancias))  # Rango automático basado en los datos

        col1, col2 = st.columns(2)

        with col1:
            option = st.radio("Seleccione el tipo de distancia", ("Un Rango", "Un Punto"))

        if option == "Un Rango":
            range_values = st.select_slider(
                "Selecciona un rango de distancias",
                options=(new_distancias),
                value=(min(distancias), max(distancias))
            )
        elif option == "Un Punto":
            distancia = st.number_input('Especifique la distancia en la que quiere calcular', value=min(distancias),
                                        step=0.1)

        # Botón para aplicar cambios
        if st.button('Aplicar cambios'):
            st.session_state.mostrar = True
            st.session_state.selected_electrones = energias_fijas  # Actualizar energías seleccionadas en la sesión
            st.session_state.selected_orbitas = orbitas
            if option == "Un Rango":
                st.session_state.selected_range = range_values
            else:
                st.session_state.selected_range = (distancia, distancia)

# Solo actualizar los átomos editados y otras variables si el botón es presionado
if st.session_state.mostrar:
    edited_atoms = []
    electrons = []
    orbits = []
    for i, (symbol, (x, y, z)) in enumerate(atoms):
        new_x = x
        new_y = y
        new_z = z
        new_electrons = 0  # Define aquí cómo se calcularían los nuevos valores de electrones y órbitas si es necesario
        new_orbits = 0
        edited_atoms.append((symbol, (new_x, new_y, new_z)))
        electrons.append(new_electrons)
        orbits.append(new_orbits)

    # Crear objeto Atoms de ASE con las propiedades editadas
    edited_atoms_ase = Atoms(
        symbols=[atom[0] for atom in edited_atoms],
        positions=[atom[1] for atom in edited_atoms]
    )

    # Guardar el archivo XYZ editado en un nuevo archivo temporal
    edited_xyz_file = io.StringIO()
    write(edited_xyz_file, edited_atoms_ase, format='xyz')

    titulo = '<p style="font-family:monospace; color:White; font-size: 32px;">Visualización 3D</p>'
    st.markdown(titulo, unsafe_allow_html=True)

    # Crear la visualización 3D
    with st.container():
        col1, col2 = st.columns(2)

        with col1:
            edited_xyz_content = edited_xyz_file.getvalue()
    
            # Opción para ajustar el tamaño del objeto 3D
            scaled_content = edited_xyz_content # Dividir por 2 para reducir a la mitad
    
            # Llamada a speck_plot con el contenido escalado
            speck_plot(scaled_content, wbox_height="350px",wbox_width="350px")
           
        with col2:
            ''