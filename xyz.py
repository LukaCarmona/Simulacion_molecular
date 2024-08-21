'''Autor: Luka Carmona Rivas (i3B)'''

#comento el tema speck_mol al estar a la espera

# from st_speckmol import speck_plot
from pathlib import Path
from PIL import Image
from backend_molecular_210824 import calculate_outputs, write_hamiltonians
from graficador import crete_graph
import streamlit as st
import matplotlib.pyplot as plt
import json
import time
import os
import base64

#-------------------------------------- SESIONES --------------------------------------------
# Inicializar estado de sesión para 'mostrar' si no está ya establecido
if 'mostrar' not in st.session_state:
    st.session_state.mostrar = False
    st.session_state.pulsado = False
    st.session_state.selected_molecule = ''
    st.session_state.resultado = ''
    st.session_state.selected_option = 'Un Punto'
    st.session_state.selected_range = (0, 0)
    
archived = 1
# -------------------------------------- ESTILOS ---------------------------------------------
htmlpath = Path(__file__).parent / 'style.css'

if htmlpath.exists():
    with open(htmlpath) as f:
        css = f.read()
    st.markdown(f'<style>{css}</style>', unsafe_allow_html=True)
else:
    st.error(f"El archivo CSS en {htmlpath} no se encontró.")

# -------------------------------------- CARGA DE DATOS --------------------------------------
# Cargar datos desde el JSON
with open("datos.json", "r") as data:
    datos_json = json.load(data)

# ------------------------------------- SIDE BAR ---------------------------------------------
with st.sidebar:
    #imagen de logo Q4Real
    imagen_path = Path(__file__).parent / 'logotipo-web-alpha.png'
    image = Image.open(imagen_path)
    st.image(image, use_column_width=True)
    
    # Cargar las opciones de moléculas desde el JSON
    moleculas = list(datos_json.keys())
    
    #select box de molecula
    molecula = st.selectbox("Molécula", moleculas, key='molecule')
    #carga de datos de molecula en base a molecula seleccionada en el select box
    datos_molecula = datos_json[molecula]['case_1']
    
    #carga de datos de select box en base a molecula seleccionada y contenido del json
    energias_fijas = datos_molecula['Electrones_activos']
    energy = st.selectbox("Electrones activos", energias_fijas, key='energy')
    
    #carga de datos de select box en base a molecula seleccionada y contenido del json
    numeros_orbitas = datos_molecula['Orbitales_moleculares']
    orbitas = st.selectbox("Orbitales moleculares", numeros_orbitas, key='orbits')

    #asignacion de columnas para el estilo del sidebar
    col1, col2 = st.columns(2)
    with col1:
        #selector de  tipo de distancia para generar los graficos 
        option = st.radio("Seleccione el tipo de distancia", ("Un Rango", "Un Punto"), key='option')
        distancias = datos_molecula['distance']
        
    if option == "Un Rango":
        #inicializo las variables para la hora de crear el slider
        new_distancias = []
        min_distancias = min(distancias)
        max_distancias = max(distancias)
        current_value = min_distancias
        #input para poder elegir el step del grafico
        if archived == 0:
            step = st.number_input("Seleccione el step para el gráfico", value=0.3, step=0.1)
        else:
            step = 0.3
            
        #calculo el numero de valores para poder rellenar el array de "new_distancias"
        num_values = int((max_distancias - min_distancias) / step) + 1
        
        #for que rellena el array de "new_distancias" con los puntos que tendra el slider 
        for _ in range(num_values):
            new_distancias.append(round(current_value, 2))
            current_value += step
       
        #creo el slider en vase a los valores anteriormente calculados        
        range_values = st.select_slider(
            "Selecciona un rango de distancias",
            options=(new_distancias),
            value=(min_distancias, max_distancias)
        )
        
    elif option == "Un Punto":
        #creo el input de tipo numerico para pasar solo una distancia que suma en funcion del step
        distancia_min = st.number_input('Especifique la distancia en la que quiere calcular', value=min(distancias), step=st.session_state.selected_step)
    
    col1, col2 = st.columns(2)

    with col1:
        #creo y compruebo el boton donde guardo las variables a los valores que quiero 
        if st.button('Aplicar cambios'):
            st.session_state.pulsado = True
            st.session_state.mostrar = True
            st.session_state.selected_electrones = energias_fijas
            st.session_state.selected_orbitas = orbitas
            st.session_state.selected_molecule = molecula
            st.session_state.selected_option = option
            
            if option == "Un Rango":
                #guardo el step y range_values por que al ser un rango la funcion de resultado necesita parametros distintos
                st.session_state.selected_step = step
                st.session_state.selected_range = range_values
                resultado = calculate_outputs(st.session_state.selected_molecule, archived, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0], st.session_state.selected_range[1], st.session_state.selected_step)
                #guardo el resultado en una sesion para poder mantener los datos del grafico 
                st.session_state.resultado = resultado
                #llamo a la funcion que crea los hamiltonianos
                write_hamiltonians(st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, resultado[0], resultado[2])
            else:
                st.session_state.selected_range = (distancia_min, distancia_min)
                resultado = calculate_outputs(st.session_state.selected_molecule, archived, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0])
                st.session_state.resultado = resultado
                write_hamiltonians(st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, [distancia_min], resultado[2])
                  
            st.rerun()
    #si se ha pulsado el boton se crea el boton de descargar hamiltonianos
    if st.session_state.pulsado:
        with col2:
            if st.session_state.selected_option == "Un Rango": 
                #dependiendo de la opcion elegida cambia el file_path
                file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],st.session_state.selected_step]}_nl1.txt"
                #habro el documento para guardar en una variable el contenido del archivo
                with open(file_path, 'r') as download_file:
                    file_content = download_file.read()
                if os.path.exists(file_path):
                    #creo elboton de descarga con la variable que tiene contenido del hamiltonianp que tiene que descargar 
                    btn = st.download_button(
                        label="Descargar Hamiltoniano",
                        data=file_content,
                        file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],st.session_state.selected_step]}_nl1.txt",
                        mime='text/plain'
                    )
                else:
                    st.write("No se ha podido crear el archivo")
            else:
                #dependiendo de la opcion elegida cambia el file_path
                file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[distancia_min]}_nl1.txt"
                #habro el documento para guardar en una variable el contenido del archivo
                with open(file_path, 'r') as download_file:
                    file_content = download_file.read()
                if os.path.exists(file_path):
                    #creo elboton de descarga con la variable que tiene contenido del hamiltonianp que tiene que descargar 
                    btn = st.download_button(
                        label="Descargar Hamiltoniano",
                        data=file_content,
                        file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[distancia_min]}_nl1.txt",
                        mime='text/plain'
                    )
                else:
                    #control de errores 
                    st.write("No se ha podido crear el archivo")
           
#--------------------------------------- CONTENIDO PRINCIPAL ---------------------------------------------          

if st.session_state.mostrar:
    if st.session_state.pulsado:
        #gif de espera
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
    #setear a false la sesion de pulsado y crear arrays de datos de objeto 3d, electrones y orbitas  
    st.session_state.pulsado = False
    edited_atoms = []
    electrons = []
    orbits = []
    
    #carga de archivo xyz en base a molecula escogida
    with open(f"Archivos_xyz/{st.session_state.selected_molecule}.xyz", 'r') as file:
        xyz_content = file.read()
# ------------------------------------- CONTENIDO 3D Y GRAFICO -----------------------------------------       
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
        # col1, col2 = st.columns(2)
       
        # with col1:
        #     # res = speck_plot(xyz_content) 
        #     pass       
# ------------------------------------ GRAFICO DE RANGO ------------------------------------------------
        # with col2:
            
        #creo un array con los nombres de las moleculas que saldran en el grafico 
        x_label_list = {"LiH": "Li-H", "SnO": "Sn-O", "H2S": "H-S", "LiSH": "Li-S", "Li2S": "Li-S"}
        #meto en una variable la molecula seleccionada
        selected_molecule = st.session_state.selected_molecule
        
        #funcion que recorre el array de nombres de moleculas y asigna la correspondiente
        def x_label_content():
            for key in x_label_list:
                if key.startswith(selected_molecule):
                    x_label = x_label_list[key] 
                    
                    return x_label
                
        if st.session_state.selected_option == "Un Rango" and len(st.session_state.resultado[1]) > 1:
            if st.session_state.selected_range is not None:
                
                energias_completas = st.session_state.selected_electrones
                distancia_inicio, distancia_fin = st.session_state.selected_range
                
                #guardo en variables los datos de las sesiones de respuesta
                energias = st.session_state.resultado[1][0]
                distancias = st.session_state.resultado[0]
                hartree_fall = st.session_state.resultado[1][1]
                exacto = st.session_state.resultado[1][2]                     
                
                #llamo a la funcion que da nombre al x_label
                x_label = x_label_content()

                #genero el grafico
                if len(distancias) == len(energias):
                    crete_graph(st.session_state.selected_molecule, distancias, hartree_fall, energias, exacto, distancia_fin, distancia_inicio)
                else:
                    st.error("Las listas de distancias y energías no tienen la misma longitud.")
                    
#-------------------------------------- GRAFICO DE PUNTO -----------------------------------------------------------                        
        else:
            energias_completas = st.session_state.selected_electrones
            distancia_inicio, distancia_fin = st.session_state.selected_range
            
            #guardo en variables los datos de las sesiones de respuesta
            energias = st.session_state.resultado[1][0]
            distancias = st.session_state.resultado[0]
                        
            #llamo a la funcion que da nombre al x_label
            x_label = x_label_content()
                
            #genero el grafico
            if len(distancias) == len(energias):
                crete_graph(st.session_state.selected_molecule, distancias, None, energias, None, distancia_fin, distancia_inicio)
            else:
                st.error("Las listas de distancias y energías no tienen la misma longitud.")

else:
    #mensaje de presentacion de la pagina
    var1 = __file__
    titulo = '<h1 style="color: #ad44ff; padding: 10px;">Esto es una página web para visualizar moléculas y ver su comportamiento además de poder ajustarlas y ver qué sucedería</h1>'
    st.markdown(titulo, unsafe_allow_html=True)
