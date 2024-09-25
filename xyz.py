'''Autor: Luka Carmona Rivas (i3B)'''

#comento el tema speck_mol al estar a la espera

# from st_speckmol import speck_plot
from pathlib import Path
from PIL import Image
from backend_molecular_210824 import calculate_outputs, write_hamiltonians
from graficador import create_graph
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
    st.session_state.resultado = 0
    st.session_state.selected_option = 'Un Punto'
    st.session_state.selected_range = (0, 0)
    st.session_state.archived_type = 'Archivo'
    
archived_type = 1
# -------------------------------------- ESTILOS ---------------------------------------------
htmlpath = Path(__file__).parent / 'style.css'

if htmlpath.exists():
    with open(htmlpath) as f:
        css = f.read()
    st.markdown(f'<style>{css}</style>', unsafe_allow_html=True)
else:
    st.error(f"El archivo CSS en {htmlpath} no se encontró.")

def aplicar_cambios():
    if st.button('Aplicar cambios'):
        st.session_state.pulsado = True
        st.session_state.mostrar = True
        st.session_state.selected_electrones = energias_fijas
        st.session_state.selected_orbitas = orbitas
        st.session_state.selected_molecule = molecula
        st.session_state.selected_option = option

        if option == "Un Rango":
            st.session_state.selected_range = (range_values[0], range_values[-1])
            print("distancia min", st.session_state.selected_range[0], "\n distancia max", st.session_state.selected_range[1])
            # guardo el step y range_values porque al ser un rango la función de resultado necesita parámetros distintos
            st.session_state.selected_range = range_values
            resultado = calculate_outputs(st.session_state.selected_molecule, archived_type, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0], st.session_state.selected_range[1], st.session_state.selected_step)
            # guardo el resultado en una sesión para poder mantener los datos del gráfico 
            st.session_state.resultado = resultado
            
            print("VALOREEES", st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, resultado[0], resultado[2])
    
            # llamo a la función que crea los hamiltonianos
            write_hamiltonians(st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, resultado[0], resultado[2])
    
        else:
            st.session_state.selected_range = (distancia_min, distancia_min)
            resultado = calculate_outputs(st.session_state.selected_molecule, archived_type, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0])
            st.session_state.resultado = resultado
            write_hamiltonians(st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, [distancia_min], resultado[2])
    
        st.rerun()

# -------------------------------------- CARGA DE DATOS --------------------------------------
# Cargar datos desde el JSON
with open("datos.json", "r") as data:
    datos_json = json.load(data)

# ------------------------------------- SIDE BAR ---------------------------------------------
with st.sidebar:
    imagen_path = Path(__file__).parent / "logotipo-web-alpha.png"
    image = Image.open(imagen_path)

    # Encode the image to base64
    buffer = BytesIO()
    image.save(buffer, format="PNG")
    img_str = base64.b64encode(buffer.getvalue()).decode()

    # HTML for clickable image
    image_html = f"""
    <a href="" target="_self">
        <img src="data:image/png;base64,{img_str}" style="width: 100%;">
    </a>
    """
    st.markdown(image_html, unsafe_allow_html=True)

    archived_type = st.selectbox("**Ejecutar en**", ["Simulación local", "Archivo", "Ordenadores cuánticos online"], index=1, key='archived')
    st.session_state.archived_type = archived_type
    
    if archived_type == "Simulación local":
        archived_type = 0
        st.markdown("""<span style='color: yellow;'>Las simulaciones pueden tardar al ser calculadas al momento</span>""", unsafe_allow_html=True)
    elif archived_type == "Archivo":
        archived_type = 1
    else:
        archived_type = 2
    # Cargar las opciones de moléculas desde el JSON     
    if archived_type == 0:
        keys = list(datos_json.keys())
        keys_to_remove = keys[-2:]
        for key in keys_to_remove:
            del datos_json[key]
            
    moleculas = list(datos_json.keys())
    # select box de molécula
    molecula = st.selectbox("**Molécula**", moleculas, key='molecule')
    # carga de datos de molécula en base a molécula seleccionada en el select box
    datos_molecula = datos_json[molecula]['case_1']
    
    # carga de datos de select box en base a molécula seleccionada y contenido del json
    energias_fijas = datos_molecula['Electrones_activos']
    if archived_type == 0:
        energy = st.selectbox("**Electrones activos**", energias_fijas, key='energy_local')
    else:
        # energy = st.selectbox("**Electrones activos**", energias_fijas[0], key='energy')
        st.write("**Electrones activos**"+": "+ str(energias_fijas[0]))
        energy = energias_fijas[0]

    # carga de datos de select box en base a molécula seleccionada y contenido del json
    numeros_orbitas = datos_molecula['Orbitales_moleculares']
    if archived_type == 0:
        orbitas = st.selectbox("**Orbitales moleculares**", numeros_orbitas, key='orbitas_local')
    else:
        # orbitas = st.selectbox("**Orbitales moleculares**", numeros_orbitas[0], key='orbitas')
        st.write("**Orbitales moleculares**"+": "+ str(numeros_orbitas[0]))
        orbitas = numeros_orbitas[0]
    # asignación de columnas para el estilo del sidebar
    col1, col2 = st.columns(2)
    with col1:
        # selector de tipo de distancia para generar los gráficos 
        option = st.radio("**Seleccione el tipo de distancia**", ("Un Rango", "Un Punto"), key='option')
        distancias = datos_molecula['distance']
        
    if option == "Un Rango":
        # inicializo las variables para la hora de crear el slider
        new_distancias = []
        min_distancias = min(distancias)
        max_distancias = max(distancias)
        current_value = min_distancias
        # input para poder elegir el step del gráfico
        
        if archived_type == 0:
            step = st.number_input("**Seleccione el step para el gráfico**", min_value=0.1, max_value=1.0, value=0.3, step=0.1)
        else:
            step = 0.3
            
        num_values = round((distancias[1] / step)-distancias[0])
        # st.write(num_values)

        # Array de distancias con los valores calculados
        new_distancias = []
        current_value = min_distancias
        
        for _ in range(num_values):
            new_distancias.append(round(current_value, 2))
            current_value += step
      
        # Creación del slider en base a los valores calculados
        range_values = st.select_slider(
            "**Selecciona un rango de distancias**",
            options=new_distancias,
            value=(new_distancias[0], new_distancias[-1])
        )
        # step = st.number_input("Seleccione el step para el gráfico", min_value=0.3, max_value=3.0, value=0.3, step=0.1)     

        # Cálculo dinámico del step para que haya 12 distancias
        # step = (max_distancias - min_distancias) / (num_values_fixed - 1)
        st.session_state.selected_step = step

        # Cálculo del número de valores
        # num_values = num_values_fixed
      
       
        # print("rango valores", range_values)
    elif option == "Un Punto":
        # creo el input de tipo numérico para pasar solo una distancia que suma en función del step
        distancia_min = st.number_input('**Especifique la distancia en la que quiere calcular**', min_value=distancias[0], max_value=max(distancias), value=min(distancias), step=st.session_state.selected_step)
    
    col1, col2 = st.columns(2)
    
    with col1:
        # creo y compruebo el botón donde guardo las variables a los valores que quiero
        # if archived_type == 0:
        #     if option == "Un Punto":
        #         if distancia_min < 0.4 or distancia_min > max(distancias):
        #             st.error(f"Por favor seleccione una distancia válida entre 0.4 y {max(distancias)}.")
        #         else:
        #             aplicar_cambios()
        #     else:
        #         if 0.3 <= step <= 3.0 and step % 0.1 == 0:
        #             aplicar_cambios()
        
        #         else:
        #             st.error("Por favor seleccione un valor válido entre 0.3 y 3.0.")
        # else:
        #     if option == "Un Punto":
        #         if distancia_min < 0.4 or distancia_min > max(distancias):
        #             st.error(f"Por favor seleccione una distancia válida entre 0.4 y {max(distancias)}.")
        #         else:
        #             aplicar_cambios()
        #     else:
        #         aplicar_cambios()
        aplicar_cambios()

    #si se ha pulsado el boton se crea el boton de descargar hamiltonianos
    if st.session_state.pulsado:
        with col2:
            if st.session_state.selected_option == "Un Rango": 
                #dependiendo de la opcion elegida cambia el file_path
                if archived_type == 0:
                    file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],round(st.session_state.selected_step,1)]}_nl1.txt"
                else:
                    file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],round(st.session_state.selected_step,1)]}_nl1.txt"
                #habro el documento para guardar en una variable el contenido del archivo
                with open(file_path, 'r') as download_file:
                    file_content = download_file.read()
                if os.path.exists(file_path):
                    #creo elboton de descarga con la variable que tiene contenido del hamiltonianp que tiene que descargar 
                    btn = st.download_button(
                        label="Descargar Hamiltoniano",
                        data=file_content,
                        file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[st.session_state.selected_range[0],st.session_state.selected_range[1],round(st.session_state.selected_step,1)]}_nl1.txt",
                        mime='text/plain'
                    )
                else:
                    st.write("No se ha podido crear el archivo")
            else:
                #dependiendo de la opcion elegida cambia el file_path
                print([round(distancia_min,1)])
                file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[round(distancia_min,1)]}_nl1.txt"
                #habro el documento para guardar en una variable el contenido del archivo
                with open(file_path, 'r') as download_file:
                    file_content = download_file.read()
                if os.path.exists(file_path):
                    #creo elboton de descarga con la variable que tiene contenido del hamiltonianp que tiene que descargar 
                    btn = st.download_button(
                        label="Descargar Hamiltoniano",
                        data=file_content,
                        file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{st.session_state.selected_orbitas}_dist{[round(distancia_min,1)]}_nl1.txt",
                        mime='text/plain'
                    )
                else:
                    #control de errores 
                    st.write("No se ha podido crear el archivo")
           
#--------------------------------------- CONTENIDO PRINCIPAL ---------------------------------------------          
if st.session_state.mostrar:
    if st.button("🏚", key="home_button", help="Home", use_container_width=True):
        st.session_state.mostrar = False
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
                # st.write(energias)
                distancias = st.session_state.resultado[0]
                hartree_fall = st.session_state.resultado[1][1]
                exacto = st.session_state.resultado[1][2]                     
                
                #llamo a la funcion que da nombre al x_label
                x_label = x_label_content()

                #genero el grafico
                # st.write("AAAAAA",len(distancias))
                # st.write("BBBBBB",len(energias))
                
                if len(distancias) == len(energias):
                    create_graph(archived_type, st.session_state.selected_option, st.session_state.selected_molecule, distancias, hartree_fall, energias, exacto, distancia_fin, distancia_inicio)
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
                create_graph(archived_type, st.session_state.selected_option, st.session_state.selected_molecule, distancias, None, energias, None, distancia_fin, distancia_inicio)
            else:
                st.error("Las listas de distancias y energías no tienen la misma longitud.")

else:
    #mensaje de presentacion de la pagina
    var1 = __file__
    titulo = '<h1 style="color: #ad44ff; padding: 10px;">Esto es una página web para visualizar moléculas y ver su comportamiento además de poder ajustarlas y ver qué sucedería</h1>'
    st.markdown(titulo, unsafe_allow_html=True)
    