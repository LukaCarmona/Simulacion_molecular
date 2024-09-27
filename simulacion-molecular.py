'''Autor: Luka Carmona Rivas (i3B)'''

#comento el tema speck_mol al estar a la espera

# from st_speckmol import speck_plot
from pathlib import Path
from PIL import Image
from backend_molecular_250924 import calculate_outputs, write_hamiltonians
from graficador import create_graph
import streamlit as st
import matplotlib.pyplot as plt
import json
import time
import os
import base64
from io import BytesIO

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
            # print("distancia min", st.session_state.selected_range[0], "\n distancia max", st.session_state.selected_range[1])
            # guardo el step y range_values porque al ser un rango la función de resultado necesita parámetros distintos
            st.session_state.selected_range = range_values
            resultado = calculate_outputs(st.session_state.selected_molecule, archived_type, energy, st.session_state.selected_orbitas, st.session_state.selected_range[0], st.session_state.selected_range[1], st.session_state.selected_step)
            # guardo el resultado en una sesión para poder mantener los datos del gráfico 
            st.session_state.resultado = resultado
            
            # print("VALOREEES", st.session_state.selected_molecule, energy, st.session_state.selected_orbitas, resultado[0], resultado[2])
    
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
    image.save(buffer, format="PNG",help="Home")
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
    st.write("---------------------------------------------------------------------------------------------------------------------------------------------")
    molecula = st.selectbox("**Molécula**", moleculas, key='molecule')
    # carga de datos de molécula en base a molécula seleccionada en el select box
    datos_molecula = datos_json[molecula]['case_1']
    
    if archived_type == 0:
        if molecula == "LiH":
            st.markdown("""<span style='color: yellow;'>Las simulaciones pueden tardar al ser calculadas al momento</span>""", unsafe_allow_html=True)
        if molecula == "SnO":
            st.markdown("""<span style='color: yellow;'>Las simulaciones pueden tardar unos 2 minutos al ser calculadas al momento</span>""", unsafe_allow_html=True)
        if molecula == "H2S":
            st.markdown("""<span style='color: yellow;'>Las simulaciones pueden tardar de 3 a 4 minutos al ser calculadas al momento</span>""", unsafe_allow_html=True)
    
    
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
            step = st.number_input("**Seleccione el step para el gráfico**", min_value=0.1, max_value=1.0, value=0.3, step=0.1, format="%.1f")
        else:
            step = 0.3
            
        step = round(step,1)
        # st.write(step)
        num_values = round((distancias[1]-distancias[0])/ step)+1
        # st.write(num_values)

        # Array de distancias con los valores calculados
        new_distancias = []
        current_value = min_distancias
        
        for _ in range(num_values):
            new_distancias.append(round(current_value, 2))
            current_value += step
            
        # st.write(new_distancias)
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
        distancia_min = st.number_input('**Especifique la distancia en la que quiere calcular**', min_value=distancias[0], max_value=max(distancias), value=min(distancias), step=st.session_state.selected_step, format="%.1f")
        # distancia_min = round(distancia_min / 0.3) * 0.3
        distancia_min = round(distancia_min, 1)
        # st.write(distancia_min)
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
                # print([round(distancia_min,1)])
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
    col1, col2 = st.columns([1, 4])  # Ajusta el tamaño de las columnas según sea necesario

    # Colocar el botón "HOME" en la primera columna
    with col1:
        if st.button("HOME", key="home_button", help="Home", use_container_width=True):
            st.session_state.mostrar = False
            st.rerun()
    
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
                # st.write("Resultado 0", st.session_state.resultado[0])               
                # st.write("Resultado 1", st.session_state.resultado[1])
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
    titulo = '''
        <div style="background-color: #0E1117; padding: 20px; border-radius: 10px; box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);">
            <h1 style="color: #ad44ff; padding: 10px;">¡Bienvenidos a la Demo de Simulación Molecular!</h1>
            <div style="max-height: 400px; overflow-y: auto; padding: 10px; background-color: #0E1117; border: 1px solid #0E1117; border-radius: 10px;">
                <p style="color:#ad44ff">
                En esta demostración, exploraremos las simulaciones moleculares aprendiendo qué parámetros necesitamos para realizarlas 
                y los resultados que obtenemos en energías. <br><br> 
                Las simulaciones moleculares son representaciones de sistemas compuestos por átomos y moléculas. 
                Estas son calculadas por ordenadores utilizando la computación cuántica, que tiene el potencial de conseguir resultados 
                más precisos y de moléculas más complejas que con computación clásica. <br><br>
                Por esta razón, en Q4Real hemos hecho un código que compara los resultados conseguidos con el algoritmo Variational 
                Quantum Eigensolver (VQE) de computación cuántica con resultados exactos que sólo se pueden calcular clásicamente 
                para sistemas sencillos. <br><br> 
                Y ahora... No te quedes solo con la teoría, ¡entra y juega con las moléculas! <br><br><br>
                <strong>¿Cómo funciona?</strong><br>
                <strong>Escoge el tipo de ejecución:</strong><br>
                <ul style="margin-left: 10px">
                    <li><strong>Archivo:</strong> Para ver los resultados de la forma más rápida. Estos ya han estado calculados anteriormente.</li>
                    <li><strong>Simulación local:</strong> Para hacer los cálculos que tú quieras, en las distancias que tú especifiques. Puede tardar varios minutos en verse los resultados.</li>
                </ul>
                <br><br>
                <strong>Definición del sistema:</strong><br>
                <ul style="margin-left: 10px">
                    <li><strong>Molécula:</strong> Escoge la molécula que quieras simular.</li>
                    <li><strong>Espacio activo:</strong> Selecciona los electrones activos y los orbitales moleculares con los que quieras simular el sistema, 
                    o bien déjalos en la selección inicial para que los cálculos sean más rápidos.</li>
                    <li><strong>Distancias:</strong></li>
                    <ul style="margin-left: 10px">
                        <li><strong>Una sola distancia:</strong> Cálculo del valor de energía de la molécula en esta distancia. Se grafica el proceso de convergencia 
                        de la energía con el algoritmo VQE, la energía en función del número de iteraciones.</li>
                        <li><strong>Un rango de distancias:</strong> Cálculo de las energías por cada distancia dentro del rango introducido. Se grafica la energía 
                        en función de la distancia entre los átomos de la molécula.</li>
                    </ul>
                </ul>
                </p>
            </div>
        </div>
        '''

        # Mostrar el contenido formateado en Streamlit
    st.markdown(titulo, unsafe_allow_html=True)