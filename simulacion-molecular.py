'''Autor: Luka Carmona Rivas (i3B)'''

#comento el tema speck_mol al estar a la espera

# from st_speckmol import speck_plot
from pathlib import Path
from PIL import Image
from backend_molecular_171024 import calculate_outputs, write_hamiltonians
from graficador import create_graph
from streamlit_pdf_viewer import pdf_viewer
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
    st.session_state.selected_molecule = 'LiH'
    st.session_state.resultado = 0
    st.session_state.selected_option = 'Un rango de distancias'
    st.session_state.selected_range = (0, 0)
    st.session_state.archived_type = 'Archivo'
    
archived_type = 1
st.session_state.selected_step = 0.3
# -------------------------------------- ESTILOS ---------------------------------------------
htmlpath = Path(__file__).parent / 'style.css'

if htmlpath.exists():
    with open(htmlpath) as f:
        css = f.read()
    st.markdown(f'<style>{css}</style>', unsafe_allow_html=True)
else:
    st.error(f"El archivo CSS en {htmlpath} no se encontró.")

def aplicar_cambios():
    if st.button('▶', help="Aplicar cambios",use_container_width=True):
        st.session_state.pulsado = True
        st.session_state.mostrar = True
        st.session_state.calculated_molecule = molecula
        st.session_state.selected_electrones = energias_fijas
        st.session_state.selected_orbitas = orbitas
        st.session_state.selected_option = option
            
        if option == "Un rango de distancias":
            st.session_state.selected_range = (range_values[0], range_values[-1])
            # print("distancia min", st.session_state.selected_range[0], "\n distancia max", st.session_state.selected_range[1])
            # guardo el step y range_values porque al ser Un rango de distancias la función de resultado necesita parámetros distintos
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
def texto_correcto(selected_molecule):
    if selected_molecule == "LiH":
        return "Li-H"
    elif selected_molecule == "SnO":
        return "Sn-O"
    elif selected_molecule == "LiSH":
        return "Li-S"
    elif selected_molecule == "H2S":
        return "S-H"
    elif selected_molecule == "Li2S":
        return "Li-S"
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

    # archived_type = st.selectbox("**Ejecutar en**", ["Simulación local", "Archivo", "Ordenadores cuánticos online"], index=1, key='archived')
    titulo = '<h3 style="color: #FFFFFF; margin-bottom: -70px;">Ejecutar en</h3>'
    st.markdown(titulo, unsafe_allow_html=True)
    archived_type = st.selectbox("", ["Archivo", "Simulación local"], index=0, key='archived',help="Para más rapidez: Archivo. Para tener más control sobre los parámetros: Simulación Local. ")
    # archived_type = st.selectbox("**Ejecutar en**", ["Archivo", "Simulación local"], index=1, key='archived',help="tipo de ejecucion del programa")
    st.session_state.archived_type = archived_type
    if archived_type == "Simulación local":
        archived_type = 0
        st.session_state.archived_type = archived_type
        
    elif archived_type == "Archivo":
        archived_type = 1
        st.session_state.archived_type = archived_type

    st.write(st.session_state.archived_type)
    st.write(archived_type)

    # Cargar las opciones de moléculas desde el JSON     
    if st.session_state.archived_type == 0:
        keys = list(datos_json.keys())
        keys_to_remove = keys[-2:]
        for key in keys_to_remove:
            del datos_json[key]
            
    moleculas = list(datos_json.keys())
    # select box de molécula
    container1 = st.container(border=True)
    with container1:
        titulo = '<h3 style="color: #FFFFFF; margin-bottom: -70px;">Molécula</h3>'
        st.markdown(titulo, unsafe_allow_html=True)
        molecula = st.selectbox("", moleculas, key='molecule')
        st.session_state.selected_molecule = molecula
        
        # Carga de datos de molécula en base a la molécula seleccionada
        datos_molecula = datos_json[molecula]['case_1']
        
        # Advertencias sobre tiempos de simulación dependiendo de la molécula
        if st.session_state.archived_type == 0:
            if molecula == "LiH":
                st.markdown("<span style='color: yellow;'>Las simulaciones pueden tardar al ser calculadas al momento</span>", unsafe_allow_html=True)
            if molecula == "SnO":
                st.markdown("<span style='color: yellow;'>Las simulaciones pueden tardar unos 2 minutos al ser calculadas al momento</span>", unsafe_allow_html=True)
            if molecula == "H2S":
                st.markdown("<span style='color: yellow;'>Las simulaciones pueden tardar de 3 a 4 minutos al ser calculadas al momento</span>", unsafe_allow_html=True)
        
        # Electrones activos
        energias_fijas = datos_molecula['Electrones_activos']
        if st.session_state.archived_type == 0:
            titulo = '<h3 style="color: #FFFFFF; margin-bottom: -70px;">Electrones activos</h3>'
            st.markdown(titulo, unsafe_allow_html=True)
            energy = st.selectbox("", energias_fijas, key='energy_local', help="Selección del espacio activo.")
        else:
            st.write("**Electrones activos**" + ": " + str(energias_fijas[0]))
            energy = energias_fijas[0]

        # Orbitales moleculares
        numeros_orbitas = datos_molecula['Orbitales_moleculares']
        if st.session_state.archived_type == 0:
            titulo = '<h3 style="color: #FFFFFF; margin-bottom: -70px;">Orbitales moleculares</h3>'
            st.markdown(titulo, unsafe_allow_html=True)
            orbitas = st.selectbox("", numeros_orbitas, key='orbitas_local', help="Selección del espacio activo.")
        else:
            st.write("**Orbitales moleculares**" + ": " + str(numeros_orbitas[0]))
            orbitas = numeros_orbitas[0]
        
        # Selección de distancias
        st.write("---------------------------------------------------------------------------------------------------------------------------------------------")
        titulo = '<h3 style="color: #FFFFFF; margin-bottom: -70px; max-width: 1000px; width: 270px;">Selección de distancias</h3>'
        st.markdown(titulo, unsafe_allow_html=True, help="Escoge “Una sola distancia” para graficar el proceso de convergencia con el algoritmo de VQE y “Un rango de distancias” para la energía en función de las distancias seleccionadas entre los átomos de la molécula.")
        option = st.radio("", ("Una sola distancia", "Un rango de distancias"), key='option')
        
        distancias = datos_molecula['distance']
        
        # Si se selecciona un rango de distancias
        if option == "Un rango de distancias":
            new_distancias = []
            min_distancias = min(distancias)
            max_distancias = max(distancias)
            
            if st.session_state.archived_type == 0:
                step = st.number_input("Seleccione el intervalo entre distancias: ", min_value=0.1, max_value=1.0, value=0.3, step=0.1, format="%.1f")
            else:
                step = 0.3
            
            step = round(step, 1)
            num_values = round((distancias[1] - distancias[0]) / step) + 1
            
            # Crear el array de distancias
            current_value = min_distancias
            for _ in range(num_values):
                new_distancias.append(round(current_value, 2))
                current_value += step

            labels = {distancia: f"{distancia} Å" for distancia in new_distancias}
            molecule_text = texto_correcto(st.session_state.selected_molecule)

            range_values = st.select_slider(
                f'**Distancia {molecule_text} (Å):**',
                options=new_distancias,
                value=(new_distancias[0], new_distancias[-1]),
                format_func=lambda x: labels[x]  # Aplicar etiquetas
            )
            
            if range_values[0] == range_values[1]:
                option = "Una sola distancia"
                distancia_min = round(range_values[0], 1)
            
            st.session_state.selected_step = step

        # Si se selecciona una sola distancia
        elif option == "Una sola distancia":
            if st.session_state.archived_type == 1:
                molecule_text = texto_correcto(st.session_state.selected_molecule)
                array_distancias = [distancias[0]]
                for i in range(10):
                    array_distancias.append(round(array_distancias[i] + 0.3, 1))
                distancia_min = st.selectbox(f'**Distancia {molecule_text} (Å):**', array_distancias, index=0)
            else:
                molecule_text = texto_correcto(st.session_state.selected_molecule)
                distancia_min = st.number_input(f'**Distancia {molecule_text} (Å):**', min_value=distancias[0], max_value=max(distancias), value=min(distancias), step=st.session_state.selected_step, format="%.1f")
                distancia_min = round(distancia_min, 1)

        # Aplicar cambios si el botón ha sido pulsado
        if st.session_state.pulsado == True and st.session_state.mostrar == True:
            col1, col2 = st.columns(2)
            with col1:
                aplicar_cambios()
        else:
            aplicar_cambios()
        
        # Botón para descargar hamiltonianos
        if st.session_state.pulsado:
            with col2:
                if option == "Un rango de distancias":
                    file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{orbitas}_dist{[range_values[0], range_values[1], round(st.session_state.selected_step, 1)]}_nl1.txt"
                    with open(file_path, 'r') as download_file:
                        file_content = download_file.read()
                    if os.path.exists(file_path):
                        btn = st.download_button(
                            label="💾",
                            help="Descargar Hamiltoniano",
                            use_container_width=True,
                            data=file_content,
                            file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{orbitas}_dist{[range_values[0], range_values[1], round(st.session_state.selected_step, 1)]}_nl1.txt",
                            mime='text/plain'
                        )
                    else:
                        st.write("No se ha podido crear el archivo")
                else:
                    file_path = f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{orbitas}_dist{[round(distancia_min, 1)]}_nl1.txt"
                    with open(file_path, 'r') as download_file:
                        file_content = download_file.read()
                    if os.path.exists(file_path):
                        btn = st.download_button(
                            label="💾",
                            help="Descargar Hamiltoniano",
                            use_container_width=True,
                            data=file_content,
                            file_name=f"{st.session_state.selected_molecule}_hamiltonians_ae{energy}_mo{orbitas}_dist{[round(distancia_min, 1)]}_nl1.txt",
                            mime='text/plain'
                        )
                    else:
                        st.write("No se ha podido crear el archivo")

#--------------------------------------- CONTENIDO PRINCIPAL --------------------------------------------- 
if st.session_state.archived_type != archived_type:
    pass
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
        if st.session_state.selected_molecule == "Li2s":
            titulo = '<h1 style="color: #ad44ff; padding: 10px; margin-left: 150px;">Molécula Li<sub>2</sub>S</h1>'
        elif st.session_state.selected_molecule == "H2s":
            titulo = '<h1 style="color: #ad44ff; padding: 10px; margin-left: 150px;">Molécula H<sub>2</sub>S</h1>'
        else:
            titulo = f'<h1 style="color: #ad44ff; padding: 10px; margin-left: 150px;">Molécula {st.session_state.calculated_molecule}</h1>'

        # Mostrar el título
        st.markdown(titulo, unsafe_allow_html=True)
            
        # Título secundario según la opción seleccionada
        if st.session_state.selected_option == "Un rango de distancias":
            titulo2 = f'<p style="color: #FFFFFF; padding: 1px; margin-left: 10px; font-size: 24px; width: 800;">Energía de la molécula en función de la geometría</p>'
        else:
            titulo2 = f'<p style="color: #FFFFFF; padding: 1px; margin-left: 17px; font-size: 24px;">Convergencia del optimizador del VQE</p>'

        st.markdown(titulo2, unsafe_allow_html=True)
    
    
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
                
        energias_completas = st.session_state.selected_electrones
        distancia_inicio, distancia_fin = st.session_state.selected_range
        
        #guardo en variables los datos de las sesiones de respuesta
        energias = st.session_state.resultado[1][0]
        min_energia = round(min(energias),3)
        # st.write(energias)
        distancias = st.session_state.resultado[0]
        
        # st.write("Resultado 0", st.session_state.resultado[0])               
        # st.write("Resultado 1", st.session_state.resultado[1])
        #llamo a la funcion que da nombre al x_label
        x_label = x_label_content()

        #genero el grafico
        # st.write("AAAAAA",len(distancias))
        # st.write("BBBBBB",len(energias))
        
        if st.session_state.selected_option == "Un rango de distancias" and len(st.session_state.resultado[1]) > 1:
            if st.session_state.selected_range is not None:
                hartree_fall = st.session_state.resultado[1][1]
                exacto = st.session_state.resultado[1][2]      
                if len(distancias) == len(energias):
                    create_graph(st.session_state.archived_type, st.session_state.selected_option, st.session_state.calculated_molecule, distancias, hartree_fall, energias, exacto, distancia_fin, distancia_inicio)
                    st.markdown(titulo, unsafe_allow_html=True)
                else:
                    st.error("Las listas de distancias y energías no tienen la misma longitud.")
                
#-------------------------------------- GRAFICO DE PUNTO -----------------------------------------------------------                        
        else:
            if len(distancias) == len(energias):
                create_graph(st.session_state.archived_type, st.session_state.selected_option, st.session_state.calculated_molecule, distancias, None, energias, None, distancia_fin, distancia_inicio)
                titulo = f'<p style="color: #ffffff;margin-left: 93px; font-size: 22px;">⚪ Energía mínima: {min_energia} Ha</p>'
                st.markdown(titulo, unsafe_allow_html=True)
            else:
                st.error("Las listas de distancias y energías no tienen la misma longitud.")

else:
    #mensaje de presentacion de la pagina
    image_1 = "H2S_image.png"
    image_2 = "SnO_image.png"
    
    var1 = __file__
    titulo = '''
        <div style="background-color: #0E1117; padding: 20px; border-radius: 10px; box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);">
            <h1 style="color: #ad44ff; padding: 10px;">¡Bienvenid@s a la Demo de Simulación Molecular!</h1>
            <div style="padding: 10px; background-color: #0E1117; border: 1px solid #0E1117; border-radius: 10px;">
                <p style="color:#ffffff">
                En esta demostración, exploraremos las <span style="color: #32c7af;">simulaciones moleculares</span> aprendiendo qué parámetros necesitamos para realizarlas 
                y los resultados que obtenemos en energías. <br><br> 
                Las simulaciones moleculares son representaciones de sistemas compuestos por átomos y moléculas. 
                Estas son calculadas por ordenadores utilizando la <span style="color: #32c7af;">computación cuántica</span>, que tiene el potencial de conseguir resultados 
                más precisos y de moléculas más complejas que con computación clásica. <br><br>
                Por esta razón, en Q4Real hemos hecho un código que compara los resultados conseguidos con el algoritmo <span style="color: #32c7af;">Variational 
                Quantum Eigensolver (VQE)</span> de computación cuántica con resultados exactos que sólo se pueden calcular clásicamente 
                para sistemas sencillos. <br><br> 
                Y ahora... No te quedes solo con la teoría, ¡entra y juega con las moléculas! <br><br>
            </div>
        </div>
    '''
    st.markdown(titulo, unsafe_allow_html=True)
    col1, col2 = st.columns([0.5,0.5])
    with col1:
        st.image(image_1, width=350)
    with col2:
        st.image(image_2, width=450)
                
    titulo2 = '''
        <div style="background-color: #0E1117; padding: 20px; border-radius: 10px; box-shadow: 0 4px 8px rgba(0, 0, 0, 0.1);">
            <div style="padding: 10px; background-color: #0E1117; border: 1px solid #0E1117; border-radius: 10px;">
                <h3 style="color: #ad44ff">¿Cómo funciona?</h3><br>
                <strong>Escoge el tipo de ejecución:</strong>
                <ul style="margin-left: 10px">
                    <li><strong>Archivo:</strong> Para ver los resultados de la forma más rápida. Estos ya han estado calculados anteriormente.</li>
                    <li><strong>Simulación local:</strong> Para hacer los cálculos que tú quieras, en las distancias que tú especifiques. Puede tardar varios minutos en verse los resultados.</li>
                </ul>
                <br>
                <strong>Definición del sistema:</strong><br>
                <ul style="margin-left: 10px">
                    <li><strong>Molécula:</strong> Escoge la molécula que quieras simular.</li>
                    <li><strong>Espacio activo:</strong> Selecciona los electrones activos y los orbitales moleculares con los que quieras simular el sistema, 
                    o bien déjalos en la selección inicial para que los cálculos sean más rápidos.</li>
                    <li><strong>Distancias:</strong></li>
                    <ul style="margin-left: 10px">
                        <li><strong>Una sola distancia:</strong> Cálculo del valor de energía de la molécula en esta distancia. Se grafica el proceso de convergencia 
                        de la energía con el algoritmo VQE, la energía en función del número de iteraciones.</li>
                        <li><strong>Un rango de distancias de distancias:</strong> Cálculo de las energías por cada distancia dentro del rango introducido. Se grafica la energía 
                        en función de la distancia entre los átomos de la molécula.</li>
                    </ul>
                </ul>
                </p>
            </div>
        </div>
    '''

        # Mostrar el contenido formateado en Streamlit
    st.markdown(titulo2, unsafe_allow_html=True)
    # Q4Real-QMatter24-poster

    # Ruta del archivo PDF
    pdf_file_path = "Q4Real-QMatter24-poster.pdf"

    # Verificar si el archivo existe
    st.write("**Aquí os mostramos el póster presentado en la cuarta edición de la conferencia internacional de Quantum Matter:**")
    try:
        with open(pdf_file_path, "rb") as pdf_file:
            PDFbyte = pdf_file.read()
    except FileNotFoundError:
        st.error("El archivo PDF no se encontró. Asegúrate de que la ruta sea correcta.")
        st.stop()
    except Exception as e:
        st.error(f"Error al abrir el archivo PDF: {e}")
        st.stop()

    # Mostrar el PDF utilizando streamlit_pdf_viewer
    pdf_viewer(PDFbyte)

    # Botón para descargar el PDF
    st.download_button(
        label="Descargar PDF", 
        data=PDFbyte, 
        use_container_width=True,
        file_name="archivo.pdf", 
        mime="application/pdf"
    )