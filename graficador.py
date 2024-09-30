# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 08:59:16 2024

@author: BICARILU
"""
import matplotlib.pyplot as plt
import streamlit as st

def create_graph(archived_type, option, selected_molecule, distancias, hartree_fall, energias, exacto, distancia_fin, distancia_inicio):
    
    plt.figure(facecolor='#0E1117')
    plt.rcParams.update({'font.size': 18, 'axes.facecolor': '#0E1117', 'figure.facecolor': '#0E1117', 'axes.edgecolor': 'white', 'axes.labelcolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white', 'legend.facecolor': '#0E1117', 'legend.edgecolor': 'white'})
    
    fig = plt.figure(figsize=(10, 7), facecolor='#0E1117')
    ax = fig.add_subplot(111, facecolor='#0E1117')
    ax.ticklabel_format(useOffset=False, style='plain') 
    
    if archived_type == 1:
        if option == "Un rango de distancias":
            ax.plot(distancias, hartree_fall, label = 'Hartree-Fock', linestyle = '--', linewidth=2, color = "white")
            ax.plot(distancias, energias, label='VQE ideal', marker='o', linewidth=0, color='#32C7AF', markersize=10)
            
            if selected_molecule != "Li2S":
                ax.plot(distancias, exacto, label = 'Exacto', linestyle='-', color = '#AD44FF')

            if st.session_state.selected_molecule == "LiH":
                molecula = "Li-H"
            elif st.session_state.selected_molecule == "SnO":
                molecula = "Sn-O"
            elif st.session_state.selected_molecule == "LiSH":
                molecula = "Li-S"
            elif st.session_state.selected_molecule == "H2S":
                molecula = "S-H"
            elif st.session_state.selected_molecule == "Li2S":
                molecula = "Li-S"
                
            ax.set_xlabel('$R_{'+molecula+'}$ [Å]')
            ax.set_ylabel('Energía [Ha]')
            ax.legend(labelcolor='#FFFFFF') 
            
            step = (distancia_fin - distancia_inicio) / 10 
            
            st.pyplot(fig)
        else: 
            ax.plot(distancias, energias, linestyle='-', color='darkviolet', label='VQE ideal')
        
            ax.set_xlabel('Iteraciones')
            ax.set_ylabel('Energía [Ha]')
            
            step = (distancia_fin - distancia_inicio) / 10 
            
            st.pyplot(fig)
            
    elif archived_type == 0:
        # st.write("EXACTO: ", exacto)
        # st.write("DISTANCIAS: ", distancias)
        if option == "Un rango de distancias":
            # if exacto != None: 
            #     ax.plot(distancias, hartree_fall, label = 'Hartree-Fock', linestyle = '--', linewidth=2, color = "white")
                
            ax.plot(distancias, energias, label='VQE ideal', marker='o', linewidth=0, color='#32C7AF', markersize=10)
            
            if selected_molecule != "Li2S" and exacto != None:
                ax.plot(distancias, exacto, label = 'Exacto', linestyle='-', color = '#AD44FF')
            
            if st.session_state.selected_molecule == "LiH":
                molecula = "Li-H"
            elif st.session_state.selected_molecule == "SnO":
                molecula = "Sn-O"
            elif st.session_state.selected_molecule == "H2S":
                molecula = "S-H"
                
            ax.set_xlabel('$R_{'+molecula+'}$ [Å]')
            ax.set_ylabel('Energía [Ha]')
            ax.legend(labelcolor='#FFFFFF') 
            
            step = (distancia_fin - distancia_inicio) / 10 
            
            st.pyplot(fig)
        else:
            ax.plot(distancias, energias, linestyle='-', color='darkviolet', label='VQE ideal')
        
            ax.set_xlabel('Iteraciones')
            ax.set_ylabel('Energía [Ha]')
            
            step = (distancia_fin - distancia_inicio) / 10 
            
            st.pyplot(fig)
