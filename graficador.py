# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 08:59:16 2024

@author: BICARILU
"""
import matplotlib.pyplot as plt
import streamlit as st

def crete_graph(selected_molecule, distancias, hartree_fall, energias, exacto, distancia_fin, distancia_inicio):
    
    plt.figure(facecolor='#0E1117')
    plt.rcParams.update({'font.size': 18, 'axes.facecolor': '#0E1117', 'figure.facecolor': '#0E1117', 'axes.edgecolor': 'white', 'axes.labelcolor': 'white', 'xtick.color': 'white', 'ytick.color': 'white', 'legend.facecolor': '#0E1117', 'legend.edgecolor': 'white'})
    
    fig = plt.figure(figsize=(10, 7), facecolor='#0E1117')
    ax = fig.add_subplot(111, facecolor='#0E1117')
    ax.ticklabel_format(useOffset=False, style='plain') 
 
    if exacto != None:
        ax.plot(distancias, hartree_fall, label = 'Hartree-Fock', linestyle = '--', linewidth=2, color = "white")
        ax.plot(distancias, energias, label='VQE ideal', marker='o', linewidth=0, color='#32C7AF', markersize=10)
        
        if selected_molecule != "Li2S":
            ax.plot(distancias, exacto, label = 'Exacto', linestyle='-', color = '#AD44FF')
        
        ax.set_xlabel('$R_{'+st.session_state.selected_molecule+'}$ [Å]')
        ax.set_ylabel('Energía [Ha]')
        ax.legend(labelcolor='#FFFFFF') 
        
        step = (distancia_fin - distancia_inicio) / 10 
        
        st.pyplot(fig)
    else: 
        ax.plot(distancias, energias, linestyle='-', color='darkviolet', label='VQEideal')
    
        ax.set_xlabel('Iteraciones')
        ax.set_ylabel('Energía [Ha]')
        
        step = (distancia_fin - distancia_inicio) / 10 
        
        st.pyplot(fig)