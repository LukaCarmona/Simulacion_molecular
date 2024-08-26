import qiskit
import pyscf
import qiskit_nature
# import qiskit_ibm_runtime
import qiskit_algorithms

import quantumsymmetry as qs


from pickle import FALSE
'''
Creation: 04-07-2024
Autor: Marina Ristol Roura
Title: backend-molecular-notebook

Notes:
 - This notebook is made to create the python file (.py) to create the outputs for the Streamlit app. This is the final and definitive program.
'''

import demo_background0 as db
import numpy as np
import time


def calculate_outputs(name_mol: str, archived: int, active_electrons:int , molecular_orbitals:int, distance_min: float, distance_max: float = None, step: float = None, theta_min: float = None, theta_max: float = None, theta_step: float = None, nlayers:int = 1 , init_params = None, backend = None, session = None) -> tuple[list, list, list, str]:
  '''
  Function that takes all the user inputs and gives the proper outputs to build the results in the Streamlit app.

  NOTES ON THE INPUTS
  -------------------
  - Archived:
      0 - Not archived. Compute the VQE as a quantum simulation (ideal).
      1 - It is archived. Don't need to compute anythin. Just read the results.
      2 - Not archived. Compute the VQE with quantum hardware (??????)
  - Distance:
      If the user only wants to compute one distance, the distance_min variable will be the one with the computed distance.
  - Theta:
      If it's none, is defined as the default value. If
  '''

  #Distance definition
  if distance_max is None or step is None:
    distance = [distance_min]
    if not (distance_max is None and step is None):
      raise ValueError("If it's computed only one distance, the variables 'distance_max' and 'step' should be 'None'.")
  else:
    distance = list(np.arange(distance_min, distance_max+(step)/2, step))

  geometry = [distance]

  # Theta definition
  if (theta_max is None or theta_step is None) and (theta_min != None):
    theta = [theta_min]
    geometry.append(theta)
    if not (theta_max is None and theta_step is None):
      raise ValueError("If it's computed only one theta value, the variables 'theta_max' and 'theta_step' should be 'None'. ")
  elif theta_min != None and theta_max != None and theta_step!= None:
    theta = list(np.arange(theta_min, theta_max+theta_step, theta_step))
    geometry.append(theta)


  # ---------- RETRIEVING THE ARCHIVED RESULTS ----------

  if archived == 1:
    energies_indexs = [1,2,3]
    # 1. LiH
    if name_mol == "LiH":
      if np.any(np.isclose(distance_min, [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8])):
        distance, energy, hamiltonians = db.read_outputs(name_mol, distance, "content/LiH_energies-all_dist[0.5,3.9,0.3].txt", "content/LiH_hamiltonians_dist[0.5,3.9,0.3].txt", energies_indexs)
      else:
        raise ValueError("The minimum distance for archived LiH must be one of this list: [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8]")

    # 2. SnO
    elif name_mol == "SnO":
      if np.any(np.isclose(distance_min, [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7])):
        distance, energy, hamiltonians = db.read_outputs(name_mol, distance, "content/SnO_energies-all_dist[0.7,3.8,0.3]_nl-1.txt", "content/SnO_hamiltonians_dist[0.7,3.8,0.3].txt", energies_indexs)
      else:
        raise ValueError("The minimum distance for archived SnO must be one of this list: [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7]")

    # 3. H2S
    elif name_mol == "H2S":
      if np.any(np.isclose(distance_min, [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1])):
        distance, energy, hamiltonians = db.read_outputs(name_mol, distance, "content/H2S_energies-all_dist[0.5,4.2,0.3]_theta-92.1-nl-1_mit-1000.txt", "content/H2S_hamiltonians_dist[0.5,3.9,0.3]_theta-92.1.txt", energies_indexs)
      else:
        raise ValueError("The minimum distance for archived H2S must be one of this list: [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1]")

    # 4. LiSH
    elif name_mol == "LiSH":
      if np.any(np.isclose(distance_min, [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7])):
        distance, energy, hamiltonians = db.read_outputs(name_mol, distance, "content/LiSH_energies-all_dist-s-li[0.7,3.8,.03]_dist-s-h-1.25312_theta-90_nl-2_mit-5000.txt", "content/LiSH_hamiltonians_dist-s-li[0.7,3.8,.03]_dist-s-h-1.25312_theta-90_nl-2.txt", energies_indexs)
      else:
        raise ValueError("The minimum distance for archived LiSH must be one of this list: [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7]")

    elif name_mol == "Li2S":
      if np.any(np.isclose(distance_min, [0.4, 0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4])):
        distance, energy, hamiltonians = db.read_outputs(name_mol, distance, "content/Li2S_energies-all_dist[0.4,3.5,0.3]_theta-180_nl-1_mit-1000.txt", "content/H2S_hamiltonians_dist[0.5,3.9,0.3]_theta-92.1.txt", energies_indexs)
      else:
        raise ValueError("The minimum distance for archived Li2S must be one of this list: [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4]")

    elif name_mol == "LiCoO2":
      raise ValueError('This molecule is not prepared yet.')


  #  ----------  COMPUTE THE ENERGIES AT THE SPOT -------------

  elif archived == 0 or archived == 2:
    # 1. LiH
    if name_mol == "LiH":

      read_hams =[True, "content/LiH_hamiltonians_dist[0.5,3.9,0.3].txt"]

      if np.any(np.isclose(distance_min, [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8])) and (step == 0.3 or step == None):
        r_energies = db.read_energies("content/LiH_energies-all_dist[0.5,3.9,0.3].txt")

        #Indexes
        index_min = round((distance[0]-0.5)/0.3)   #Here 0.3 is the step
        index_max = round((distance[-1]-0.5)/0.3)
      else:
        read_hams[0] = False
        index_min = None

    # 2. SnO
    elif name_mol == "SnO":
      read_hams =[True, "content/SnO_hamiltonians_dist[0.7,3.8,0.3].txt"]

      if np.any(np.isclose(distance_min, [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7])) and (step == 0.3 or step == None):
        r_energies = db.read_energies("content/SnO_energies-all_dist[0.7,3.8,0.3]_nl-1.txt")

        #Indexes
        index_min = round((distance[0]-0.7)/0.3)   #Here 0.3 is the step
        index_max = round((distance[-1]-0.7)/0.3)
      else:
        read_hams[0] = False
        index_min = None

    # 3. H2S
    elif name_mol == "H2S":
      read_hams =[True, 'content/H2S_hamiltonians_dist[0.5,3.9,0.3]_theta-92.1.txt']

      if np.any(np.isclose(distance_min, [0.5, 0.8, 1.1, 1.4, 1.7, 2.0, 2.3, 2.6, 2.9, 3.2, 3.5, 3.8, 4.1])) and (step == 0.3 or step == None):
        r_energies = db.read_energies("content/H2S_energies-all_dist[0.5,4.2,0.3]_theta-92.1-nl-1_mit-1000.txt")

        #Indexes
        index_min = round((distance[0]-0.5)/0.3)   #Here 0.3 is the step
        index_max = round((distance[-1]-0.5)/0.3)
      else:
        read_hams[0] = False
        index_min = None

    # 4. LiSH
    elif name_mol == "LiSH":
      read_hams =[True, 'content/LiSH_hamiltonians_dist-s-li[0.7,3.8,.03]_dist-s-h-1.25312_theta-90_nl-2.txt']

      if np.any(np.isclose(distance_min, [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4, 3.7])) and (step == 0.3 or step == None):
        r_energies = db.read_energies("content/LiSH_energies-all_dist-s-li[0.7,3.8,.03]_dist-s-h-1.25312_theta-90_nl-2_mit-5000.txt")

        #Indexes
        index_min = round((distance[0]-0.7)/0.3)   #Here 0.3 is the step
        index_max = round((distance[-1]-0.7)/0.3)
      else:
        read_hams[0] = False
        index_min = None

    elif name_mol == "Li2S":
      read_hams =[True, 'content/Li2S_hamiltonians_dist[0.4,3.5,0.3]_theta-180_nl-1.txt']

      if np.any(np.isclose(distance_min, [0.4, 0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4])) and (step == 0.3 or step == None):
        r_energies = db.read_energies("content/Li2S_energies-all_dist[0.4,3.5,0.3]_theta-180_nl-1_mit-1000.txt")

        #Indexes
        index_min = round((distance[0]-0.4)/0.3)   #Here 0.3 is the step
        index_max = round((distance[-1]-0.4)/0.3)
      else:
        read_hams[0] = False
        index_min = None


    else:
      raise ValueError("Only the LiH, SnO, H2S and LiSH molecules can be calculated on the spot. Try to use the LiH, SnO and H2S molecules or set archived to 1.")

    print('Read_hams', read_hams)
    # ---- COMPUTE THE ENERGIES WITH QISKIT PRIMITIVES ----
    if archived == 0:
      comp_style = 0
      backend = None
      session = None
      options_runtime = None
      max_iter = [500,200]
      
    elif archived == 2:
      comp_style = 3
      options_runtime = [500, 3, 2]
      max_iter = [2,2]#[500,200] #[5,5]

      #[options.execution.shots, options.optimization_level, options.resilience_level]
      #error_mitigation = 2      #Twirled Readout Error eXtinction (TREX) measurement twirling + Zero Noise Extrapolation (ZNE) and gate twirling
      

    if init_params == None:
      distance, energy_vqe, hamiltonians, ttime = db.compute_now(name_mol, nlayers, [active_electrons, molecular_orbitals], geometry,max_iter,
                                                    read_hams, computational_style = comp_style, options_runtime = options_runtime, backend = backend, session = session)
    else:
      distance, energy_vqe, hamiltonians, ttime = db.compute_now(name_mol, nlayers, [active_electrons, molecular_orbitals], geometry,max_iter,
                                                    read_hams, computational_style = comp_style, options_runtime = options_runtime, backend = backend, session = session,  init_params = init_params)

    if len(distance)>1:
      if index_min != None:
        energy = [energy_vqe, r_energies[2][index_min:index_max+1], r_energies[3][index_min:index_max+1]]
      else:
        energy = [energy_vqe]

    else:
        # energy = [energy_vqe, r_energies[2][index_min], r_energies[3][index_min]]

      result = db.read_params(f'{name_mol}_dist-{distance[0]:0.2f}_log.txt')

      energy = [result[3]]
      distance = range(len(energy[0]))

  


  return distance, energy, hamiltonians


# COMMAND ----------

def write_hamiltonians(name_mol: str, active_electrons:int , molecular_orbitals:int, distance: list, hamiltonians: list, theta = None, nlayers: int = 1):
  '''
  Function that creates a file and writes the hamiltonians in it. This function is agreed with the "read_hamiltonians" in the demo_functions.py file.

  NOTES ON THE INPUTS
  -------------------
  - Distance:
      They are the distances or the thetas for which the hamiltonians are computed. They are distances if the theta variable is None. Otherwise, if
      theta is True these values are the thetas.
  - Hamiltonians:
      List of SparsePauliOp for all the hamiltons for each distance.
  - File name:
      It has to be a .txt file! It should contain the name of the molecule and other properties of the computations.

  OUTPUT
  ------
  File .txt with all the information of the hamiltonians computed for each distance/theta.
  '''

  #We create the file
  if len(distance)>1:
    if theta is None:
      file_name = f'{name_mol}_hamiltonians_ae{active_electrons}_mo{molecular_orbitals}_dist[{distance[0]:.1f}, {distance[-1]:.1f}, {(distance[1]-distance[0]):.1f}]_nl{nlayers}.txt'
    else:
      file_name = f'{name_mol}_hamiltonians_ae{active_electrons}_mo{molecular_orbitals}_theta[{distance[0]:.1f}, {distance[-1]:.1f}, {(distance[1]-distance[0]):.1f}]_nl{nlayers}.txt'

  elif len(distance)==1:
    if theta is None:
      file_name = f'{name_mol}_hamiltonians_ae{active_electrons}_mo{molecular_orbitals}_dist[{distance[0]:.1f}]_nl{nlayers}.txt'
    else:
      file_name = f'{name_mol}_hamiltonians_ae{active_electrons}_mo{molecular_orbitals}_theta[{distance[0]:.1f}]_nl{nlayers}.txt'


  f = open(file_name, "w")

  if len(distance)>1:
    for index, hamiltonian in enumerate(hamiltonians):
        print("-------------------------------------------------------")
        print("print de la longitud del hamiltoniano",len(hamiltonians))
        print("-------------------------------------------------------")
        print('Index', index)
        print(distance)
        if theta is None:
            f.write(f'\n Hamiltonian {index} for distance={distance[index]}:\n {hamiltonian}\n')
        else:
            f.write(f'\n Hamiltonian {index} for theta={distance[index]}:\n {hamiltonian}\n')
        # I need to think if all the hamiltonians are going to be full written or not.
  else:
    f.write(f'\n Hamiltonian for distance={distance[0]}:\n {hamiltonians}\n')

  f.close()
