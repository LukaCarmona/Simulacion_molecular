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


def calculate_outputs(name_mol: str, archived: int, active_electrons:int , molecular_orbitals:int, distance_min: float, distance_max: float = None, step: float = None, theta_min: float = None, theta_max: float = None, theta_step: float = None, nlayers:int = 1 ) -> tuple[list, list, list, str]:
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
    distance = distance_min
    if not (distance_max is None and step is None):
      raise ValueError("If it's computed only one distance, the variables 'distance_max' and 'step' should be 'None'.")
  else:
    distance = list(np.arange(distance_min, distance_max+step, step))

  geometry = [distance]

  # Theta definition
  if (theta_max is None or theta_step is None) and (theta_min != None):
    theta = theta_min
    geometry.append(theta)
    if not (theta_max is None and theta_step is None):
      raise ValueError("If it's computed only one theta value, the variables 'theta_max' and 'theta_step' should be 'None'. ")
  elif theta_min != None and theta_max != None and theta_step!= None:
    theta = list(np.arange(theta_min, theta_max+theta_step, theta_step))
    geometry.append(theta)


  #If we retrieve the archived results

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
      if np.any(np.isclose(distance_min, [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4])):
        distance, energy, hamiltonians = db.read_outputs(name_mol, distance, "content/SnO_energies-all_dist[0.7,3.5,0.3]_nl-1.txt", "content/SnO_hamiltonians_dist[0.7,3.8,0.3].txt", energies_indexs)
      else:
        raise ValueError("The minimum distance for archived SnO must be one of this list: [0.7, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.8, 3.1, 3.4]")

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

  elif archived == 0:


    if name_mol!= "LiH":
      # The LiH molecule is the only one that can be computed on the spot.
      raise ValueError("Only the LiH molecule can be calculated on the spot. Try to use the LiH molecule or set archived to 1.")

    else:
      raise ValueError("The program to compute the energies of LiH in the spot is being developped. We apologize for this inconvenience.")

  return distance, energy, hamiltonians

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
  if theta is None:
    file_name = f'{name_mol}_hamiltonians_ae{active_electrons}_mo{molecular_orbitals}_dist[{distance[0]:.1f}, {distance[-1]:.1f}, {(distance[1]-distance[0]):.1f}]_nl{nlayers}.txt'
  else:
    file_name = f'{name_mol}_hamiltonians_ae{active_electrons}_mo{molecular_orbitals}_theta[{distance[0]:.1f}, {distance[-1]:.1f}, {(distance[1]-distance[0]):.1f}]_nl{nlayers}.txt'

  f = open(file_name, "a")

  for index, hamiltonian in enumerate(hamiltonians):
    if theta is None:
      f.write(f'\n Hamiltonian {index} for distance={distance[index]}:\n {hamiltonian}\n')
    else:
      f.write(f'\n Hamiltonian {index} for theta={distance[index]}:\n {hamiltonian}\n')
    # I need to think if all the hamiltonians are going to be full written or not.

  f.close()

