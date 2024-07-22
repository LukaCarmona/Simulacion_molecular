#Imports
import numpy as np
from pyscf import mcscf, scf, gto
from dataclasses import dataclass

import qiskit
from qiskit_nature.second_q.mappers import ParityMapper
from qiskit_nature.second_q.problems.base_problem import BaseProblem
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit.quantum_info import SparsePauliOp
from qiskit_algorithms.utils import algorithm_globals

# import quantumsymmetry as qs

import re


#---------------------------- FUNCTIONS FILE ---------------------------------------
# Is organized as follows:
# 1. Molecule functions: choose_molecule
# 2. Reading functions: read_hamiltonians, read_energies, read_params
# 3. Classes for save the results: TerminationeChecker, VQELog
# 4. Definition optimizer: def_optimizer
# 5. Definition ansatz: def_ansatz
# 6. Transpilation: transpilation
# 7. VQE
#------------ FUNCTIONS TO DEMO -----
# 8. Read energy and hamiltonians
# 9. All computation at the spot




#Molecule functions
def choose_molecule(
    name: str, bond_length: float, bond_length2: float = None, theta: float = None, phi1=None, phi2= None, mol_orb: int = None,
    ) -> tuple[ParityMapper, BaseProblem]:

    """Define a `name` compound separated by a given `bond_length` (in ang),
    with angle `theta` (in deg) and a fixed number of `mol_orb`."""

    if name=='LiH':  # definimos la estructura del LiH
        if mol_orb is None:
            mol_orb = 3  # por defecto: 3 orbitales moleculares

        driver = PySCFDriver(
            atom="Li 0 0 0; H 0 0 {}".format(bond_length),
            basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)

    elif name=='H2S':  # definimos la estructura del H2S
        if theta is None:
            theta = 92.1  # por defecto: 92.1 deg
        if mol_orb is None:
            mol_orb = 6  # por defecto: 6 orbitales moleculares

        theta = np.deg2rad(theta)
        y_cor, z_cor = bond_length*np.sin(theta/2), bond_length*np.cos(theta/2)
        driver = PySCFDriver(
            atom="S 0 0 0; H 0 {} {}; H 0 {} {}".format(y_cor, z_cor, -y_cor, z_cor),
            basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)

    elif name=='SnO':  # definimos la estructura del SnO
        if mol_orb is None:
            mol_orb = 5  # por defecto: 5 orbitales moleculares (chequear con los de Multiverse)

        driver = PySCFDriver(
            atom="Sn 0 0 0; O 0 0 {}".format(bond_length),
            basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)

    elif name=='H2':
        if mol_orb is None:
            mol_orb = 1

        driver = PySCFDriver(
            atom="H 0 0 0; H 0 0 {}".format(bond_length),
            basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)

    elif name=='Li2S':
        if theta is None:
            theta = 180
        if mol_orb is None:
            mol_orb = 12

        theta = np.deg2rad(theta)
        y_cor, z_cor = bond_length*np.sin(theta/2), bond_length*np.cos(theta/2)

        driver = PySCFDriver(
            atom="S 0 0 0; Li 0 {} {}; Li 0 {} {}".format(y_cor, z_cor, -y_cor, z_cor),
            basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)


    elif name=='LiSH':
        if theta is None:
            theta = 85 #Me estoy inventando estos números!! Cambiarlos
        if mol_orb is None:
            mol_orb = 9
        if bond_length2 is None: #Bond_length2 es la HS
            bond_length2 = 1.336 #Esta inventado cogido mas o menos de la Ref1
        if bond_length is None: #Bond_length es la LiS
            bond_length= 1.85 #Inventado igual

        theta = np.deg2rad(theta)
        y_cor, z_cor = bond_length*np.sin(theta/2), bond_length*np.cos(theta/2)
        y_cor2, z_cor2 = bond_length2*np.sin(theta/2), bond_length2*np.cos(theta/2)

        driver = PySCFDriver(
            atom="S 0 0 0; Li 0 {} {}; H 0 {} {}".format(y_cor, z_cor, -y_cor2, z_cor2), #0, bond_length, y_cor, z_cor),
            basis="sto3g", charge=0, spin=0, unit=DistanceUnit.ANGSTROM)

    elif name == 'LiCoO2': #Valores sacados del programa Vesta
        if bond_length is None: #Bond length between O-Co
            bond_length =1.921
        if bond_length2 is None: #Bond length between O-Li
            bond_length2 = 2.16448
        if theta is None: #Angle between O-Co-O
            theta = 86.2069
        if phi1 is None: #Angle between y and z coordinates
            phi1 = 47.412
        if phi2 is None: #Angle between x and z coordinates
            phi2 = 132.9443

        theta, phi1, phi2 = np.deg2rad(theta), np.deg2rad(phi1), np.deg2rad(phi2)

        x_Co, y_Co = bond_length*np.sin(theta/2), bond_length*np.cos(theta/2)
        x_O = 2*x_Co
        x_Li, y_Li, z_Li = bond_length2*np.cos(phi2), bond_length2*np.cos(phi1), bond_length2*np.sin(phi2)

        driver = PySCFDriver(
            atom = "O 0 0 0; Co {} {} 0; O {} 0 0 ; Li {} {} {}".format(x_Co, y_Co, x_O, x_Li, y_Li, z_Li),
            basis = "sto3g", charge = 0, spin = 0, unit = DistanceUnit.ANGSTROM
        )

    else:
        raise ValueError('Given \'name\' is not in the list. Try \'LiH\', ' +
                         '\'H2S\', \'SnO\' or \'Li2S\'.')

    molecule = driver.run()
    nuc_rep = molecule.nuclear_repulsion_energy
    active_space_trafo = ActiveSpaceTransformer(num_electrons= molecule.num_particles, num_spatial_orbitals= mol_orb)
    molecule = active_space_trafo.transform(molecule)
    return molecule, driver



# Read functions

def read_hamiltonians(name_file):
    f = open(name_file, "r")
    Lines = f.readlines()

    last_line = Lines[-1]

    count = 0
    hamiltonians = []
    # Strips the newline character
    coeff_beg= False
    for line in Lines:
        if line.strip().startswith('SparsePauliOp'):
            gates = line.strip()[14:]
            gates = gates[:-1]
            gates = eval(gates)

        if coeff_beg:
            coeff += line.strip()
            # print('coeff2', coeff)

        if line.strip().startswith('coeffs'):
            coeff_beg = True
            coeff = line.strip()[7:]
            # print('Coef1', coeff)

        # print("Line{}: {}".format(count, line.strip()))
        if line == '\n': # and count>1:
            if count==0:
                count+=1

            else:
                coeff_beg = False
                coeff = coeff[:-1]
                coeff = eval(coeff)

                # Construir el hamiltonià amb count-1
                hamiltonian = SparsePauliOp(gates, coeffs = coeff)
                # print(type(hamiltonian), hamiltonian)
                hamiltonians.append(hamiltonian)
                count +=1

        if line == last_line and coeff_beg:

            coeff_beg = False
            coeff = coeff[:-1]
            coeff = eval(coeff)

            # Construir el hamiltonià amb count-1
            hamiltonian = SparsePauliOp(gates, coeffs = coeff)
            # print(type(hamiltonian), hamiltonian)
            hamiltonians.append(hamiltonian)

    f.close()
    return hamiltonians


def read_energies(name_text):
    ' This is to read the energies in the text files that I created before, where there were just the energies form different types. '
    f = open(name_text, "r")
    Lines = f.readlines()

    count = 0
    x = []
    energies_vqe = []
    energies_rhf1 = []
    energies_exact = []

    for line in Lines:
        if line.strip().startswith('BL'):
            x.append(eval(line.strip()[4:]))
        elif line.strip().startswith('Theta'):
            x.append(eval(line.strip()[7:]))

        if line.strip().startswith('E VQE:'):
            energies_vqe.append(eval(line.strip()[7:]))

        if line.strip().startswith('E HF'):
            energies_rhf1.append(eval(line.strip()[6:]))

        if line.strip().startswith('Energy exact'): #it seems that it has to go without the space
            energies_exact.append(eval(line.strip()[14:]))
            
        elif line.strip().startswith('E exact'):
            energies_exact.append(eval(line.strip()[9:]))


    f.close()
    return x, energies_vqe, energies_rhf1, energies_exact


def read_params(name_text):
    f = open(name_text, "r")
    Lines = f.readlines()
    opt_params = []
    opt_val = 0
    opt_iter = 0

    iter = 0
    vals, params = [], []

    not_all = False

    for line in Lines:

        if line.strip().startswith('Vals:'):
            now_val = eval(line.strip()[5:])
            vals.append(now_val)
            # print('Val', now_val)
            iter += 1

            if now_val < opt_val:
                opt_val = now_val
                opt_iter = iter

        if not_all:
            prova =prova_not_all + " "+line.strip()
            prova = prova.split()
            prova.pop(0)
            lastitem = prova[-1]
            lastitem = lastitem[:-2]
            prova[-1] = lastitem
            now_params = []

            for i in prova:
              if i != "" and i.isspace()== False:
                now_params.append(eval(i))
            params.append(now_params)
            if iter == opt_iter:
              opt_params = now_params
            not_all = False

        if line.strip().startswith('Params:'):
            #Text to list
            if re.search(']',line.strip()[7:]):
              prova = line.strip()[7:].split()
              prova.pop(0)
              lastitem = prova[-1]
              lastitem = lastitem[:-2]
              prova[-1] = lastitem
              now_params = []
              for i in prova:
                if i != "" and i.isspace()== False:
                  now_params.append(eval(i))
              params.append(now_params)
            else:
              prova_not_all = line.strip()[7:]
              not_all = True
            if iter == opt_iter and not_all==False:
                opt_params = now_params

    f.close()
    return opt_val, opt_params, opt_iter, vals, params


class TerminationChecker:
    global vals_saved, f_termch

    def __init__(self, N : int, slope_min: float, e_min: float):
        self.N = N
        self.values = []
        self.slope_min = slope_min
        self.e_min = e_min
        self.iter_opt = 5000
    def __call__(self, nfev, parameters, value, stepsize, accepted) -> bool:

        print('- Iteration TermCh:', len(self.values))
        print('Last value:', value)
        print('Last parameters:', parameters)

        if len(self.values) > self.N:
            last_values = self.values[-self.N:]
            pp = np.polyfit(range(self.N), last_values, 1)
            slope = pp[0]
            f_termch.write(f'Slope:{slope}\n\n')
            if abs(slope)<self.slope_min:
                print("--- We are stopping: criteria achieved ---")
                return True

        if value<self.e_min: #Clip
            self.values.append(self.e_min)
            vals_saved.append(self.e_min)
            print('>>>>>>>   We clip the value to :', self.e_min)
            print('At iteration:', len(self.values))

            f_termch.write(f'\n - Iteration {len(self.values)}:')
            f_termch.write(f'\nVals:{self.e_min}\n')
            f_termch.write(f'\nParams:{parameters}\n\n')
            self.iter_opt = len(self.values)
        else:
            self.values.append(value)
            vals_saved.append(value)

            f_termch.write(f'\n - Iteration {len(self.values)}:')
            f_termch.write(f'\nVals:{value}\n')
            f_termch.write(f'\nParams:{parameters}\n\n')

        if len(self.values)> self.iter_opt+10:
            return True


        return False

@dataclass
class VQELog:
    global f_log
    values: list
    parameters: list
    def update(self, count, parameters, mean, _metadata):
        self.values.append(mean)
        self.parameters.append(parameters)
        f_log.write(f'\n - Iteration {count}:')
        f_log.write(f'\nVals:{mean}\n')
        f_log.write(f'\nParams:{parameters}\n\n')
        print(f"---------------Running circuit {count} of ~350", end="\r", flush=False)



# Optimizer
def def_optimizer(optimizer_params):
    from qiskit_algorithms.optimizers import SLSQP, SPSA
    algorithm_globals.random_seed = 10
    if optimizer_params[0]== "SPSA":
        if (len(optimizer_params)==7):
            optimizer = SPSA(maxiter = optimizer_params[1], termination_checker=TerminationChecker(optimizer_params[3][0],optimizer_params[3][1], optimizer_params[3][2]), learning_rate= optimizer_params[4], perturbation=optimizer_params[5], blocking = optimizer_params[6])
        elif (len(optimizer_params) ==4):
            optimizer = SPSA(maxiter = optimizer_params[1], termination_checker=TerminationChecker(optimizer_params[3][0], optimizer_params[3][1], optimizer_params[3][2]))
        elif (len(optimizer_params) ==3):
            optimizer = SPSA(maxiter = optimizer_params[1])

    elif optimizer_params[0] == "SLSQP":
        optimizer = SLSQP(maxiter = optimizer_params[1])

    return optimizer

#Computation ansatz

def def_ansatz(hamiltonian, driver, CAS_user, nlayers):
    from qiskit.circuit.library import TwoLocal
    #Preparing initial state of the ansatz (Hartree Fock)
    encoding = qs.make_encoding(atom=str(driver.atom), basis='sto-3g', CAS=(CAS_user[0],CAS_user[1]))

    #Qiskit Circuit
    hf_circ = qs.HartreeFockCircuit(encoding, atom=str(driver.atom), basis='sto-3g', CAS=(CAS_user[0],CAS_user[1]))
    print(' > Initial state: \n', hf_circ)

    ansatz = TwoLocal(
        num_qubits = hamiltonian.num_qubits,
        rotation_blocks = ["ry"],
        entanglement_blocks = "cx",
        entanglement = "linear",
        reps = nlayers,
        initial_state = hf_circ)

    return ansatz

#Transpilation
def transpilation(backend, options, ansatz, hamiltonian ):
    'Backend needs to be the backend, not only the name'
    from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
    target = backend.target
    pm = generate_preset_pass_manager(target=target, optimization_level= options.optimization_level)
    ansatz_transpiled = pm.run(ansatz)
    hamiltonian_transpiled = hamiltonian.apply_layout(layout=ansatz_transpiled.layout)
    return ansatz_transpiled, hamiltonian_transpiled

# #VQE

# def def_VQE(estimator, ansatz_transpiled, hamiltonian_transpiled, optimizer, initial_pt, log ):
#     'Function that calls the VQE function in qiskit'
#     from qiskit_algorithms import VQE
#     import time
#     time0 = time.time()
#     vqe = VQE(estimator, ansatz= ansatz_transpiled, optimizer= optimizer, initial_point=initial_pt, callback= log.update)
#     result = vqe.compute_minimum_eigenvalue(operator = hamiltonian_transpiled)
#     timef = time.time()-time0


#     return result.optimal_value, result.optimal_point #Energy, parameters


#Write results

def write_results(f, index, name_mol, distance, nlayers, optimizer_params, result, backend=None, theta = [False]):
    '''Function that writes the results for each BL or angle. Index = 0 if its the first entrance, 1 if its not.
    Theta has to be= [False/True, value]. If it s False it means that the loop for is done in the distance. Then the Theta can be assigned or not. If it s True, means that the loop
    is for the thetas. '''
    if index==0:
        f.write(f'> Molecule:{name_mol}\n')
        f.write(f' >> Number of layers:{nlayers}\n')
        f.write(f' >> Optimizer:{optimizer_params[0]}, Maxiter:{optimizer_params[1]}, Seed :{optimizer_params[2]}\n')
        if len(optimizer_params)>=4:
            f.write(f' >> TerminationChecker:{optimizer_params[3]}\n')
        elif len(optimizer_params) ==6:
            f.write(f' >> Learning rate and perturbation:{optimizer_params[4]},{optimizer_params[5]}.\n')
        if backend!=None:
            f.write(f' >> Backend:{backend.name}\n')
        if len(theta)>1 and theta[0]==False:
            f.write(f' >> Theta: {theta[1]:0.2f}\n')
        if theta[0]:
            f.write(f' >> Bond length: {distance:0.2f}\n')
        f.write(f'......\n\n')
    if theta[0]:
        f.write(f' >> Theta:{theta[1]:0.2f}\n') 
    else:
        f.write(f' >> Bond length:{distance:0.2f}\n')
    f.write(f' >> Results:\n')
    f.write(f'      >>> Energy:{result.optimal_value}\n')
    f.write(f'      >>> Last parameters: {result.optimal_point}\n\n')
    return f




def read_outputs(name_mol, distance, energy_file, hamiltonian_file, energies_indexs):
  '''
  This function makes all the work for reading the energy and hamiltonians files for each molecule. 
  '''

  if name_mol == "LiH":
    #With only one distance
    if len(distance)== 1:
      if distance[0] == 0.5 :
        result = read_params("convergencias/LiH_convergencias/LiH_qh_dist-0.5_nl-1_mit-600_seed-10_a-0_c-0_termch_hanoi2104.txt")
      elif distance[0] == 0.8 or distance ==1.1:
        result = read_params(f'convergencias/LiH_convergencias/LiH_qh_dist-{distance}_nl-1_mit-150_seed-10_a-0.9_c-0.3_termch_hanoi2404.txt')
      elif distance[0] == 1.4:
        result = read_params(f'convergencias/LiH_convergencias/LiH_qh_dist-1.4_nl-1_mit-250_seed-10_a-0.9_c-0.3_termch_hanoi2104.txt')
      elif distance[0] > 1.4:
        result = read_params(f'convergencias/LiH_convergencias/LiH_qh_dist-{distance}_nl-1_mit-120_seed-10_a-0.9_c-0.3_termch_hanoi2504_EM2_2.txt')
      
      energy = [result[3]]
      #We redifine the distance as the iterations, because they will be the "x" coordinate for the streamlit plot
      distance = range(len(energy[0]))

      index_min = int((distance[0]-0.5)/0.3)
      hamiltonians = [read_hamiltonians(hamiltonian_file)[index_min]]
      
    else:
      # I need to get the indexes of the distances
      index_min = int((distance[0]-0.5)/0.3)   #Here 0.3 is the step
      index_max = int((distance[-1]-0.5)/0.3)

      #VQE Energy
      energies = read_energies(energy_file) #[:][index_min:index_max+1]
      range_energies = [items[index_min:index_max+1] for items in energies]
      energy = []
      for i in energies_indexs:
         energy.append(range_energies[i])

      hamiltonians = read_hamiltonians(hamiltonian_file)[index_min:index_max+1]
      

  if name_mol == "SnO":
    if len(distance) ==1:
      result = read_params(f'convergencias/SnO_convergencias/SnO_dist-{distance[0]:0.2f}_log.txt')

      energy = [result[3]]
      #We redifine the distance as the iterations, because they will be the "x" coordinate for the streamlit plot
      distance = range(len(energy[0]))

      index_min = int((distance[0]-0.7)/0.3)
      hamiltonians = [read_hamiltonians(hamiltonian_file)[index_min]]


    else:
      index_min = int((distance[0]-0.7)/0.3)   #Here 0.3 is the step
      index_max = int((distance[-1]-0.7)/0.3)

      energies = read_energies(energy_file) #[:][index_min:index_max+1]
      range_energies = [items[index_min:index_max+1] for items in energies]
      energy = []
      for i in energies_indexs:
         energy.append(range_energies[i])

      hamiltonians = read_hamiltonians(hamiltonian_file)[index_min:index_max+1]
    
  if name_mol == "H2S":
    if len(distance) == 1:
      result = read_params(f'convergencias/H2S_convergencias/H2S_dist-{distance[0]:0.2f}_log.txt')
      
      energy = [result[3]]
      #We redifine the distance as the iterations, because they will be the "x" coordinate for the streamlit plot
      distance = range(len(energy[0]))

      index_min = int((distance[0]-0.5)/0.3)
      hamiltonians = [read_hamiltonians(hamiltonian_file)[index_min]]


    else:
      index_min = int((distance[0]-0.5)/0.3)   #Here 0.3 is the step
      index_max = int((distance[-1]-0.5)/0.3)
      
      energies = read_energies(energy_file) #[:][index_min:index_max+1]
      range_energies = [items[index_min:index_max+1] for items in energies]
      energy = []
      for i in energies_indexs:
         energy.append(range_energies[i])

      hamiltonians = read_hamiltonians(hamiltonian_file)[index_min:index_max+1]
    
  if name_mol == "LiSH":
    if len(distance) == 1:
      pass

    else:
      index_min = int((distance[0]-0.7)/0.3)   #Here 0.3 is the step
      index_max = int((distance[-1]-0.7)/0.3)
      
      energies = read_energies(energy_file) #[:][index_min:index_max+1]
      range_energies = [items[index_min:index_max+1] for items in energies]
      energy = []
      for i in energies_indexs:
         energy.append(range_energies[i])

      hamiltonians = read_hamiltonians(hamiltonian_file)[index_min:index_max+1]

  if name_mol == "Li2S":
    if len(distance) == 1:
      pass

    else:
      index_min = int((distance[0]-0.4)/0.3)   #Here 0.3 is the step
      index_max = int((distance[-1]-0.4)/0.3)
      
      energies = read_energies(energy_file) #[:][index_min:index_max+1]
      range_energies = [items[index_min:index_max+1] for items in energies]
      energy = []
      for i in energies_indexs:
         energy.append(range_energies[i])

      hamiltonians = read_hamiltonians(hamiltonian_file)[index_min:index_max+1]

  if name_mol == "LiCoO2":
    raise ValueError('This molecule is not prepared yet.')


    
  return distance, energy, hamiltonians




def compute_now(name_mol, nlayers, CAS_user, geometry, maxiter_vals, read_hams, computational_style: int = 0, options_runtime = None, backend=None ):
  '''
  This function computes the VQE energies of the name_mol at the spot.
  With computational_style = 0, the energies are computed with the Qiskit primitives, so the results are the "Ideal VQE" energies. 


  Parameters
  ----------
  Maxiter_vals should be a list of two values = [maxiter0, maxiter1].
      maxiter0 is the maximum number of iterations for the first distance.
      maxiter1 is the maximum number of iterations for the rest of the distances.
      
  '''
  import time
  from qiskit_algorithms import NumPyMinimumEigensolver
  import copy
  from qiskit_algorithms import VQE
  
  #Initializations
  time0 = time.time()
  solver = NumPyMinimumEigensolver()
  distance = geometry[0]


  #Reading hamiltonians
  if (read_hams[0]):
    hamiltonians = read_hamiltonians(read_hams[1])
  else:
    hamiltonians = []

  #Opening the file for all the distances
  f = open(f'{name_mol}_all-dist.txt',"a" )

  

  # -------   Computations with Qiskit Primitives (Ideal VQE)  --------
  if computational_style ==0 :

    #Optimizer params
    optimizer_name = "SLSQP"  #SPSA or SLSQP
    optimizer_a = False #0.9        #Learning rate of SPSA
    optimizer_c = False #0.3        #Perturbation of SPSA

    for index, dist in enumerate(distance):
      if index == 0:
        maxiter = [maxiter_vals[0]]
      else:
        maxiter.append(maxiter_vals[1])    #Maximum number of iterations. The first one is for the first distance and the second for the other ones

      if name_mol == "LiH":
        index_min = int((distance[0]-0.5)/0.3)   #Here 0.3 is the step
        index_max = int((distance[-1]-0.5)/0.3)
        
        lih_exact_energy = [-7.033795929869496, -7.617372795046837, -7.809919870041346, -7.86159760722682, -7.858223185560653, -7.832221960355422, -7.7978862383531045, -7.762828147409343, -7.7338992836327165, -7.719239180399629, -7.723268135321398, -7.73947743537584] #Minimum energy when the optimizer stops if it's equal or less (should'nt be)
        term_min_energy = lih_exact_energy[index_min:index_max+1]

      elif name_mol == "SnO":
        index_min = int((distance[0]-0.7)/0.3)   #Here 0.3 is the step
        index_max = int((distance[-1]-0.7)/0.3)  

        sno_exact_energy = [-6025.259601757683, -6034.538103527285, -6036.664177026194, -6037.132640923791, -6037.142936507439, -6037.02269921052, -6036.883047350734, -6036.85257501038, -6036.900510190118, -6036.822903278939]
        term_min_energy = sno_exact_energy[index_min:index_max+1]

      elif name_mol == "H2S":
        index_min = int((distance[0]-0.5)/0.3)   #Here 0.3 is the step
        index_max = int((distance[-1]-0.5)/0.3)  

        h2s_exact_energy =[-389.4364195665655, -393.33049227659166, -394.2327285133028, -394.3530834531598, -394.2777598781129, -394.17977480836146, -394.11228148305946, -394.08031705929903, -394.0688206429858, -394.065078128168, -394.0638973543976, -394.0630889537324, -394.0632640762127]
        term_min_energy = h2s_exact_energy[index_min:index_max+1] 
      
      term_num = 15
      term_error = 0.0000001
      termination_checker = [term_num, term_error, term_min_energy]

      #Optimization parameters
      seed = 10

      if optimizer_name == "SPSA":
        optimizer_params0 = [optimizer_name, maxiter, seed, termination_checker, optimizer_a, optimizer_c]
      else:
        optimizer_params0 = [optimizer_name, maxiter, seed, termination_checker]


  #####################
  optimizer_params = copy.deepcopy(optimizer_params0)

  energy = []

  for ind, dist in enumerate(distance):

   # Variables initialization 
   global vals_saved
   vals_saved = []
   log = VQELog([],[])

   #Definition optimizer
   if ind == 0:
      print('Optimizer params:', optimizer_params, '\n')
      if optimizer_params[0]== "SPSA":
        termination_checker = copy.deepcopy(optimizer_params[3][2])
      maxiter = np.copy(optimizer_params[1])
   if optimizer_params[0]== "SPSA":
      optimizer_params[3][2] = copy.deepcopy(termination_checker[ind])
   optimizer_params[1] = copy.deepcopy(maxiter[ind])

   if ind==0:
      optimizer = def_optimizer(optimizer_params[:4]) #aixooo ho he canviaaat, abans [:3]
   else:
      optimizer = def_optimizer(optimizer_params)

   #Definition molecule
   if len(geometry)>1:
     molecule, driver = choose_molecule(name_mol, bond_length = dist, theta = geometry[1])
   else:
     molecule, driver = choose_molecule(name_mol, bond_length = dist)

   #Hamiltonian
   if read_hams[0]:
      hamiltonian = hamiltonians[ind+index_min]
   else:
     hamiltonian =  qs.reduced_hamiltonian(atom = str(driver.atom), basis = 'sto-3g', output_format = 'qiskit', verbose = False, show_lowest_eigenvalue= False , CAS=(CAS_user[0], CAS_user[1]))
     hamiltonians.append(hamiltonian)

    #Ansatz and Options
   if ind==0:
     ansatz = def_ansatz(hamiltonian, driver, CAS_user, nlayers)
      
     #Set options
     if computational_style!=0:
        from qiskit_ibm_runtime import QiskitRuntimeService, Session, Options, Estimator
    
        options = Options()
        #Set number of shots, optimization_level and resilience_level
        options.execution.shots = options_runtime[0]
        options.optimization_level = options_runtime[1]
        options.resilience_level = options_runtime[2]
        options.seed_simulator = 42

    #Transpilation
   if backend!=None:
      ansatz_transpiled, hamiltonian_transpiled = transpilation(backend, options, ansatz, hamiltonian)

   #Creating files for results
   global f_log
   f_log = open(f'{name_mol}_dist-{dist:0.2f}_log.txt',"a" )
   # f_dist = open(f'{name_mol}_only-dist-{dist}.txt',"a")

   if optimizer_params[0]== "SPSA":
     global f_termch
     if len(geometry)>1:
        theta = geometry[1]
        f_termch = open(f'{name_mol}_dist-{distance:0.2f}_theta-{theta:0.2f}_termch.txt',"a")
     else:
        f_termch = open(f'{name_mol}_dist-{distance:0.2f}_termch.txt',"a")

  

   # Initial point
   if ind==0:
     np.random.seed(optimizer_params[2])
     initial_pt = [np.random.uniform( - np.pi, np.pi ) for i in range(ansatz.num_parameters)]
     print(' > Random initialization. Parameters:\n', initial_pt, '\n')
   else: 
     initial_pt =  result.optimal_point #NOT RESUUUUUULT
     print(' > Parameters for this distance:\n', initial_pt, '\n')

   # VQE
   if computational_style==0:
     from qiskit.primitives import Estimator as Estimator_primitives
     estimator = Estimator_primitives(options = {"simulator": {"seed_simulator": 42}})

   elif computational_style==3:
     estimator = Estimator(session = session,options=options)#solo tenía options, lo he cambiado con Estimator V2!!

   else:
     estimator = Estimator(session=session, options=options)


   if backend!=None:

      vqe = VQE(estimator, ansatz= ansatz_transpiled, optimizer= optimizer, initial_point=initial_pt, callback= log.update)
      result = vqe.compute_minimum_eigenvalue(operator = hamiltonian_transpiled)
   else:
     vqe = VQE(estimator, ansatz= ansatz, optimizer= optimizer, initial_point=initial_pt, callback= log.update)
     result = vqe.compute_minimum_eigenvalue(operator = hamiltonian)

   print(f' ---- Experiment with bond length={dist:0.2f} has ended ---- ')
   print(f' >  Final energy:', result.optimal_value)
   print(f' >  Final parameters:', result.optimal_point)
   energy.append(result.optimal_value)
   #Closing the Log and the TerminationChecker files
   f_log.close()
   if optimizer_params[0]== "SPSA":
     f_termch.close()

   #Writing into the file of only this distance and the all distance 
   # f_dist = dp.write_results(f, 0, name_mol, dist, nlayers, optimizer_params, result, backend)
   if len(geometry)>1:
     f = write_results(f, ind, name_mol, dist, nlayers, optimizer_params, result, backend, theta = [False, geometry[1]])
   else: 
     f = write_results(f, ind, name_mol, dist, nlayers, optimizer_params, result, backend, theta = [False])
  ttime = time.time()-time0
  print('\n \n ************************ FINISHED ************************')
  print(f'---- Total ellapsed time = {(time.time()-time0):0.2f} seconds. -----')
  f.close()

  return distance, energy, hamiltonians, ttime

