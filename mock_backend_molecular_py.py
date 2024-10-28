# -*- coding: utf-8 -*-

'''
Creation: 03-07-2024
Autor: Marina Ristol Roura
Title: mock-backend-molecular-notebook

Notes:
 - This notebook is made to create the python file (.py) to try the proper functioning of the Streamlit app. In the end this file should be exchanged for the final and definitive program.
'''


# ------------------------ IMPORTS  ------------------------
from qiskit.quantum_info import SparsePauliOp
import numpy as np


# ------------------------ MAIN FUNCTION ------------------------

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

  # ---- Here the test starts ----
  if archived != 1:
    raise ValueError("This is a test! All the computations are archived! So the archived value must be 1.")

  if name_mol != "LiH":
    raise ValueError("This test is now prepared just for the LiH molecule.")

  else:
    #Distance and Energy
    distance = list(np.arange(0.5,3.9,0.3))
    energy = [[-7.033795929869496, -7.617372795046837, -7.809919870041346, -7.86159760722682, -7.858223185560653, -7.832221960355422, -7.7978862383531045, -7.762828147409343, -7.7338992836327165, -7.719239180399629, -7.723268135321398, -7.73947743537584]]



    #Read the hamiltonians (copied from the function read_hamiltonians)
    hamiltonians = [SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
               coeffs=[-6.55974945e+00+0.j,  1.58681931e-01+0.j,  1.58681931e-01+0.j,
                  -9.19064653e-03+0.j, -9.19064653e-03+0.j,  9.19064653e-03+0.j,
                  9.19064653e-03+0.j, -2.44956180e-01+0.j, -2.44956180e-01+0.j,
                  -3.58457668e-01+0.j,  1.91965777e-01+0.j,  9.19065976e-03+0.j,
                  9.19065976e-03+0.j, -9.19065976e-03+0.j, -9.19065976e-03+0.j,
                  -6.96832467e-03+0.j, -6.96832467e-03+0.j,  1.00717384e-02+0.j,
                  1.00717384e-02+0.j,  1.00717384e-02+0.j,  1.00717384e-02+0.j,
                  -1.29620254e-01+0.j,  5.53442255e-03+0.j,  5.53442255e-03+0.j,
                  5.53442255e-03+0.j,  5.53442255e-03+0.j, -5.53442255e-03+0.j,
                  -5.53442255e-03+0.j, -5.53442255e-03+0.j, -5.53442255e-03+0.j,
                  9.11904537e-03+0.j, -9.11904537e-03+0.j,  9.11904537e-03+0.j,
                  9.11904537e-03+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.00564566e+00+0.j,  1.20146971e-01+0.j,  1.20146971e-01+0.j,
                  -7.00943394e-03+0.j, -7.00943394e-03+0.j,  7.00943394e-03+0.j,
                  7.00943394e-03+0.j, -2.70089637e-01+0.j, -2.70089637e-01+0.j,
                  -3.95809074e-01+0.j,  2.13292023e-01+0.j,  7.00910689e-03+0.j,
                  7.00910689e-03+0.j, -7.00910689e-03+0.j, -7.00910689e-03+0.j,
                  -4.62861893e-03+0.j, -4.62861893e-03+0.j,  7.28004857e-03+0.j,
                  7.28004857e-03+0.j,  7.28004857e-03+0.j,  7.28004857e-03+0.j,
                  -1.27722120e-01+0.j,  5.06351642e-03+0.j,  5.06351642e-03+0.j,
                  5.06351642e-03+0.j,  5.06351642e-03+0.j, -5.06351642e-03+0.j,
                  -5.06351642e-03+0.j, -5.06351642e-03+0.j, -5.06351642e-03+0.j,
                  1.07710282e-02+0.j, -1.07710282e-02+0.j,  1.07710282e-02+0.j,
                  1.07710282e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.14018248e+00+0.j,  1.02365098e-01+0.j,  1.02365098e-01+0.j,
                  -1.02755819e-02+0.j, -1.02755819e-02+0.j,  1.02755819e-02+0.j,
                  1.02755819e-02+0.j, -2.77655755e-01+0.j, -2.77655755e-01+0.j,
                  -4.12030218e-01+0.j,  2.15066231e-01+0.j,  1.02755980e-02+0.j,
                  1.02755980e-02+0.j, -1.02755980e-02+0.j, -1.02755980e-02+0.j,
                  -4.88464760e-03+0.j, -4.88464760e-03+0.j,  6.61729111e-03+0.j,
                  6.61729111e-03+0.j,  6.61729111e-03+0.j,  6.61729111e-03+0.j,
                  -1.21015398e-01+0.j,  4.84325478e-03+0.j,  4.84325478e-03+0.j,
                  4.84325478e-03+0.j,  4.84325478e-03+0.j, -4.84325478e-03+0.j,
                  -4.84325478e-03+0.j, -4.84325478e-03+0.j, -4.84325478e-03+0.j,
                  1.05027255e-02+0.j, -1.05027255e-02+0.j,  1.05027255e-02+0.j,
                  1.05027255e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.18037576e+00+0.j,  1.01017450e-01+0.j,  1.01017450e-01+0.j,
                  -1.25500188e-02+0.j, -1.25500188e-02+0.j,  1.25500188e-02+0.j,
                  1.25500188e-02+0.j, -2.77276702e-01+0.j, -2.77276702e-01+0.j,
                  -4.22857650e-01+0.j,  2.10196426e-01+0.j,  1.25500506e-02+0.j,
                  1.25500506e-02+0.j, -1.25500506e-02+0.j, -1.25500506e-02+0.j,
                  -5.67999219e-03+0.j, -5.67999219e-03+0.j,  6.14444501e-03+0.j,
                  6.14444501e-03+0.j,  6.14444501e-03+0.j,  6.14444501e-03+0.j,
                  -1.14983173e-01+0.j,  4.79584468e-03+0.j,  4.79584468e-03+0.j,
                  4.79584468e-03+0.j,  4.79584468e-03+0.j, -4.79584468e-03+0.j,
                  -4.79584468e-03+0.j, -4.79584468e-03+0.j, -4.79584468e-03+0.j,
                  1.03491184e-02+0.j, -1.03491184e-02+0.j,  1.03491184e-02+0.j,
                  1.03491184e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.18531675e+00+0.j,  1.05973312e-01+0.j,  1.05973312e-01+0.j,
                  -1.47836570e-02+0.j, -1.47836570e-02+0.j,  1.47836570e-02+0.j,
                  1.47836570e-02+0.j, -2.73321198e-01+0.j, -2.73321198e-01+0.j,
                  -4.30952083e-01+0.j,  2.04278842e-01+0.j,  1.47836473e-02+0.j,
                  1.47836473e-02+0.j, -1.47836473e-02+0.j, -1.47836473e-02+0.j,
                  -7.07675247e-03+0.j, -7.07675247e-03+0.j,  5.72360782e-03+0.j,
                  5.72360782e-03+0.j,  5.72360782e-03+0.j,  5.72360782e-03+0.j,
                  -1.10459199e-01+0.j,  4.84581010e-03+0.j,  4.84581010e-03+0.j,
                  4.84581010e-03+0.j,  4.84581010e-03+0.j, -4.84581010e-03+0.j,
                  -4.84581010e-03+0.j, -4.84581010e-03+0.j, -4.84581010e-03+0.j,
                  1.03173275e-02+0.j, -1.03173275e-02+0.j,  1.03173275e-02+0.j,
                  1.03173275e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.17711881e+00+0.j,  1.12765963e-01+0.j,  1.12765963e-01+0.j,
                  -1.72672156e-02+0.j, -1.72672156e-02+0.j,  1.72672156e-02+0.j,
                  1.72672156e-02+0.j, -2.66918816e-01+0.j, -2.66918816e-01+0.j,
                  -4.36048501e-01+0.j,  1.97985180e-01+0.j,  1.72672048e-02+0.j,
                  1.72672048e-02+0.j, -1.72672048e-02+0.j, -1.72672048e-02+0.j,
                  -9.31031359e-03+0.j, -9.31031359e-03+0.j,  5.40424674e-03+0.j,
                  5.40424674e-03+0.j,  5.40424674e-03+0.j,  5.40424674e-03+0.j,
                  -1.07417747e-01+0.j,  4.98454858e-03+0.j,  4.98454858e-03+0.j,
                  4.98454858e-03+0.j,  4.98454858e-03+0.j, -4.98454858e-03+0.j,
                  -4.98454858e-03+0.j, -4.98454858e-03+0.j, -4.98454858e-03+0.j,
                  1.03350768e-02+0.j, -1.03350768e-02+0.j,  1.03350768e-02+0.j,
                  1.03350768e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.16555492e+00+0.j,  1.19883684e-01+0.j,  1.19883684e-01+0.j,
                  -2.00129906e-02+0.j, -2.00129906e-02+0.j,  2.00129906e-02+0.j,
                  2.00129906e-02+0.j, -2.58484108e-01+0.j, -2.58484108e-01+0.j,
                  -4.38265651e-01+0.j,  1.91329172e-01+0.j,  2.00129834e-02+0.j,
                  2.00129834e-02+0.j, -2.00129834e-02+0.j, -2.00129834e-02+0.j,
                  -1.27896738e-02+0.j, -1.27896738e-02+0.j,  5.23310199e-03+0.j,
                  5.23310199e-03+0.j,  5.23310199e-03+0.j,  5.23310199e-03+0.j,
                  -1.06112232e-01+0.j,  5.21078916e-03+0.j,  5.21078916e-03+0.j,
                  5.21078916e-03+0.j,  5.21078916e-03+0.j, -5.21078916e-03+0.j,
                  -5.21078916e-03+0.j, -5.21078916e-03+0.j, -5.21078916e-03+0.j,
                  1.03474114e-02+0.j, -1.03474114e-02+0.j,  1.03474114e-02+0.j,
                  1.03474114e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.15579444e+00+0.j,  1.26935932e-01+0.j,  1.26935932e-01+0.j,
                  -2.29211760e-02+0.j, -2.29211760e-02+0.j,  2.29211760e-02+0.j,
                  2.29211760e-02+0.j, -2.47903982e-01+0.j, -2.47903982e-01+0.j,
                  -4.37737326e-01+0.j,  1.84241849e-01+0.j,  2.29211717e-02+0.j,
                  2.29211717e-02+0.j, -2.29211717e-02+0.j, -2.29211717e-02+0.j,
                  -1.82847504e-02+0.j, -1.82847504e-02+0.j,  5.22993330e-03+0.j,
                  5.22993330e-03+0.j,  5.22993330e-03+0.j,  5.22993330e-03+0.j,
                  -1.07178384e-01+0.j,  5.52442501e-03+0.j,  5.52442501e-03+0.j,
                  5.52442501e-03+0.j,  5.52442501e-03+0.j, -5.52442501e-03+0.j,
                  -5.52442501e-03+0.j, -5.52442501e-03+0.j, -5.52442501e-03+0.j,
                  1.03080181e-02+0.j, -1.03080181e-02+0.j,  1.03080181e-02+0.j,
                  1.03080181e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.15146831e+00+0.j,  1.33718185e-01+0.j,  1.33718185e-01+0.j,
                  -2.54668370e-02+0.j, -2.54668370e-02+0.j,  2.54668370e-02+0.j,
                  2.54668370e-02+0.j, -2.34656242e-01+0.j, -2.34656242e-01+0.j,
                  -4.34169338e-01+0.j,  1.76749651e-01+0.j,  2.54668382e-02+0.j,
                  2.54668382e-02+0.j, -2.54668382e-02+0.j, -2.54668382e-02+0.j,
                  -2.68433728e-02+0.j, -2.68433728e-02+0.j,  5.38170508e-03+0.j,
                  5.38170508e-03+0.j,  5.38170508e-03+0.j,  5.38170508e-03+0.j,
                  -1.11454676e-01+0.j,  5.91728684e-03+0.j,  5.91728684e-03+0.j,
                  5.91728684e-03+0.j,  5.91728684e-03+0.j, -5.91728684e-03+0.j,
                  -5.91728684e-03+0.j, -5.91728684e-03+0.j, -5.91728684e-03+0.j,
                  1.01861672e-02+0.j, -1.01861672e-02+0.j,  1.01861672e-02+0.j,
                  1.01861672e-02+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.15468883e+00+0.j,  1.39785780e-01+0.j,  1.39785780e-01+0.j,
                  -2.65136274e-02+0.j, -2.65136274e-02+0.j,  2.65136274e-02+0.j,
                  2.65136274e-02+0.j, -2.18729038e-01+0.j, -2.18729038e-01+0.j,
                  -4.27084875e-01+0.j,  1.69506986e-01+0.j,  2.65136274e-02+0.j,
                  2.65136274e-02+0.j, -2.65136274e-02+0.j, -2.65136274e-02+0.j,
                  -3.90530602e-02+0.j, -3.90530602e-02+0.j,  5.63662696e-03+0.j,
                  5.63662696e-03+0.j,  5.63662696e-03+0.j,  5.63662696e-03+0.j,
                  -1.19263131e-01+0.j,  6.35254434e-03+0.j,  6.35254434e-03+0.j,
                  6.35254434e-03+0.j,  6.35254434e-03+0.j, -6.35254434e-03+0.j,
                  -6.35254434e-03+0.j, -6.35254434e-03+0.j, -6.35254434e-03+0.j,
                  9.97660915e-03+0.j, -9.97660915e-03+0.j,  9.97660915e-03+0.j,
                  9.97660915e-03+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.16410906e+00+0.j,  1.44492500e-01+0.j,  1.44492500e-01+0.j,
                  -2.51218889e-02+0.j, -2.51218889e-02+0.j,  2.51218889e-02+0.j,
                  2.51218889e-02+0.j, -2.01869028e-01+0.j, -2.01869028e-01+0.j,
                  -4.16854622e-01+0.j,  1.63924940e-01+0.j,  2.51219211e-02+0.j,
                  2.51219211e-02+0.j, -2.51219211e-02+0.j, -2.51219211e-02+0.j,
                  -5.38946562e-02+0.j, -5.38946562e-02+0.j,  5.92468170e-03+0.j,
                  5.92468170e-03+0.j,  5.92468170e-03+0.j,  5.92468170e-03+0.j,
                  -1.29409802e-01+0.j,  6.76179778e-03+0.j,  6.76179778e-03+0.j,
                  6.76179778e-03+0.j,  6.76179778e-03+0.j, -6.76179778e-03+0.j,
                  -6.76179778e-03+0.j, -6.76179778e-03+0.j, -6.76179778e-03+0.j,
                  9.69259432e-03+0.j, -9.69259432e-03+0.j,  9.69259432e-03+0.j,
                  9.69259432e-03+0.j]),
                SparsePauliOp(['III', 'ZIZ', 'ZZI', 'IZX', 'ZZX', 'ZXZ', 'IXZ', 'IIZ', 'IZI', 'ZII', 'IZZ', 'ZIX', 'IIX', 'IXI', 'ZXI', 'ZYY', 'IYY', 'XIZ', 'XZZ', 'XII', 'XZI', 'ZZZ', 'YYI', 'XXI', 'YYZ', 'XXZ', 'XZX', 'XIX', 'YZY', 'YIY', 'YXY', 'XYY', 'XXX', 'YYX'],
                              coeffs=[-7.17520444e+00+0.j,  1.47374071e-01+0.j,  1.47374071e-01+0.j,
                  -2.17920419e-02+0.j, -2.17920419e-02+0.j,  2.17920419e-02+0.j,
                  2.17920419e-02+0.j, -1.86884965e-01+0.j, -1.86884965e-01+0.j,
                  -4.05317606e-01+0.j,  1.60898762e-01+0.j,  2.17921648e-02+0.j,
                  2.17921648e-02+0.j, -2.17921648e-02+0.j, -2.17921648e-02+0.j,
                  -6.86585745e-02+0.j, -6.86585745e-02+0.j,  6.19062082e-03+0.j,
                  6.19062082e-03+0.j,  6.19062082e-03+0.j,  6.19062082e-03+0.j,
                  -1.39260525e-01+0.j,  7.07976980e-03+0.j,  7.07976980e-03+0.j,
                  7.07976980e-03+0.j,  7.07976980e-03+0.j, -7.07976980e-03+0.j,
                  -7.07976980e-03+0.j, -7.07976980e-03+0.j, -7.07976980e-03+0.j,
                  9.36971732e-03+0.j, -9.36971732e-03+0.j,  9.36971732e-03+0.j,
                  9.36971732e-03+0.j])]

    xyz_data = '''27
Li    2.041500    2.041500    2.041500
Li    2.041500    0.000000    0.000000
Li    2.041500    0.000000    4.083000
Li    2.041500    4.083000    0.000000
Li    2.041500    4.083000    4.083000
Li    0.000000    2.041500    0.000000
Li    0.000000    2.041500    4.083000
Li    4.083000    2.041500    0.000000
Li    4.083000    2.041500    4.083000
Li    0.000000    0.000000    2.041500
Li    0.000000    4.083000    2.041500
Li    4.083000    0.000000    2.041500
Li    4.083000    4.083000    2.041500
H    0.000000    0.000000    0.000000
H    0.000000    0.000000    4.083000
H    0.000000    4.083000    0.000000
H    0.000000    4.083000    4.083000
H    4.083000    0.000000    0.000000
H    4.083000    0.000000    4.083000
H    4.083000    4.083000    0.000000
H    4.083000    4.083000    4.083000
H    0.000000    2.041500    2.041500
H    4.083000    2.041500    2.041500
H    2.041500    0.000000    2.041500
H    2.041500    4.083000    2.041500
H    2.041500    2.041500    0.000000
H    2.041500    2.041500    4.083000
'''


  return distance, energy, hamiltonians, xyz_data



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
