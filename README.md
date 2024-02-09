This repository contains the scripts and functions that simulate the dynamics from the papers:

This code can be redistributed and/or modified
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
 
This code is distributed by the authors in the hope that it will be 
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

If you use this code please cite the following paper:


(1) ”Dirac synchronization is rhythmic and explosive”
L. Calmon, J. G. Restrepo, J. Torres, G. Bianconi, Communications Physics 5.1 (2022): 1-17.

(2) ”Local Dirac synchronization on networks” (2023)
L. Calmon, S. Krishnagopal, G. Bianconi, Chaos 33 (3) (2023)


Files: 

ddt_fc.m: function implementing the dynamical equation from paper (1) on fully connected networks.
ddt_II_fc.m: function implementing the dynamical equation from paper (2) on fully connected networks.
ddt_II.m: function implementing the dynamical equation from paper (2) on sparse networks.


Topological_synchronization.m: script that numerically integrates the dynamics defined in paper (1), on a fully connected network -- uses ddt_fc.m.

Dirac_synchronization_fully_connected.m: script that numerically integrates the dynamics defined in paper (2), on a fully connected network -- uses ddt_II_fc.m.

Dirac_synchronization_sparse.m: script that numerically integrates the dynamics defined in paper (2), on a Poisson network -- uses ddt_II.m. Can be edited to include any adjacency matrix.


