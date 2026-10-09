(1.0) About the file  "MDH.f "
The fortran code, entitled "MDH.f " , uses the boundary element method(BEM) and molecular Debye-Huckel(MDH) theory to calculate the electrostatic energy of a  diatomic solute solvated  in an electrolyte solution.
The electrostatic response of the electrolyte solution is described by two DH-like response modes.
The two effective inverse Debye length of the solution are  k1=(0.5846+0.7673*i)/Angstrom and  k2=(0.5846-0.7673*i)/Angstrom. 
k1 and k2 are a pair of  conjugate complex numbers.

(1a)The  "MDH.f" file directly evaluate the electrostatic energy contribbution  of  one DH mode with k1=(0.5846+0.7673*i)/Angstrom . 
This electrostatic energy contribution is  a complex number, e.g.,  u1= c1 +c2*i. 
 The second DH mode with  effective Debye parameter k2=(0.5846-0.7673*i)/Angstrom would lead to
an electrostatic energy u2=  c1 -c2*i. Then the total electrostatic energy from MDH theory with two DH response modes  reads ue= u1+ u2= 2*c1= 2*Re(u1). 

The code "MDH.f" needs to be compiled  together with the LAPACK package to obtain an excutive file.

(1b)Four input files are needed for the BEM electrostatic energy calculation using "MDH.f " : 
the file "DA4.crg" contains the location and  charge number of each atom of the complex solute
the file "DA4.vert" contains the vertice information of the molecular surface
the file "DA4.face" contains the face information of the molecular surface
the file "DA4.curvature" contains the curvature information of the molecular surface.

(1c) Both the "DA4.vert" file and "DA4.face" file are generated with Sanner's msms code(e.g., https://csb.vanderbilt.edu/comp/soft/msms/) .
"DA4.xyzr" is the input file for msms code.

The file "DA4.curvature" is generated with the fortran code "curvature.f90".

(1e)For a solute with fixed solute geommetry but tunable site charges ,  the  "DA4.vert", "DA4.face"  and  "DA4.curvature"  files keep the same 
while the  "DA4.crg" needs to be modified according to the site charge numbers.
e.g., the "DA4.crg" file is used for a diatomic solute with site charges (z1=-1, z2=0).
The "DA4b.crg" file is used for a diatomic solute with site charges (z1=-1, z2=1). 
When the code "MDH.f" or "DH.f" is used to study a diatomic solute with site charges (z1=-1, z2=1),  please copy the content of "DA4b.crg" file to  "DA4.crg"
or modify the site charge numbers of "DA4.crg".

(1.1) About the file  "DH.f "
The fortran code, entitled "DH.f " , uses the boundary element method(BEM) to calculate the electrostatic energy of a  solute solvated  in an electrolyte solution,
for a DH response mode with Debye parameter kD= 0.805 /Angstrom.
Only three input files ("DA4.crg", "DA4.vert", "DA4.face" ) are  needed for the electrostatic energy calculation.

(2) About the file  "curvature.f90"
 The fortran code, entitled "curvature.f90", is used to calculate the Gaussian curvature and mean curvature of the molecular surface of a diatomic solute.

Three input files ("DA4.xyzr", "DA4.vert" and "DA4.face" ) are needed for the "curvature.f90" code.
The output file "DA4.curvature" contains the curvature information( Gaussian curvature and mean curvature ) on the molecular surface. 

(3) About Sanner's msms code
Both the "filename.vert" file and "filename.face" file are generated with Sanner's msms code(e.g., https://csb.vanderbilt.edu/comp/soft/msms/) 
using the "filename.xyzr" file as input.
The  "filename.xyzr" file, which contains the location and hard sphere radius of the atoms of the complex solute.
e.g., for probe radius 1.00 Angstrom and vertice number density 2.0 vertice/Angstrom^2, 
one may use the msms command "./msms -if DA4.xyzr -prob 1.00 -den 2.0 -of DA4" to generate  "DA4.vert" and "DA4.face" with "DA4.xyzr" as input.

If one want to use the above code to calculate the electrostatic energy of other solutes with different geometry and site charges,
then one needs to update the related input files.








