<h1>Physically Plausible Animation of 2D Plasma Balls</h1>

Code to accompany bachelor thesis of this title. Parallel computing of the poison equation using the finite element method, implemented in Cuda C and C++. Results written to vtk files.

Requirements: 
- C++14
- Cuda C libraries
- QR solver by simonschoelly
  
<h2>Abstract</h2>
<img width="30%" align="right" alt="Shows an a photograph of a Plasma ball." src="https://user-images.githubusercontent.com/44576195/187530672-32a3c359-6b72-4bdf-b1f3-c2db5a2888ae.jpg"><p>
In this bachelor thesis, filaments are modelled for a virtual plasma ball. These filamentary discharges are of the Dielectric Barrier Discharge (DBD) type. We continue to expand upon and explore a modified Dielectric Breakdown Model (DBM) from the
preceding advanced practical. The modifications draw from existing DBD models to
increase realism in filamentary structure.
We then use aspects of this modified DBM to develop a novel filamentary growth
model. In this model we use a Finite Elements Method (FEM) to solve the Poisson
equation for the electric potential. We then introduce a "virtual electrode" to mitigate
the impact of oscillations that appear in our numerical solution, as as we introduce
filaments and filament growth.
Reintroducing a residual charge between individual executions, allows for the addition
of temporally coherent branching, as well as, a mimicking of the filament life cycle.
The resulting model is implemented by solving the system of linear equations from the
FEM on the GPU to allow for higher resolution modelling before some rudimentary
visualisation</p>
Abstract and images taken from {Thyssen, K. and Sadlo, F., n.d. Physically Plausible Animation of 2D Plasma Balls. Undergraduate Thesis. Ruprecht Karls University Heidelberg.}
  
The full thesis and accompanying presentation is available on request.

<h2>Results (Sneak peek)</h2>

<img width="33%" align="center" alt="Oscillations in electric field." src="https://user-images.githubusercontent.com/44576195/188334042-6095fbad-0fb3-4c2c-b12e-d7f9666426ed.png"><img width="30%" align="center" alt="Branching algorith diagram." src="https://user-images.githubusercontent.com/44576195/188334283-61f9b3f7-8e0a-4483-aec3-9fc5544f0108.png"><img width="33%" align="center" alt="Electric field within a virtual plasma ball as filaments grow." src="https://user-images.githubusercontent.com/44576195/188333509-ba388b07-c0c4-4c17-a565-3c32eddab41e.png">
</br>
<em>(left) Oscillations in electric field due to numerical errrors in poisson field calculation, a difficulty overcome in process (explored in section 4.2.5 Filaments as Boundaries)</em>
<em> (center) Filament growth stepping algorithm diagram (4.2.7 Variation of Step Size)</em>
<em> (right) Electric field within a virtual plasma ball as filaments grow (final) </em>
  
https://user-images.githubusercontent.com/44576195/188333980-7d03ea86-83d5-447e-88fa-2341d9335647.mp4

<em>Video of simulated virtual plasma ball</em>


Results are generated as a series of .vtk files, these can be opened using Paraview or similar programs
