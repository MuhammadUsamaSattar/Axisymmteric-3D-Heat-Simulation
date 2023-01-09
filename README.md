# Axisymmteric-3D-Heat-Simulation (Computational Heat Transfer)
Simulation of heat transfer in an axisymmetric mesh using Finite Difference Method.
Pre-defined cases are solved with the code to demonstrate usability of FDM for solution of PDES. Both implicit and explicit schemes have been employed to generate results. The boundary condition is assumed to be that the outer wall of the cylinder is at 35 C.
A case for convective heat transfer is also simualted with air temperature of 25 C and convective heat transfer of 15 W/Km2.

Contour plots are generated to show temperature at each node at different times. Additionally, a plot is generated for temperature variation as measured by a thermocouple placed at the center of the axisymmetric mesh at 0.25m above the bottom surface.

# Results
### Implicit
![Implicit at 86400s](https://user-images.githubusercontent.com/60822455/208429636-61c25c82-94fc-4289-a09f-8fa5f57e781e.png)
![Implicit Thermocouple](https://user-images.githubusercontent.com/60822455/208429965-9ed7e1b6-6b74-48d7-a306-3c5ad8e0478b.png)

### Explicit
![Explicit at 86400s](https://user-images.githubusercontent.com/60822455/208429683-7642e54d-2b9d-443c-884f-e836c328f742.png)
![Explict Thermocouple](https://user-images.githubusercontent.com/60822455/208430015-59a4b948-a7c9-4120-8ee7-ce8ac1e7b064.png)

### Convection
![Convection at 86400s](https://user-images.githubusercontent.com/60822455/208429713-06f41fd9-356a-493e-8eda-59a3fb61c23c.png)
![Convection Thermocouple](https://user-images.githubusercontent.com/60822455/208430032-ef9e0e87-fbe2-46bd-adb8-b8ae6df13f3f.png)
