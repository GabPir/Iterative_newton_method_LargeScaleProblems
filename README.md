# Newton Method with Finite differences (from scratch)
## Project of "Numerical methods for Large Scale Problems" Course 2022/2023
### Pt1: low dimensionality case
<br>
This code aims to minimize the Rosenbrock's function starting from two generic initial points:
<br>

![immagine](https://github.com/user-attachments/assets/2dac2735-cc07-4683-9894-cb0b02bf0567)

As a minimization method, the Newton method was applied, where the descent direction is iteratively computed by solving the linear system:

![immagine](https://github.com/user-attachments/assets/f4cb1f94-adcc-4daf-8890-df875e796b60)

The Newton method with exact resolution was compared with the following methods:
- Finite difference approximation: to compute the Hessian and gradient of the second-order local model.
- Hessian matrix correction: in case the matrix is not positive definite (e.g., numerical errors).
- Local search with Armijo's conditions: for determining the descent step.
<br>
The results (in terms of gradient norm) were:

![immagine](https://github.com/user-attachments/assets/b41cc44a-f93b-44e0-a549-d9844dd2fb69)
![immagine](https://github.com/user-attachments/assets/df5f0012-b2ee-4431-b8fa-e0d9e3260747)

<br>
###Pt2: high dimensionality case
<br>
he following functions were considered to test the algorithms with non-trivial problems:

![immagine](https://github.com/user-attachments/assets/786d5ab8-fe86-41ce-932c-06f7d022893f)

The results obtained were:
<br>
![immagine](https://github.com/user-attachments/assets/8392bd44-46ca-484d-bc95-caf7c383922b)
<br>
