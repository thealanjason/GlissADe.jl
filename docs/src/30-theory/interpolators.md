# Interpolation Schemes

*When you need values at faces in a cell-centred solver.*

## Introduction

The current implementation of the finite area method is cell centered. Even though it has some benefits in the discretization of the non-linear equations of the Savage-Hutter model, in some ways it has its disadvantages. One of them is the calculation of face values. The presence of face values in the discretization of equations is what allows for coupling of cells, which then results in accurate representation of flow physics. Without the coupling of cells, all cells are independent and no flow ever occurs.

Popular schemes from commercial solvers like OpenFOAM, Ansys Fluent, and Star CCM use different kinds of interpolation schemes, such as:

- Linear/Central Differencing
- Upwind Differencing
- Linear Upwind Differencing
- Hybrid Schemes
- QUICK

QUICK is a quadratic interpolation scheme. Even though it's more accurate, it is much less preferred because of its non-linearity. Here's a brief description of each of the interpolation schemes available in the library.

## Linear/Central Differencing

This might be the most commonly used interpolation scheme. Since it is a linear scheme and depends on both the cell ``P`` and its neighbouring cell ``N``, it is a second order scheme. This scheme is preferred because of its high accuracy, but it has some disadvantages: it is unbounded and causes unphysical numerical oscillations. Thus, this scheme is generally used for interpolating diffusion terms and often avoided for discretizing convection terms.

The central differencing interpolation scheme is simply:

```math
\phi_f = f_x\;\phi_N + (1-f_x)\;\phi_P \quad , \quad  f_x = \frac{\lvert x_f - x_P \rvert}{\lvert x_N - x_P \rvert} \quad , \quad 0 \leq f_x \geq 1
```

This scheme can't identify the direction of flux (flow), making it less reliable for simulating convection-dominated flows.

## Upwind Differencing

This scheme is dependent on the direction of flux. Let ``F_f`` be the mass flow rate through a face, expressed as:

```math
F_f = \rho_f\;A_f\;(U_f \cdot \hat{n})
```

where ``\hat{n}`` is the unit normal pointing away from the face. Then,

```math
F_f =  \begin{cases}
      \geq 0 & \text{mass flow out of cell} \\
      >0 & \text{mass flow into cell} \\
   \end{cases}
```

and the interpolation scheme is:

```math
\phi_f = \begin{cases}
          \phi_P & F_f \geq 0 \\
          \phi_N & F_f < 0 \\
          \end{cases}
```

Since this scheme doesn't vary linearly between cells, it is first order and can cause truncation errors. Even though it is less accurate, this scheme is stable, and is generally used to generate an initial stable solution which is then corrected using higher order interpolations.

## Linear Upwind Differencing

An extension to the upwind differencing scheme, this scheme lies between accuracy and stability. It also depends on the direction of the flux, as:

```math
\phi_f = \begin{cases}
          \phi_P + \nabla\phi_P \cdot r & F_f \geq 0 \\
          \phi_N + \nabla\phi_N \cdot r & F_f < 0 \\
          \end{cases}
```

Even though this scheme depends somewhat linearly between cells and their neighbours, it is second order. But in cases of a high gradient, there can be situations where ``\phi_P \leq \phi_f \leq \phi_N`` is violated, causing numerical oscillations. To avoid this problem, gradient limiters are introduced:

```math
\phi_f = \begin{cases}
          \phi_P + \lambda\nabla\phi_P \cdot r & F_f \geq 0 \\
          \phi_N + \lambda\nabla\phi_N \cdot r & F_f < 0 \\
          \end{cases}
```

Notice the ``\lambda`` in the equations. It lies in the range ``\lbrack 0, 1 \rbrack`` and controls how much influence the gradient has on the face value. The value of ``\lambda`` is typically chosen such that the following holds:

```math
\phi_P \leq \phi_f \leq \phi_N
```

if ``P \rightarrow N`` is the direction of the flux. In case the direction of flux is ``N\rightarrow P``, the inequalities are reversed.

## Hybrid Schemes

Hybrid schemes are nothing more than the linear combination of the primitive interpolation schemes described above. They can be expressed as:

```math
\phi_f = \Psi\;\phi_{UD} + (1-\Psi)\;\phi_{CD}
```

where ``\Psi`` is the chosen blending function. This scheme tries to retain the accuracy of the central differencing scheme and the stability of the upwind scheme.

The blending function ``\Psi`` is generally chosen such that:

- Near areas of high numerical oscillation,

```math
\Psi \rightarrow 1
```

  i.e. use upwind differencing

- Near areas of less numerical oscillation,

```math
\Psi \rightarrow 0
```

  i.e. use central differencing
