# Jacobians and the DCM Derivation

## A note on jacobians during integral transformations

Jacobians arise wherever there are multidimensional integral transformations. In the current case, transformations are done for depth averaging, i.e. the following splitting of integrals:

```math
\int_V = \int_S \int_0^{h} \lvert J \rvert
```

Notice the requirement to compute the Jacobian ``J``. However, for surfaces with mild curvatures, the discrete elements are planar. The planar geometry of such meshes has an added advantage: both the global and local coordinate systems are Cartesian. This means that the Jacobian used for coordinate systems, expressed as:

```math
J = \begin{bmatrix} \frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \frac{\partial f_1}{\partial x_3} \\ \frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \frac{\partial f_2}{\partial x_3} \\ \frac{\partial f_3}{\partial x_1} & \frac{\partial f_3}{\partial x_2} & \frac{\partial f_3}{\partial x_3} \end{bmatrix}
```

can be solved analytically. Here ``\mathbf{f}`` and ``\mathbf{x}`` represent the local and global coordinate systems. Since both the global and local coordinate systems are of the same type, they can be expressed as:

```math
f_1 = J_{11} x_1 + J_{12} x_2 + J_{13} x_3
```

```math
f_2 = J_{21} x_1 + J_{22} x_2 + J_{23} x_3 \\
```

```math
f_3 = J_{31} x_1 + J_{32} x_2 + J_{33} x_3
```

or more conveniently,

```math
    \begin{bmatrix} f_1 \\ f_2 \\ f_3 \end{bmatrix} = \begin{bmatrix} J_{11} & J_{12} & J_{13} \\ J_{21} & J_{22} & J_{23} \\ J_{31} & J_{32} & J_{33} \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \\ x_3 \end{bmatrix}
```

where ``\mathbf{J}`` is the Jacobian. But how do we calculate the Jacobian? We only have the numerical values of the coordinate vectors, not their variable representations. Here's where the Cartesian nature of both coordinate systems comes in handy. Because both the local and global coordinate systems are Cartesian, the local coordinate system can be represented as a 3D rotation of the global system, expressed simply as:

```math
f = R x
```

where ``\mathbf{R}`` is expressed as:

```math
\mathbf{R} = \begin{bmatrix} f_1 \cdot x_1 & f_1 \cdot x_2 & f_1 \cdot x_3 \\ f_2 \cdot x_1 & f_2 \cdot x_2 & f_2 \cdot x_3 \\ f_3 \cdot x_1 & f_3 \cdot x_2 & f_3 \cdot x_3 \end{bmatrix}
```

the direction cosine matrix (DCM). Since DCMs are orthogonal, their determinant is ``1``. This means that the Jacobian of the coordinate transformation also has determinant ``1``, which lets us ignore the transformation determinants when transforming integrals.
