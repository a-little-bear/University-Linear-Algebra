# Chapter 51: Quaternion Matrices

<div class="context-flow" markdown>

**Prerequisites**: Matrix Algebra (Ch02) · Eigenvalues (Ch06) · Clifford Algebra (Ch50) · Matrix Groups (Ch55)

**Chapter Outline**: Algebraic Structure of Quaternions $\mathbb{H}$ → Non-commutativity and its Challenges → Complex Representation $\chi(A)$ → Left vs. Right Eigenvalues → Diagonalization & SVD of Quaternion Matrices → Unitary Quaternions & the Symplectic Group $Sp(n)$ → Applications: 3D Rotation (Avoiding Gimbal Lock), Signal Processing, and Quaternion Neural Networks

**Extension**: Quaternion matrices provide the optimal tool for describing 3D rotations and attitude; they extend complex linear algebra into the non-commutative realm and are key to understanding Symplectic Geometry and Compact Lie Groups.

</div>

Complex numbers brought linear algebra into the realm of 2D rotations, but **Quaternions** ($\mathbb{H}$) extend this power to 3D space. Because quaternion multiplication is non-commutative ($ij \neq ji$), the theory of quaternion matrices exhibits a landscape starkly different from classical linear algebra. This chapter establishes rigorous methods for handling such non-commutative algebraic structures and demonstrates their immense power in computer graphics and attitude control.

---

## 51.1 Foundations of Quaternions $\mathbb{H}$

!!! definition "Definition 51.1 (Quaternions)"
    A quaternion $q$ is expressed in the form $q = a + bi + cj + dk$, where $a, b, c, d \in \mathbb{R}$.
    The fundamental imaginary units satisfy:
    $$i^2 = j^2 = k^2 = ijk = -1$$
    **Core Property**: $ij = k, ji = -k$ (Non-commutativity).

!!! definition "Definition 51.2 (Conjugate and Norm)"
    - **Conjugate**: $\bar{q} = a - bi - cj - dk$.
    - **Norm**: $|q| = \sqrt{q\bar{q}} = \sqrt{a^2 + b^2 + c^2 + d^2}$.
    - **Inverse**: $q^{-1} = \bar{q} / |q|^2$ (if $q \neq 0$).

---

## 51.2 Quaternion Matrices and Complex Representation

!!! definition "Definition 51.3 (Complex Representation)"
    Any quaternion $q = \alpha + \beta j$ ($\alpha, \beta \in \mathbb{C}$) can be represented as a $2 \times 2$ complex matrix:
    $$\chi(q) = \begin{pmatrix} \alpha & \beta \\ -\bar{\beta} & \bar{\alpha} \end{pmatrix}$$
    For an $n \times n$ quaternion matrix $A = A_1 + A_2 j$, its **complex adjoint matrix** $\chi(A)$ is a $2n \times 2n$ complex matrix:
    $$\chi(A) = \begin{pmatrix} A_1 & A_2 \\ -\bar{A}_2 & \bar{A}_1 \end{pmatrix}$$
    **Significance**: This allows us to transform non-commutative quaternion operations into well-established complex matrix operations.

---

## 51.3 Eigenvalues: Left and Right

!!! note "Warning: The Cost of Non-commutativity"
    Since $Aq = \lambda q$ and $Aq = q\lambda$ are not equivalent over $\mathbb{H}$, quaternion matrices have two types of eigenvalues:
    1.  **Left Eigenvalues**: satisfy $Ax = \lambda x$. These lack good properties and are rarely used.
    2.  **Right Eigenvalues**: satisfy $Ax = x\lambda$. These are the focus of research and possess clear geometric meaning.

!!! theorem "Theorem 51.1 (Properties of Right Eigenvalues)"
    An $n \times n$ quaternion matrix $A$ has exactly $n$ right eigenvalues (in the sense of equivalence classes). If $\lambda$ is a right eigenvalue, then for any non-zero quaternion $u$, $u^{-1}\lambda u$ is also a right eigenvalue (forming an eigenvalue orbit).

---

## Exercises

1.  **[Basics] Calculate the value of $ij + ji$.**
    ??? success "Solution"
        $k + (-k) = 0$. This reflects the anti-commutativity of the imaginary units.

2.  **[Conjugate] Prove $\overline{q_1 q_2} = \bar{q}_2 \bar{q}_1$.**
    ??? success "Solution"
        Expand using the definition. Note that due to non-commutativity, the order of multiplication must be reversed under conjugation.

3.  **[Representation] Write the $2 \times 2$ complex representation matrix for the unit $j$.**
    ??? success "Solution"
        $j = 0 + 1j \implies \alpha=0, \beta=1$. Thus $\chi(j) = \begin{pmatrix} 0 & 1 \\ -1 & 0 \end{pmatrix}$.

4.  **[Rotation] What rotation does a unit quaternion $q = \cos(\theta/2) + \mathbf{u} \sin(\theta/2)$ represent?**
    ??? success "Solution"
        It represents a rotation by angle $\theta$ around the unit axis $\mathbf{u}$.

5.  **[Eigenvalue] Prove: If $\lambda$ is a right eigenvalue of $A$, then its conjugate $\bar{\lambda}$ must belong to an eigenvalue orbit of $A$.**
    ??? success "Solution"
        Utilizing the complex adjoint matrix $\chi(A)$. The eigenvalues of $\chi(A)$ appear in conjugate pairs, which correspond to the right eigenvalue orbits of the original quaternion matrix.

6.  **[Unitary] What is a unitary quaternion matrix (or Symplectic matrix)?**
    ??? success "Solution"
        A quaternion matrix satisfying $A^* A = I$. The group of such matrices is known as the compact symplectic group $Sp(n)$.

7.  **[Determinant] Why is there no simple determinant definition for quaternion matrices?**
    ??? success "Solution"
        Because multiplication is non-commutative, the order of terms in an expansion cannot be uniquely determined. Typically, the determinant of the real representation or the Dieudonné determinant is used.

8.  **[Inverse] Find the inverse of $q = 1 + i$.**
    ??? success "Solution"
        $|q|^2 = 1^2 + 1^2 = 2$. Thus $q^{-1} = (1-i)/2 = 0.5 - 0.5i$.

9.  **[Application] Why are quaternions preferred over Euler angles in satellite attitude control?**
    ??? success "Solution"
        Because quaternions provide a continuous coverage of the sphere, avoiding the singularities (Gimbal lock) associated with Euler angles at specific orientations.

****

??? success "Solution"
    

## Chapter Summary

Quaternion matrices are powerful tools for handling high-dimensional rotation and symmetry:

1.  **Dimensional Leap**: By sacrificing the commutative law, quaternions gain the ability to perform compact rotation arithmetic in 3D and 4D spaces, serving as the logical extension of complex numbers in geometric descriptive power.
2.  **Bridging via Adjoints**: Complex adjoint matrices transform difficult non-commutative computations into established linear algebra algorithms, ensuring numerical stability for quaternion calculus.
3.  **Group Foundations**: The theory of quaternion matrices is the natural gateway to Symplectic Geometry and Lie groups (e.g., $SU(2) \cong Sp(1)$), revealing deep connections between operator algebra and physical spatial symmetries.
