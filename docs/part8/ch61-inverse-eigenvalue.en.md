# Chapter 61: Inverse Eigenvalue Problems (IEP)

<div class="context-flow" markdown>

**Prerequisites**: Eigenvalues (Ch06) · Matrix Analysis (Ch14) · Positive Definite Matrices (Ch16) · Numerical Algorithms (Ch22)

**Chapter Outline**: Motivation for Inverse Eigenvalue Problems → Challenges of Existence and Construction → Symmetric IEP (Inverse Application of Schur-Horn) → Jacobi Matrix IEP (Construction from Two Spectra) → Non-negative IEP (The NIEP and Suleimanova Conditions) → Structured IEP (Toeplitz and Banded Constraints) → Additive and Multiplicative IEP → Numerical Solutions: Quasi-Newton and Lifting Methods → Applications: Vibration System Design (Spring-Mass Models), Seismic Imaging, and Pole Placement in Control Theory

**Extension**: If eigenvalue computation is "deriving frequencies from physical structures," then IEP is "designing physical structures based on target frequencies"; it is a core algebraic task in engineering design and parameter identification.

</div>

In much of linear algebra, our task is to compute the eigenvalues of a given matrix. In engineering design, however, the problem is often reversed: we want a system to have specific resonance frequencies (eigenvalues) and need to determine the system's parameters (matrix entries). This is the **Inverse Eigenvalue Problem** (IEP). This chapter explores how to reconstruct matrices satisfying specific structural constraints from fragmentary spectral data.

---

## 61.1 What is an Inverse Eigenvalue Problem?

!!! definition "Definition 61.1 (IEP)"
    Given a set of scalars $\{\lambda_1, \ldots, \lambda_n\}$ and a matrix family $\mathcal{M}$, find a matrix $A \in \mathcal{M}$ such that its eigenvalues are precisely these scalars.
    - **Additive IEP**: Find a diagonal matrix $D$ such that $A+D$ has the given spectrum.
    - **Multiplicative IEP**: Find a diagonal matrix $D$ such that $DA$ has the given spectrum (a variant of the pole placement problem).

---

## 61.2 Typical IEP Types

### 61.2.1 Jacobi Matrix IEP
!!! theorem "Theorem 61.1 (Dual Spectrum Construction)"
    A symmetric tridiagonal matrix (Jacobi matrix) is uniquely determined by its full spectrum and the spectrum of the submatrix obtained by deleting the last row and column. This is the basis for the reconstruction theory of vibrational chains.

### 61.2.2 Non-negative IEP (NIEP)
!!! challenge "The NIEP Puzzle"
    Determining which sets of complex numbers can be the spectrum of some non-negative matrix is a difficult problem involving inverse constraints of Perron-Frobenius theory. It remains a frontier in matrix theory, though the Suleimanova criterion provides some sufficient conditions.

---

## 61.3 Pole Placement: The Control Theory IEP

!!! technique "Pole Placement"
    In a control system $\dot{x} = Ax + Bu$, we use feedback $u = -Kx$ to place the eigenvalues of the closed-loop system $A-BK$ at specific locations in the left half-plane.
    **Algebraic Essence**: Find matrix $K$ such that the characteristic polynomial of $A-BK$ matches a target polynomial.

---

## 61.4 Numerical Solution Methods

!!! algorithm "Algorithm 61.1 (Iterative Lifting Method)"
    1.  Initialize a matrix $A_0 \in \mathcal{M}$.
    2.  Calculate the deviation between current and target eigenvalues.
    3.  Utilize the sensitivity matrix of eigenvalues with respect to entries (the Jacobian, see Ch47B) to correct entries via Newton-like iterations.

---

## Exercises

1.  **[Basics] A $2 \times 2$ symmetric matrix has eigenvalues 1 and 3, and its diagonal entries are both 2. Find the matrix.**
    ??? success "Solution"
        $\operatorname{tr}(A) = 2+2=4 = 1+3$. $\det(A) = 4 - a_{12}^2 = 1 \cdot 3 = 3 \implies a_{12}^2 = 1 \implies a_{12} = \pm 1$.
        The matrix is $\begin{pmatrix} 2 & 1 \\ 1 & 2 \end{pmatrix}$ or $\begin{pmatrix} 2 & -1 \\ -1 & 2 \end{pmatrix}$.

2.  **[Schur-Horn] If target eigenvalues are $(10, 0)$, can the diagonal entries be set to $(6, 4)$?**
    ??? success "Solution"
        Yes, because $(6, 4) \prec (10, 0)$, satisfying the majorization requirement for symmetric matrices.

3.  **[NIEP] Determine if $\{2, -1, -1\}$ can be the spectrum of a non-negative matrix.**
    ??? success "Solution"
        Yes. This set is trace-zero. A $3 \times 3$ circulant matrix like $\begin{pmatrix} 0 & 1 & 1 \\ 1 & 0 & 1 \\ 1 & 1 & 0 \end{pmatrix}$ has eigenvalues $2, -1, -1$.

4.  **[Uniqueness] Why are IEPs typically non-unique?**
    ??? success "Solution"
        Because eigenvalues do not contain information about eigenvectors (rotation). Unless strict structural constraints are applied (like Jacobi matrices), many similar matrices share the same spectrum.

5.  **[Structure] Can an $n \times n$ Toeplitz matrix be uniquely determined by its eigenvalues?**
    ??? success "Solution"
        Generally no. Eigenvalues provide only $n$ pieces of information, while a Toeplitz matrix has $2n-1$ parameters.

6.  **[Vibration] Briefly describe the application of IEP in bridge monitoring.**
    ??? success "Solution"
        By measuring the vibration frequencies (eigenvalues) of a bridge under excitation, one can back-calculate changes in the stiffness matrix to locate structural damage.

7.  **[Polynomial] How do you write the target characteristic polynomial for a target spectrum $\{\lambda_1, \lambda_2\}$?**
    ??? success "Solution"
        $p(\lambda) = (\lambda - \lambda_1)(\lambda - \lambda_2) = \lambda^2 - (\lambda_1+\lambda_2)\lambda + \lambda_1\lambda_2$.

8.  **[Companion] Is the companion matrix a trivial construction for an IEP?**
    ??? success "Solution"
        Yes. It fills in the coefficients of the target polynomial directly, providing the simplest non-symmetric construction for any given spectrum.

9.  **[Sensitivity] Why do repeated eigenvalues make IEPs numerically difficult?**
    ??? success "Solution"
        Eigenvalues are not Fréchet differentiable at points of multiplicity (see Ch47B), making gradient-based methods fail or become extremely unstable.

10. **[Control] If system $(A, B)$ is not controllable, can arbitrary pole placement be achieved?**
    ??? success "Solution"
        No. Only the poles of the controllable part can be moved; uncontrollable poles remain fixed.

## Chapter Summary

Inverse eigenvalue problems are the "blueprints" of mathematical modeling:

1.  **From Abstract to Real**: IEP establishes the path from abstract spectral targets back to concrete physical parameters, bridging theoretical physics and engineering.
2.  **Art of Constraint**: The core of IEP lies in the balance between "structural constraints" and "spectral data"; an unconstrained problem is trivial (companion matrix), while an over-constrained one is unsolvable.
3.  **Return to Stability**: Through numerical IEP algorithms, we can finely adjust local elements of an operator to achieve global stability, demonstrating linear algebra's agency in complex system optimization.
