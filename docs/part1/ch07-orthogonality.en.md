# Chapter 07: Orthogonality and Least Squares

<div class="context-flow" markdown>

**Prerequisites**: Inner Product Basics · Vector Spaces (Ch4) · Systems of Equations (Ch1)

**Chapter Outline**: Inner Product and Norm → Definition of Orthogonality → Orthogonal Bases and Matrices → Gram-Schmidt Process → QR Decomposition → Orthogonal Projection → Least Squares Method → Normal Equations

**Extension**: Orthogonality ensures the stability of numerical computations and is the geometric soul of signal processing (Fourier Transform).

</div>

Orthogonality introduces the concepts of "angle" and "distance" to linear spaces. Among all bases, orthogonal bases are favored for their computational simplicity and numerical stability. The Least Squares method is the optimal strategy for handling "unsolvable" systems in the real world.

---

## 07.1 Orthogonality and Projection

!!! definition "Definition 07.1 (Orthogonal Vectors)"
    Two vectors $u, v$ are orthogonal if their inner product is zero: $\langle u, v \rangle = u^T v = 0$.

!!! theorem "Theorem 07.5 (Normal Equations)"
    The optimal solution to the linear least squares problem $\min \|Ax - b\|^2$ satisfies the normal equations:
    $$A^T A \hat{x} = A^T b$$

---

## Exercises

1. **[Fundamentals] Verify if $v_1 = (1, 1)^T$ and $v_2 = (1, -1)^T$ are orthogonal.**
   ??? success "Solution"
       $v_1 \cdot v_2 = 1(1) + 1(-1) = 1 - 1 = 0$. Yes, they are orthogonal.

2. **[Normalization] Normalize the vector $v = (3, 4)^T$.**
   ??? success "Solution"
       Norm $\|v\| = \sqrt{3^2 + 4^2} = 5$.
       Unit vector $u = v/5 = (0.6, 0.8)^T$.

3. **[Orthogonal Matrix] If $Q$ is an orthogonal matrix, prove $\|Qx\| = \|x\|$.**
   ??? success "Solution"
       $\|Qx\|^2 = (Qx)^T (Qx) = x^T Q^T Q x$.
       Since $Q$ is orthogonal, $Q^T Q = I$.
       Thus $\|Qx\|^2 = x^T I x = x^T x = \|x\|^2$. Orthogonal transformations preserve length.

4. **[Projection] Find the projection of vector $b = (1, 2, 3)^T$ onto the direction of $a = (1, 1, 1)^T$.**
   ??? success "Solution"
       Projection formula $p = \frac{a \cdot b}{a \cdot a} a$.
       $a \cdot b = 1+2+3 = 6$; $a \cdot a = 1+1+1 = 3$.
       Thus $p = (6/3)a = 2a = (2, 2, 2)^T$.

5. **[Gram-Schmidt] Orthogonalize $v_1 = (1, 1, 0), v_2 = (1, 0, 1)$.**
   ??? success "Solution"
       1. $u_1 = v_1 = (1, 1, 0)$.
       2. $u_2 = v_2 - \frac{v_2 \cdot u_1}{u_1 \cdot u_1} u_1 = (1, 0, 1) - \frac{1}{2}(1, 1, 0) = (0.5, -0.5, 1)$.
       The orthogonal basis is $\{(1, 1, 0), (0.5, -0.5, 1)\}$.

6. **[Least Squares] Given $Ax=b$ has no solution, with $A = \begin{pmatrix} 1 \\ 1 \end{pmatrix}, b = \begin{pmatrix} 1 \\ 3 \end{pmatrix}$. Find the least squares solution $\hat{x}$.**
   ??? success "Solution"
       Normal equation $A^T A \hat{x} = A^T b$.
       $A^T A = (1, 1) \begin{pmatrix} 1 \\ 1 \end{pmatrix} = 2$.
       $A^T b = (1, 1) \begin{pmatrix} 1 \\ 3 \end{pmatrix} = 4$.
       Thus $2\hat{x} = 4 \implies \hat{x} = 2$.

7. **[QR Intro] If $A = QR$ ($Q$ orthogonal, $R$ upper triangular), prove $A^T A = R^T R$.**
   ??? success "Solution"
       $A^T A = (QR)^T (QR) = R^T Q^T Q R = R^T I R = R^T R$. This shows how QR decomposition simplifies normal equations.

8. **[Pythagorean] Prove: If $u \perp v$, then $\|u+v\|^2 = \|u\|^2 + \|v\|^2$.**
   ??? success "Solution"
       $\|u+v\|^2 = (u+v) \cdot (u+v) = u \cdot u + 2(u \cdot v) + v \cdot v$.
       Since $u \cdot v = 0$, the expression reduces to $\|u\|^2 + \|v\|^2$.

9. **[Orthogonal Complement] If $W$ is the $xy$-plane in $\mathbb{R}^3$, find its orthogonal complement $W^\perp$.**
   ??? success "Solution"
       $W^\perp$ is the set of all vectors perpendicular to the $xy$-plane, which is the $z$-axis.
       $W^\perp = \{ (0, 0, z)^T : z \in \mathbb{R} \}$.

10. **[Error Vector] Prove: The least squares error vector $e = b - A\hat{x}$ is orthogonal to the column space of $A$.**
    ??? success "Solution"
        $A^T e = A^T(b - A\hat{x}) = A^Tb - A^TA\hat{x}$.
        From the normal equation $A^TA\hat{x} = A^Tb$, so $A^T e = 0$.
        This means $e$ is orthogonal to every column of $A$.

## Chapter Summary

Orthogonality is the "vertical" aesthetic of linear algebra:

1. **Projection Logic**: Least Squares is essentially finding the orthogonal projection of a target vector onto the feasible subspace.
2. **Decomposition Power**: QR decomposition transforms a general matrix into a combination of orthogonal rotation and scaling.
3. **Geometric Metric**: The inner product establishes a rigid structure in space, allowing us to quantify approximation errors.
