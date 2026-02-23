# Chapter 46A: Operator Monotone and Operator Convex Functions

<div class="context-flow" markdown>

**Prerequisites**: Positive Definite Matrices (Ch16) · Matrix Functions (Ch13) · Matrix Inequalities (Ch18)

**Chapter Outline**: Löwner Partial Order → Definition of Operator Monotone Functions → Löwner-Heinz Inequality → Furuta Inequality → Löwner's Theorem (Pick Functions & Integral Representations) → Löwner Matrices → Operator Convex and Concave Functions → Hansen-Pedersen Characterization → Jensen's Operator Inequality → Lieb's Concavity Theorem → Strong Subadditivity of Quantum Entropy

**Extension**: Operator monotone functions are the mathematical foundation for quantum relative entropy and fidelity in quantum information science.

</div>

In the scalar world, $a \ge b > 0 \Rightarrow a^2 \ge b^2$. However, for matrices, $A \succeq B \succ 0$ does **not** imply $A^2 \succeq B^2$. This observation leads to the study of **operator monotone functions**—functions that preserve the Löwner partial order across all dimensions. Löwner's Theorem (1934) provides a complete characterization of these functions through their analytic continuation to the upper half-plane.

---

## 46A.1 Core Concepts

!!! definition "Definition 46A.2 (Operator Monotone)"
    A continuous function $f: I 	o \mathbb{R}$ is **operator monotone** if for all $n$ and all Hermitian matrices $A, B$ with spectra in $I$:
    $$A \succeq B \quad \Rightarrow \quad f(A) \succeq f(B).$$

!!! theorem "Theorem 46A.3 (Löwner-Heinz Inequality)"
    The power function $f(t) = t^r$ is operator monotone on $[0, \infty)$ if and only if $0 \le r \le 1$.

---

## Exercises

****
??? success "Solution"
         Note that $B^{-1} - A^{-1} = B^{-1}(A - B)A^{-1}$. Let $C = A-B \succeq 0$.
    $B^{-1} - A^{-1} = B^{-1} C A^{-1}$. This product is not obviously PSD.
    However, $B \preceq A \Rightarrow A^{-1/2} B A^{-1/2} \preceq I \Rightarrow \lambda_{\max}(A^{-1/2} B A^{-1/2}) \le 1$.
    The eigenvalues of $A^{1/2} B^{-1} A^{1/2}$ are the inverses, so $\lambda_{\min} \ge 1$, which means $A^{1/2} B^{-1} A^{1/2} \succeq I$.
    Multiplying by $A^{-1/2}$ on both sides gives $B^{-1} \succeq A^{-1}$.

****
??? success "Solution"
         For $z = -1 + i$ (in the upper half-plane), $z^2 = (-1+i)^2 = 1 - 1 - 2i = -2i$.
    The imaginary part is $-2 < 0$. Thus, $t^2$ is not a Pick function.

****
??? success "Solution"
         $L_{22} = f'(4) = \frac{1}{2\sqrt{4}} = 0.25$.
    $L_{12} = L_{21} = \frac{\sqrt{4}-\sqrt{1}}{4-1} = \frac{2-1}{3} = 1/3$.
    $L_2 = \begin{pmatrix} 0.5 & 1/3 \ 1/3 & 0.25 \end{pmatrix}$. $\det L_2 = 0.125 - 1/9 > 0$. It is positive definite.

****
??? success "Solution"
         $(\lambda A + (1-\lambda)B)^2 = \lambda^2 A^2 + (1-\lambda)^2 B^2 + \lambda(1-\lambda)(AB + BA)$.
    The inequality reduces to checking if $\lambda(1-\lambda)(A-B)^2 \succeq 0$, which is true since $(A-B)^2 \succeq 0$ for any Hermitian $A, B$.

****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ****
??? success "Solution"
    ## Chapter Summary

This chapter explores the preservation of matrix structure under functional maps:


****: Defined operator monotonicity and identified the barrier at $t^r$ for $r > 1$.

****: Leveraged Löwner's theorem to link operator monotonicity to the theory of Pick functions.

****: Developed operator convexity and Jensen-type inequalities for matrix operators.

****: Applied these results to establish the joint concavity of trace functions, leading to the proof of Strong Subadditivity in quantum mechanics.
