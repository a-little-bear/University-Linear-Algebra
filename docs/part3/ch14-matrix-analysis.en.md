# Chapter 14: Matrix Analysis and Convergence

<div class="context-flow" markdown>

**Prerequisites**: Matrix Norms (Ch15) · Matrix Functions (Ch13) · Real Analysis

**Chapter Outline**: Sequences of Matrices → Convergence in Norm → Matrix Series → Geometric Series of Matrices → Gelfand's Formula for Spectral Radius → Norm-independence of Limits → Continuity of Matrix Operations → Matrix Power Series and Stability

**Extension**: Matrix analysis provides the rigorous justification for iterative algorithms and the long-term behavior of discrete dynamical systems.

</div>

Matrix analysis extends the concepts of limits, continuity, and convergence from real numbers to the space of matrices. Because all norms on a finite-dimensional space are equivalent, the choice of matrix norm does not affect whether a sequence converges, only the rate of convergence. The central result of this chapter is **Gelfand's Formula**, which provides the link between the spectral radius and the behavior of high matrix powers. This theory is essential for understanding the stability of iterative methods and the existence of matrix power series.

---

## 14.1 Sequences and Series of Matrices

!!! definition "Definition 14.1 (Convergence of a Sequence)"
    A sequence of matrices $\{A_k\}$ converges to a matrix $A$ (denoted $\lim_{k \to \infty} A_k = A$) if $\lim_{k \to \infty} \|A_k - A\| = 0$ for any matrix norm $\|\cdot\|$.

!!! theorem "Theorem 14.1 (Gelfand's Formula)"
    For any square matrix $A$ and any matrix norm $\|\cdot\|$, the spectral radius $\rho(A)$ satisfies:
    $$\rho(A) = \lim_{k \to \infty} \|A^k\|^{1/k}$$

---

## Exercises

1. **[Fundamentals] Prove that a sequence $\{A_k\}$ converges iff its entries $(A_k)_{ij}$ converge individually.**
   ??? success "Solution"
       This follows from the equivalence of norms. The max-norm $\|A\|_{\max} = \max |a_{ij}|$ is a valid norm, and convergence in this norm is precisely the entry-wise convergence.

2. **[Spectral Radius] If $\rho(A) < 1$, show that $\lim_{k \to \infty} A^k = 0$.**
   ??? success "Solution"
       By Gelfand's formula, for large $k$, $\|A^k\|^{1/k} \approx \rho(A) < 1$. Thus $\|A^k\| \le (\rho(A) + \epsilon)^k$. As $k \to \infty$, the RHS goes to zero, so $A^k \to 0$.

3. **[Geometric Series] When does the series $\sum_{k=0}^\infty A^k$ converge? What is its sum?**
   ??? success "Solution"
       The series converges if and only if $\rho(A) < 1$. In this case, the sum is $(I - A)^{-1}$. This is the matrix version of the scalar geometric series.

4. **[Norm Bound] Prove that $\rho(A) \le \|A\|$ for any consistent matrix norm.**
   ??? success "Solution"
       Let $(\lambda, v)$ be an eigenpair. $Av = \lambda v \implies \|Av\| = \|\lambda v\| = |\lambda| \|v\|$. By the consistency property $\|Av\| \le \|A\| \|v\|$, we have $|\lambda| \|v\| \le \|A\| \|v\|$. Dividing by $\|v\| \neq 0$ gives $|\lambda| \le \|A\|$.

5. **[Invertibility] Show that if $\|A\| < 1$, then $(I-A)$ is non-singular.**
   ??? success "Solution"
       Since $\rho(A) \le \|A\| < 1$, the geometric series $I + A + A^2 + \dots$ converges. Its sum is $(I-A)^{-1}$, confirming that the inverse exists.

6. **[Continuity] Prove that the matrix product $(A, B) \mapsto AB$ is a continuous map.**
   ??? success "Solution"
       $\|A B - A_0 B_0\| = \|A(B-B_0) + (A-A_0)B_0\| \le \|A\|\|B-B_0\| + \|A-A_0\|\|B_0\|$. As $A \to A_0$ and $B \to B_0$, the terms on the right go to zero.

7. **[Exponential] Why does the series $\sum \frac{1}{k!} A^k$ always converge for any $A$?**
   ??? success "Solution"
       Using the ratio test on the norms: $\sum \frac{\|A\|^k}{k!} = e^{\|A\|} < \infty$. Absolute convergence in the norm implies convergence in the matrix space.

8. **[Invariance] Does the choice of norm affect the limit of a sequence?**
   ??? success "Solution"
       No. In finite-dimensional spaces, all norms are equivalent ($c\|A\|_1 \le \|A\|_2 \le C\|A\|_1$). If a sequence converges in one norm, it converges in all others to the same limit.

9. **[Iterative Methods] How does Gelfand's formula relate to the convergence of the power method?**
   ??? success "Solution"
       The power method's error decays as $(\lambda_2/\lambda_1)^k$. Gelfand's formula generalizes this by showing that the dominant eigenvalue (spectral radius) dictates the asymptotic growth or decay of the operator.

10. **[Trace] Is the trace a continuous function?**
    ??? success "Solution"
        Yes. The trace is a linear combination of the entries. Since linear functions on finite-dimensional spaces are always continuous, $A_k \to A \implies \operatorname{tr}(A_k) \to \operatorname{tr}(A)$.

## Chapter Summary

This chapter establishes the analytic foundation for operator limits:

1. **Norm Equivalence**: Demonstrated that convergence in matrix spaces is independent of the coordinate system or metric choice.
2. **Asymptotic Dominance**: Established Gelfand's formula as the law governing the long-term behavior of matrix powers.
3. **Stability Condition**: Positioned the spectral radius condition $\rho(A) < 1$ as the universal requirement for system contraction and series convergence.
4. **Functional Continuity**: Validated matrix operations as continuous mappings, justifying the use of limits in linear algebra.
