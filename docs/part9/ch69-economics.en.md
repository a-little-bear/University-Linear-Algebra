# Chapter 69: Applications of Linear Algebra in Economics

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch1) · Matrix Operations (Ch2) · Eigenvalues (Ch6) · Non-negative Matrices (Ch17) · M-Matrices (Ch38)

**Chapter Outline**: Leontief Input-Output Model $\to$ Hawkins-Simon Conditions $\to$ Sraffa Model $\to$ Input-Output Structure $\to$ Linear Exchange Model $\to$ Game Theory $\to$ Von Neumann Minimax Theorem $\to$ Linear Programming Duality

**Extension**: Input-output analysis (Leontief Nobel Prize) is a standard tool for macroeconomic policy formulation; matrix methods in game theory form the basis for market equilibrium analysis.

</div>

A core problem in economics is understanding the interdependence between sectors of an economy. Leontief's input-output model abstracts this dependence as a consumption matrix, revealing the algebraic essence of production activities.

---

## 69.1 Input-Output Theory

!!! definition "Definition 69.1 (Leontief Model)"
    Let $C$ be the consumption matrix, $x$ the total output, and $d$ the external demand. The equilibrium equation is:
    $$x = Cx + d \implies (I-C)x = d$$
    where $(I-C)^{-1}$ is known as the Leontief inverse matrix, representing the multiplier effect of the economy.

!!! theorem "Theorem 69.3 (Hawkins-Simon Conditions)"
    For any non-negative demand $d \ge 0$, a non-negative output $x \ge 0$ exists if and only if all leading principal minors of $I-C$ are positive.

---

## Exercises

1. **[Fundamentals] Why must the spectral radius of the consumption matrix $C$ in the Leontief model be less than 1?**
   ??? success "Solution"
       Only if $\rho(C) < 1$ can the series $(I-C)^{-1} = I + C + C^2 + \dots$ converge and yield a non-negative result. Economically, this represents a "productive" system where the internal resources required to produce 1 unit of output are strictly less than 1 unit.

2. **[Calculation] Given $C = \begin{pmatrix} 0.5 & 0.4 \\ 0.2 & 0.5 \end{pmatrix}$. Verify the Hawkins-Simon conditions.**
   ??? success "Solution"
       $I-C = \begin{pmatrix} 0.5 & -0.4 \\ -0.2 & 0.5 \end{pmatrix}$.
       Leading principal minors: $D_1 = 0.5 > 0$; $D_2 = \det(I-C) = 0.25 - 0.08 = 0.17 > 0$. Both are positive, so the conditions are satisfied.

3. **[Economic Intuition] Explain the meaning of the elements in the $j$-th column of the Leontief inverse matrix $(I-C)^{-1}$.**
   ??? success "Solution"
       The $j$-th column describes the necessary increase in total output across all sectors of the economy required to satisfy a 1-unit increase in the final demand for the $j$-th sector's products.

4. **[Game Theory] Prove: In a two-person zero-sum game, if the payoff matrix $A$ is skew-symmetric ($A^T = -A$), the value of the game is 0.**
   ??? success "Solution"
       By the minimax theorem, $v = \max_x \min_y x^T A y$. If $A$ is skew-symmetric, then $x^T A x = 0$ for any strategy $x$. Since the game is perfectly symmetric, neither player can have a non-zero expected gain in equilibrium.

5. **[Shadow Prices] In linear programming, what does a positive shadow price for a particular resource constraint imply?**
   ??? success "Solution"
       It implies that the resource is a bottleneck. Increasing the capacity of that resource by 1 unit will increase the objective value (e.g., total profit) by the amount of the shadow price.

6. **[Income Mobility] How is "class stagnation" measured using the eigenvalues of an income mobility matrix $P$?**
   ??? success "Solution"
       The rate of convergence to the steady-state distribution is determined by $|\lambda_2|^k$. As $\lambda_2 \to 1$, the time required to reshuffle the initial distribution increases, representing higher levels of socioeconomic rigidity.

7. **[Sraffa Model] Why are wage rates and profit rates negatively correlated in the Sraffa price model?**
   ??? success "Solution"
       Because the total value produced is fixed by the technology matrix $A$. This value must be split between labor (wages) and capital (profits). The equation $\mathbf{p}^T[I - (1+r)A] = w\mathbf{l}^T$ strictly defines this tradeoff.

8. **[Calculation] Find the mixed strategy equilibrium for $A = \begin{pmatrix} 1 & -1 \\ -1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       Due to symmetry, the optimal strategies are $x^* = [0.5, 0.5]^T$ and $y^* = [0.5, 0.5]^T$. The value of the game is $0.5(1) + 0.5(-1) = 0$.

9. **[Perron-Frobenius] What ensures the uniqueness of the equilibrium output ratio in a closed economic model?**
   ??? success "Solution"
       The Perron-Frobenius theorem for irreducible stochastic matrices guarantees that the eigenvector corresponding to the eigenvalue 1 is unique up to scaling and is strictly positive.

10. **[Structural Analysis] In matrix terms, what signifies an economic crisis in an input-output system?**
    ??? success "Solution"
        A crisis is signified by $\rho(C) \ge 1$ or the failure of the Hawkins-Simon conditions. Algebraically, this results in the Leontief inverse having negative entries or diverging, meaning the economy consumes more than it produces.

## Chapter Summary

This chapter demonstrates the four pillars of linear algebra in economics: the Leontief model establishes the self-consistency of production; the Sraffa model establishes the duality between price and distribution; the minimax theorem governs game equilibrium; and Markov chains characterize the evolution of economic structures.
