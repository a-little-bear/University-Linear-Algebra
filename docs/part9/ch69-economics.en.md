# Chapter 69: Applications of Linear Algebra in Economics

<div class="context-flow" markdown>

**Prerequisites**: Linear Equations (Ch1) · Matrix Operations (Ch2) · Eigenvalues (Ch6) · Non-negative Matrices (Ch17) · M-Matrices (Ch38)

**Chapter Outline**: Leontief Input-Output Model $	o$ Hawkins-Simon Conditions $	o$ Sraffa Model $	o$ Input-Output Structure $	o$ Linear Exchange Model $	o$ Game Theory $	o$ Von Neumann Minimax Theorem $	o$ Linear Programming Duality

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
    For any non-negative demand $d \ge 0$, a non-negative output $x \ge 0$ exists if and only if all principal minors of $I-C$ are positive.

---

## 69.2 Game Theory and Duality

!!! theorem "Theorem 69.7 (Von Neumann Minimax Theorem)"
    A two-person zero-sum game always possesses a mixed-strategy equilibrium. Solving for this equilibrium is equivalent to solving a pair of dual linear programming problems.

---

## Exercises

1. **[Fundamentals] Why must the spectral radius of the consumption matrix $C$ in the Leontief model be less than 1?**
   ??? success "Solution"
       Only if $ho(C) < 1$ can the series $(I-C)^{-1} = I + C + C^2 + \dots$ converge and yield a non-negative result. Economically, this represents a "productive" system where internal consumption is less than total output.

2. **[Calculation] Given $C = \begin{pmatrix} 0.5 & 0.4 \ 0.2 & 0.5 \end{pmatrix}$. Verify the Hawkins-Simon conditions.**
   ??? success "Solution"
       $I-C = \begin{pmatrix} 0.5 & -0.4 \ -0.2 & 0.5 \end{pmatrix}$.
       Principal minors: $0.5 > 0$; $\det(I-C) = 0.25 - 0.08 = 0.17 > 0$. The conditions are satisfied.

3. **[Economic Intuition] Explain the meaning of the elements in the $j$-th column of the Leontief inverse matrix $(I-C)^{-1}$.**
   ??? success "Solution"
       The column describes the required increase in total output across all sectors of the economy when the final demand for the $j$-th sector increases by 1 unit.

4. **[Game Theory] If the payoff matrix $A$ is skew-symmetric, what is the value of the game?**
   ??? success "Solution"
       The value of the game must be 0. Since the positions of both players are perfectly symmetric, the expected payoffs under mixed-strategy equilibrium cancel each other out.

5. **[Shadow Prices] In linear programming, what does a positive shadow price for a particular resource imply?**
   ??? success "Solution"
       It implies that the resource is a bottleneck for production. Increasing the supply of that resource by one unit will increase the total profit by the value of the shadow price.

6. **[Income Mobility] What does a smaller second eigenvalue $\lambda_2$ of an income mobility matrix represent?**
   ??? success "Solution"
       It represents stronger social mobility, meaning that initial wealth disparities are reshuffled more quickly, and the system converges to its steady-state distribution faster.

7. **[Sraffa] In the Sraffa price model, what determines the upper bound of the profit rate $r$?**
   ??? success "Solution"
       It is determined by the spectral radius $ho(A)$ of the input-output matrix $A$. The maximum profit rate is $R = 1/ho(A) - 1$.

8. **[Calculation] Solve for the equilibrium of the payoff matrix $A = \begin{pmatrix} 1 & -1 \ -1 & 1 \end{pmatrix}$.**
   ??? success "Solution"
       The mixed strategy equilibrium is $x^* = [0.5, 0.5]^T, y^* = [0.5, 0.5]^T$, and the value of the game is 0.

9. **[Perron] Why is the equilibrium output ratio of a closed economic system unique?**
   ??? success "Solution"
       Because the matrix of a closed system is an irreducible stochastic matrix. According to the Perron-Frobenius theorem, the positive eigenvector corresponding to the eigenvalue 1 is unique up to scaling.

10. **[Structural Analysis] How is an economic crisis expressed in the language of matrices?**
    ??? success "Solution"
        An economic crisis corresponds to $ho(C) \ge 1$. This means that internal consumption exceeds production limits, causing the Leontief inverse matrix to produce negative values or diverge, leading to a collapse of the algebraic structure.

## Chapter Summary

This chapter demonstrates the four pillars of linear algebra in economics: the Leontief model establishes the self-consistency of production; the Sraffa model establishes the duality between price and distribution; the minimax theorem governs game equilibrium; and Markov chains characterize the evolution of economic structures.
