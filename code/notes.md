AML project notes
=================

survival model
$$h(t) = \lambda_0(t) \exp(\beta X)$$

Integrated hazard
$$H(t) = \int_0^t\lambda_0(t) dt\exp(\beta X) = H_0(t) \exp(\beta X)$$

\begin{align}
Pr[T < t] &= 1/h(t) \int_0^t \exp(-h(t) dt)\cr
	&= 1/h(t) \exp(-H(t))\cr
	&= 1/\lambda_0(t) \exp(-H_0(t) \times \exp(\beta X)) / \exp(\beta X)
\end{align}