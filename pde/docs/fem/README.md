# Метод конечных элементов

## Эллиптические уравнения

Уравнение эллиптического типа можно записать с следующей форме:

$$
\nabla(\lambda(u,x)\nabla u(x)) + \gamma(u,x)u(x) = f(u,x)
$$

Для задачи Коши, даются следующие краевые условия:

$$
u|_{S_1} = u_g(x)
$$

$$
\left.\lambda(u,x)\frac{\partial u}{\partial n}\right\vert_{S_2} = \theta(u,x)
$$

$$
\left.\lambda(u,x)\frac{\partial u}{\partial n}\right\vert_{S_3}
    + \beta(u,x)(u|_{S_3}
    - u_{\beta}(u,x)) = 0
$$

### Вариационная постановка

Эквивалентная вариационная постановка в форме уравнения Галёркина:

$$
\begin{align*}
\begin{split}
&\int_{\Omega}\lambda(u,x)\nabla{u}\nabla{v_0}d\Omega
    + \int_{\Omega}\gamma(u,x)uv_0d\Omega
    + \int_{S_3}\beta(u,x)uv_0dS \\
&= \int_{\Omega}f(u,x)v_0d\Omega
        + \int_{S_2}\theta(u,x)v_0dS
        + \int_{S_3}\beta(u,x)u_{\beta}(u,x)v_0dS, \forall v_0 \in H_0^1
\end{split}
\end{align*}
$$

Раскладывая функции $u$ и $v_0$ по базису, переходим к конечноэлементной СЛАУ:

$$
\begin{align*}
\begin{split}
&\sum_{j=1}^{n}(\int_{\Omega}\lambda(u,x)\nabla{\psi_j}\nabla{\psi_i}d\Omega
    + \int_{\Omega}\gamma(u,x)\psi_{j}\psi_{i}\Omega
    + \int_{S_3}\beta(u,x)\psi_{j}\psi_{i}dS)q_j \\
&= \int_{\Omega}f(u,x)\psi_{i}d\Omega
        + \int_{S_2}\theta(u,x)\psi_{i}dS
        + \int_{S_3}\beta(u,x)u_{\beta}(u,x)\psi_{i}dS
\end{split}
\end{align*}
$$

или, в более удобной форме записи:

$$
A_{ij} = G_{ij} + M^{\gamma}_{ij} + M^{S_3}_{ij}
$$

$$
b_{i} = b^{\Omega}_{i} + b^{S_2}_{i} + b^{S_3}_{i}
$$

Параметры нужно расложить по базисным функциям на конечном элементе. Процесс типовой. Для $\lambda$:

$$
\lambda(u,x) = \sum_{k=1}^{m}\hat{\lambda}_k(u)\hat{\psi}_k
$$

Выражения для компонент локальных матриц и векторов зависят от размерности пространства задачи и выбранного базиса. В случае, если параметры зависят от решения $u$, задача считается [нелинейной](non/README.md).
