## Билинейные базисные функции

$$
X_1(x) = \frac{x_{p+1} - x}{h_x},\ X_2(x) = \frac{x - x_p}{h_x}
$$

$$
Y_1(y) = \frac{y_{s+1} - y}{h_y},\ Y_2(x) = \frac{y - y_s}{h_y}
$$

$$
    \psi_i = X_{\mu(i)}Y_{\nu(i)}
$$

---

$$
\begin{align}
\begin{split}

&\int_{\hat{\Omega}}\psi_{k}\nabla\psi_{j}\nabla\psi_{i} d\Omega = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}Y_{\nu(k)}
    (
        \frac{\partial\psi_{j}}{\partial x}\frac{\partial\psi_{i}}{\partial x} + 
        \frac{\partial\psi_{j}}{\partial y}\frac{\partial\psi_{i}}{\partial y}
    )
    dxdy = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}Y_{\nu(k)}
    (
        \frac{\partial X_{\mu(j)}}{\partial x}Y_{\nu(j)}
        \frac{\partial X_{\mu(i)}}{\partial x}Y_{\nu(i)} + 
        \frac{\partial Y_{\nu(j)}}{\partial y}X_{\mu(j)}
        \frac{\partial Y_{\nu(i)}}{\partial y}X_{\mu(i)}
    )
    dxdy = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}Y_{\nu(k)}
    \frac{\partial X_{\mu(j)}}{\partial x}Y_{\nu(j)}
    \frac{\partial X_{\mu(i)}}{\partial x}Y_{\nu(i)}
    dxdy + \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}Y_{\nu(k)}
    \frac{\partial Y_{\nu(j)}}{\partial y}X_{\mu(j)}
    \frac{\partial Y_{\nu(i)}}{\partial y}X_{\mu(i)}
    dxdy = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}
    \frac{\partial X_{\mu(j)}}{\partial x}
    \frac{\partial X_{\mu(i)}}{\partial x}
    dx
\int_{\hat{\Omega}}
    Y_{\nu(k)}
    Y_{\nu(j)}
    Y_{\nu(i)}
    dy + \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}
    X_{\mu(j)}
    X_{\mu(i)}
    dx
\int_{\hat{\Omega}}
    Y_{\nu(k)}
    \frac{\partial Y_{\nu(j)}}{\partial y}
    \frac{\partial Y_{\nu(i)}}{\partial y}
    dy

\end{split}
\end{align}
$$

$$
\begin{align}
\begin{split}

&\int_{\hat{\Omega}}\psi_{k}\psi_{j}\psi_{i} d\Omega = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}Y_{\nu(k)}
    X_{\mu(j)}Y_{\nu(j)}
    X_{\mu(i)}Y_{\nu(i)}
    dxdy = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}
    X_{\mu(j)}
    X_{\mu(i)}
    dx
\int_{\hat{\Omega}}
    Y_{\nu(k)}
    Y_{\nu(j)}
    Y_{\nu(i)}
    dy

\end{split}
\end{align}
$$


$$
\begin{align}
\begin{split}

&\int_{\hat{\Omega}}\psi_{k}\psi_{i} d\Omega = \\
&\int_{\hat{\Omega}}
    X_{\mu(k)}Y_{\nu(k)}
    X_{\mu(i)}Y_{\nu(i)}
    dxdy = \\

&\int_{\hat{\Omega}}
    X_{\mu(k)}
    X_{\mu(i)}
    dx
\int_{\hat{\Omega}}
    Y_{\nu(k)}
    Y_{\nu(i)}
    dy

\end{split}
\end{align}
$$

$$
\begin{align}
\begin{split}

&\int_{\hat{S_x}}\psi_{k}\psi_{j}\psi_{i} dS_x = \\
&\int_{\hat{S_x}}
    X_{\mu(k)}Y_{\nu(k)}
    X_{\mu(j)}Y_{\nu(j)}
    X_{\mu(i)}Y_{\nu(i)}
    dx = \\
&Y_{\nu(k)}Y_{\nu(j)}Y_{\nu(i)}
    \int_{\hat{S_x}}
        X_{\mu(k)}
        X_{\mu(j)}
        X_{\mu(i)}
        dx

\end{split}
\end{align}
$$

$$
\begin{align}
\begin{split}

&\int_{\hat{S_x}}\psi_{k}\psi_{i} dS_x = \\
&\int_{\hat{S_x}}
    X_{\mu(k)}Y_{\nu(k)}
    X_{\mu(i)}Y_{\nu(i)}
    dx = \\
&Y_{\nu(k)}Y_{\nu(i)}\int_{\hat{S_x}}
    X_{\mu(k)}
    X_{\mu(i)}
    dx

\end{split}
\end{align}
$$

$$
\begin{align}
\begin{split}

G^{nx}_{abc} = \int_{\hat{\Omega}}
    X_{a}
    \frac{\partial X_{b}}{\partial x}
    \frac{\partial X_{c}}{\partial x}
    dx = \frac{1}{h_x}
    \begin{pmatrix}
        \frac{1}{2} & -\frac{1}{2} \\
        -\frac{1}{2} & \frac{1}{2}
    \end{pmatrix},
    \begin{pmatrix}
        \frac{1}{2} & -\frac{1}{2} \\
        -\frac{1}{2} & \frac{1}{2}
    \end{pmatrix};\ a,b,c = 1,2

\end{split}
\end{align}
$$

$$
\begin{align}
\begin{split}

M^{nx}_{abc} = \int_{\hat{\Omega}}
    X_{a}
    X_{b}
    X_{c}
    dx = h_x
    \begin{pmatrix}
        \frac{1}{4} & \frac{1}{12} \\
        \frac{1}{12} & \frac{1}{12}
    \end{pmatrix},
    \begin{pmatrix}
        \frac{1}{12} & \frac{1}{12} \\
        \frac{1}{12} & \frac{1}{4}
    \end{pmatrix};\ a,b,c = 1,2

\end{split}
\end{align}
$$

$$
\begin{align}
\begin{split}

M^{x}_{ab} = \int_{\hat{\Omega}}
    X_{a}
    X_{b}
    dx = h_x
    \begin{pmatrix}
        \frac{1}{3} & \frac{1}{6} \\
        \frac{1}{6} & \frac{1}{3}
    \end{pmatrix};\ a,b = 1,2

\end{split}
\end{align}
$$

---

$$
\begin{align}
\hat{G}_{ij} = \sum_{k=1}^{4}\hat{\lambda}_{k}(
    G^{nx}_{\mu(k)\mu(j)\mu(i)}
    M^{ny}_{\nu(k)\nu(j)\nu(i)} +
    M^{nx}_{\mu(k)\mu(j)\mu(i)}
    G^{ny}_{\nu(k)\nu(j)\nu(i)}
)
\end{align}
$$

$$
\begin{align}
\hat{M}^{\gamma}_{ij} = \sum_{k=1}^{4}\hat{\gamma}_{k}(
    M^{nx}_{\mu(k)\mu(j)\mu(i)}
    M^{ny}_{\nu(k)\nu(j)\nu(i)}
)
\end{align}
$$

$$
\begin{align}
\hat{M}^{S^3_x}_{ij} = \sum_{k=1}^{2}\hat{\beta}_{k}(
    Y_{\nu(k)}Y_{\nu(j)}Y_{\nu(i)}
    M^{nx}_{\mu(k)\mu(j)\mu(i)}
)
\end{align}
$$

$$
\begin{align}
\hat{b}^{\Omega}_{i} = \sum_{k=1}^{4}\hat{f}_{k}(
    M^{x}_{\mu(k)\mu(i)}
    M^{y}_{\nu(k)\nu(i)}
)
\end{align}
$$

$$
\begin{align}
\hat{b}^{S^2_x}_{i} = \sum_{k=1}^{2}\hat{\theta}_{k}(
    Y_{\nu(k)}Y_{\nu(i)}
    M^{x}_{\mu(k)\mu(i)}
)
\end{align}
$$

$$
\begin{align}
\hat{b}^{S^3_x}_{i} = 
    \sum_{k=1}^{2}\hat{\beta}_{k}
    \sum_{j=1}^{2}\hat{u}_{\beta,j}(
        Y_{\nu(k)}Y_{\nu(j)}Y_{\nu(i)}
        M^{nx}_{\mu(k)\mu(j)\mu(i)}
    )
\end{align}
$$
