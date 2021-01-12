# Expansión Multipolar

Código que calcula simbólicamente los términos de la expansión multipolar del potencial y el campo eléctrico hasta el orden $n=10$. 

## Teoría
En el límite $|\vec{r}|\gg |\vec{r}^{\,\prime}|$, un campo electrostático general se puede escribir como:
$$\phi(\vec{r}) = \frac{1}{4\pi \varepsilon_0} \sum_{n=0}^\infty \frac{(-1)^n}{n!} \,Q_{i_{1} \cdots i_{n}}\,\partial_{i_{1}} \cdots \partial_{i_{n}}\frac{1}{r},
$$
donde se ha definido el momento multipolar de orden $n$ como el tensor de rango $n$:
$$
Q_{i_1 \cdots i_n} := \int_V \rho(\vec{r})\,x_{i_1} \dots x_{i_n}\,dV.
$$

Así, el potencial electrostático se ha descompuesto en una suma de términos de distinto orden en la expansión multipolar:
$$
\phi(\vec{r}) = \sum_{n=0}^\infty \phi^{(n)}(\vec{r}),
$$

donde $\phi^{(n)}(\vec{r})$ es la contribución multipolar de orden $n$, definida por
$$
\phi^{(n)}(\vec{r}) := \frac{1}{4\pi\varepsilon_0}\,\frac{(-1)^n}{n!}\,Q_{i_{1} \cdots i_{n}}\,\partial_{i_{1}} \cdots \partial_{i_{n}}\frac{1}{r}.
$$

Como $\vec{E}(\vec{r}) = -\vec{\nabla} \phi(\vec{r})$, entonces la expansión del campo eléctrico es:
$$
E_{i}(\vec{r}) = \sum_{n=0}^\infty E_{i}^{(n)}(\vec{r}),
$$

donde  $E_{i}^{n}(\vec{r}) = -\partial_{i} \phi^{n}(\vec{r})$ es el $n$-ésimo término de la expansión.

##Requisitos

* Python 3 o superior.
* IPython.
* Sympy.
* Numpy.
* Matplotlib.
* Ipywidgets.

##Estructura de código
A modo de ejemplo, se comparte una serie de jupyter-notebooks donde se ha aplicado el código para determinar el potencial y campo eléctrico de un paralelepípedo de lados $a$, $b$ y $c$ centrado con relación al origen del sistema de coordenadas, el cual se encuentra uniformemente cargado con una densidad volumétrica de carga conocida.

Las clases y funciones importantes se definen en el [módulo principal](https://github.com/nataly-nicole/MultipolarExpansion/blob/main/multipolar_expansion.py). Los jupyter-notebooks son:
* [Expansion_Multipolar - Calcula y Guarda.ipynb](https://github.com/nataly-nicole/MultipolarExpansion/blob/main/Expansion_Multipolar%20-%20Calcula%20y%20Guarda.ipynb).
* [Expansion_Multipolar - Lee y Despliega.ipynb](https://github.com/nataly-nicole/MultipolarExpansion/blob/main/Expansion_Multipolar%20-%20Lee%20y%20Despliega.ipynb).
* [Expansion_Multipolar - Lee y Gráfica.ipynb](https://github.com/nataly-nicole/MultipolarExpansion/blob/main/Expansion_Multipolar%20-%20Lee%20y%20Gráfica.ipynb).
