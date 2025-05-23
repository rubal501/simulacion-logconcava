#import "@preview/polylux:0.4.0": *
#import "@preview/physica:0.9.5": *
#import "@preview/metropolis-polylux:0.1.0" as metropolis
#import metropolis: new-section, focus
#import "@preview/lovelace:0.3.0": *
#import "@preview/cades:0.3.0": qr-code



#show: metropolis.setup.with(
  text-font: "New Computer Modern",
  math-font: "New Computer Modern Math",
  code-font: "Fira Code",
  text-size: 23pt,
  footer: [Simulacion estocastica], // defaults to none
)

#let dt = $dd(t)$
#set text(lang:"es")
#set math.equation(numbering: "1.")

#let sigma2 = $sigma^2$
#let Mplus = $M^(+)$
#let Mminus = $M^(-)$

#let uniform = $op("Uniform")$
#let ll = $cal(l)$


#slide[
  #set page(header: none, footer: none, margin: 3em)

 
  #text(size: 1.3em)[
    *Métodos de muestreo para densidades log cóncavas*
  ]


  #metropolis.divider
  
  #set text(size: .8em, weight: "light")
  Javier Roberto Rubalcava Cortes

  Jan 16, 2025

]

#slide[
  = Agenda

  #metropolis.outline
]

#new-section[Simulación usando *aceptación rechazo*]

#slide[
- El metodo de aceptación rechazo nos permite simular muestras de una densidad $f$ usando 
  una densidad envolvente $g$ que satisface en todo el domino $f(x) <= K g(x)$. Se generan 
  muestras candidatas usando $g$ y se aceptan si:
  $
  U <= f(x')/g(x') quad "donde " U ~ uniform([0,1]), x' tilde g/K
  $

// TODO: Poner la figura donde se ve la imagen

#pagebreak()
- La eficiencia de AR depende de la elección de $g$, una buena $g$ cumple:
  - Queremos que $g/K$ sea fácil de muestrear
  - Queremos que $g$ pueda ser evaluada de forma eficiente
  - Queremos que la probabilidad de aceptacion sea lo mas alta posible

- *NO es fácil obtener una densidad $g$ que cumpla todos los puntos de arriba*
]

#new-section[Desigualdades para densidades *LC*]

#slide[
  = Desigualdades para la densidad exacta

  - Suponemos que la densidad $f$ es log-cóncava con dominio $D = RR$.
  - La primera desigualdad que se expone en @devroyeInequalitiesSimulationMethods solo requiere de la 
    localización de la moda $m$ :
  $
  f(x) <= g_1(x) := M min(1,exp(1- abs(x-m) M)) quad "donde" f(m) = M
  $<eq:desigualdad1>

  #pagebreak()
  - Si ademas tenemos acceso a la varianza y la media $sigma2, mu$ entonces 
    se propone la siguiente desigualdad:
  $
  f(x) <= g_2(x) := cases(
    1/sigma quad "si " abs(x-mu) <= (1 + sqrt(3)) sigma,
    1/(abs(x-mu) - sigma sqrt(3)) quad "si " (sqrt(3)+ sqrt(12)) sigma >= abs(x-mu) >= (1 + sqrt(3)) sigma,
    1/(sigma sqrt(12)) exp(3/2 - (abs(x-mu))/(sigma sqrt(12)) ) quad "si " abs(x-mu) >= (sqrt(3)+ sqrt(12)) sigma,
  )
  $ 

  #pagebreak()
  - Si solo tenemos acceso a la media $mu$ podemos tomar la siguiente desigualdad:
  $
  f(x) <= g_3(x) := cases(
    f(mu) e sqrt(3) quad "si " abs(x-mu) <= (1 + 1/(e sqrt(3)))/(f(mu)),
    1/(abs(x-mu) - 1/(f(mu))) quad "si " 2/(f(mu)) >= abs(x-mu) >=  (1 + 1/(e sqrt(3)))/(f(mu)),
    f(mu) exp(2 - (abs(x-mu))/f(mu) ) quad "si " abs(x-mu) >= 2/(f(mu))
  )
  $
#pagebreak()
  - Las ultimas dos desigualdades son resultado de la primera y de las siguientes desigualdades para 
    densidades *LC*
    $
    1/(e sqrt(3)) <= sigma f(mu) <= 1 quad 
    1/sqrt(12) <= sigma f(m) <= 1
    $

#pagebreak()
  - Las areas de las funciones envoventes $g_1 , g_2, g_3$ son:
  $
    integral_RR g_1(t) dt &= 4 \
    integral_RR g_3(t) dt &= 2(1 + sqrt(3)) + 2(log( sqrt(12))) + sqrt(12)/sqrt(3) approx 9.94 \
    integral_RR g_3(t) dt &= 6 + 2 e sqrt(3) + 2 log(sqrt(3)) approx 16 \
  $
  - Tienen area constante sin importar como sea $f$
#pagebreak()
  #figure(
  image("figures/cotas1.svg", width: 50%),
  caption: [
    Funciones envolventes definidas para para la densidad exacta $f$.
  ],
)

]

#slide[
  = Desigualdades para funciones proporcionales
  - Tomamos una función $h_1 prop f$ definida como $h_1(x) = (f(x))/f(m)$ de modo que $h_1 (m) = 1$
  - Si se conoce una constante $Mminus in RR$ que cumple
    $f(m)=M > Mminus$ entonces usado la desigualdad de la @eq:desigualdad1 se obtiene la 
    siguiente función envolvente para $h_1$
    $
    h_1 <= g_4(x) = min(1, exp(1- abs(x-m)Mminus ) )
    $
    El área bajo la función envolvente $g_4$ es proporcional a $M/Mminus$.
#pagebreak()
  - Si solo se tiene acceso $sigma2$ entonces se propone una envolvente la siguiente para $h_1$:
  $
    h_1 <= g_6 := min(1, exp( 1 - ( abs(x-m)/(sigma sqrt(12)))))
    $
    El área bajo esta envolvente es de $4 sigma sqrt(12)$.

#pagebreak()
#figure(
  image("figures/cotas2.svg", width: 50%),
  caption: [
    Funciones envolventes definidas para la función proporcional $h_1 prop f$.
  ],
)
]
#slide[
- Si se calibra una función $h_2 prop f$ como $h_2 (x) = f(x)/f(mu)$ y se tiene acceso a la $mu, sigma2$ de $f$ entonces usando la
  la siguiente cota inferior: $abs(x-m)>= abs(z-mu) - sigma sqrt(3)$, y 
  el hecho de que  $h_2(m) <= h_2(mu) e sqrt(3)$ entonces
  se propone la siguiente envolvente para $h_2$: 
$
h_2 <= g_5 (x) := e sqrt(3) min(1,exp(3/2 - abs(x-mu)/(sigma sqrt(12))))
$
#pagebreak()
#figure(
  image("figures/cotas3.svg", width: 50%),
  caption: [
    Funciones envolventes definidas para la función proporcional $h_2 prop f$.
  ],
)
]

#new-section[Algoritmos de caja negra]

#slide[
  - Usando las funciones envolventes expuestas previamente se proponen varios algoritmos basados en AR,
    por ejemplo usando la desigualdad de la @eq:desigualdad1 se propone:

#let ber = $op("Bernoulli")$
#let expo = $op("Exponencial")$
#show figure: set block(breakable: true)
#figure(
  kind: "algorithm",
  supplement: [Algoritmo],
  pseudocode-list(booktabs: true, title: "Algoritmo 1")[
  + $M = f(m)$
  + $"Accept" = "false"$
  + *do*
    + Tomar $B tilde ber(1/2)$
    + Tomar un signo $S in {-1,1}$ de forma aleatoria
    + Tomar $U tilde uniform(0,1)$
    + *if* $B = 1$ *then*
      + Tomar $V tilde uniform(0,1)$
      + $X = m + (S V)/M$
      + $"Accept" = U M <= f(X)$
    + *else*  *then*
      + Tomar $E tilde expo(1)$
      + $X = m + (S (1+E))/M$
      + $"Accept" = U M exp(-E) <= f(X)$
    + *end*

  + *while* accept  *is false*
  + *return* $X$
]
)]

#new-section[Adaptive rejection sampling]

#slide[
== Envolvente construida con tangentes 
  Las funciones envolventes presentadas previamente dependen de tener accesos
  a
  la media o la varianza o la moda de la densidad $f$.
  - En esta
    sección consideramos una función envolvente cuya construcción solo depende
    de la log concavidad de $f$. Denotamos $ll (x) := log(f)$, como $f$ es
    log-cóncava entonces $ll$ es cóncava. 
  - Por la concavidad sabemos que toda
    recta tangente a $ll$ se mantiene encima de la curva. Si tenemos dos puntos
    $x_1 < x_2$ tales que $ll'(x_1) > 0 quad ll'(x_2) < 0 $ entonces podemos
    definir la siguiente función lineal a trozos:
  $
  u (x) &:= cases(
    ll(x_1) + (x- x_1) ll'(x_1) quad "si" x in (-oo, z ),
    ll(x_2) + (x- x_2) ll'(x_2) quad "si" x in ( z, oo ),
  )\
  z &:= (ll(x_(2)) - x_(2) ll'(x_(2))  - ll(x_(1)) + x_(1) ll'(x_(1))  )/(ll'(x_(2))  - ll(x_1)) 
  $<eq:envelope1>
    Y como $ll$ es cóncava se tiene que $ll <=  u $. 
  - Por la concavidad de $ll$
    el segmento de recta secante que conecta a los puntos $x_1, x_2$ estará bajo la curva de $ll$.
    Definimos la siguiente función lineal a trozos:
  $
  l (x) := cases(
    - oo &quad "Si" x in.not [x_1, x_2],
    ll(x_1) + ( ll(x_2) - ll(x_1))/(x_2 - x_1) (x - x_1) &quad "si" x in [x_1, x_2]
  )
  $<eq:squeezing1>
  - La concavidad de $ll$ garantiza que $ll(x) > l(x)$.
  - A partir de estas dos funciones 
    lineales podemos acotar la densidad $f$ de la siguiente forma:
    $
    exp(l) <= f <= exp(u):= g(x)
    $

#pagebreak()
#let normal = $op("Normal")$
#figure(
  image("figures/2point-envelope.svg", width: 50%),
  caption: [
    Ejemplo de función envolvente $u$ y función de squezzing $l$ para una densidad $f prop normal(0,1)$
    con la construcción dada en @eq:envelope1 y @eq:squeezing1 usando $x_1 = -2 , x_2 = 2$.
  ],
)

#pagebreak()
  - El área bajo cada uno de los segmentos de la función envolvente $g$ sera:

    $
    A_1 = integral_(-oo)^(z_1) exp(u(x)) d x = (exp(ll'(x_1) z_1 ))/(ll'(x_1))exp(b_1) \
    A_2 = integral_(z_1)^(oo) exp(u(x)) d x = (exp(ll'(x_1) z_1 ))/(ll'(x_1))exp(b_1) \
    A = A_1 + A_2
    $
#pagebreak()
  - Definimos las distribuciones de los segmentos de $g/A$:
$
G_1 (x) = (integral_(-oo)^(x) g(x) )/(A_1) = 1 - exp(ll'(x_1) (x - z_1)) \
G_2 (x) = (integral_(z_1)^(x) g(x) )/(A_2) = 1 - exp(ll'(x_2) (x - z_1)) \
$
#let Ginv = $G^(-1)$

#pagebreak()
    Invertimos cada una de las funciones:
$
Ginv_1 (u) = z_1 + (log(1-u) )/(ll'(x_1)) \
Ginv_2 (u) = z_1 + (log(1-u) )/(ll'(x_2)) \
$
  - Si queremos tomar una muestra $x tilde g/A$ basta con seleccionar una 
    una de las componentes con probabilidad $p_1 = (A_1)/A , p_2 = 1 - p_1$ y tomar $U tilde uniform(0,1)$
    entonces $Ginv_i (U) = x$.

#pagebreak()
==  Generalización de la construcción
  - Usando la construcción de la función envolvente a partir de las rectas
    tangentes y de las secantes de $ll$ se propone el método de aceptación
    rechazo adaptativo (ARS) en @gilksAdaptiveRejectionSampling1992 el cual
    consiste de lo siguiente

  - Definimos el conjunto de abscisas $T = x_1, x_2 ,... , x_k$ con la
    condición 
    de $h'(x_1) > 0, quad h'(x_k) <0 $ y $x_1 < x_2 < ... < x_k$.
  - A partir de este conjunto 
    se extiende la construcción de la función envolvente $u$ con respecto al conjunto $T$.
#pagebreak()
  - Primero se calculan los puntos de intersección de las rectas tangentes 
    inducidas por cada punto de $T$ : 
$
z_j = (ll(x_(j+1)) - x_(j+1) ll'(x_(j+1))  - ll(x_(j)) + x_(j) ll'(x_(j))  )/(ll'(x_(j+1))  - ll(x_j)) quad "para" j=1,2,...,k-1 \
$
  - Si el dominio de $f$ es $RR$ entonces asignamos $z_0 = -oo quad z_k = oo$,
    si el dominio esta acotado entonces 
    $z_0 = min(D) quad z_k = max(D)$
    A partir de estos puntos se extiende la construccion de $u$ en la @eq:envelope1 
    para $T$:
$
u_T (x) := cases(
  ll(x_1) + (x- x_1) ll'(x_1) quad "si" x in (z_0, z_1 ),
  ll(x_2) + (x- x_2) ll'(x_2) quad "si" x in ( z_1, z_2 ),
  dots.v,
  ll(x_k) + (x- x_k) ll'(x_k) quad "si" x in ( z_(k-1), z_k ),
)
$
#pagebreak()
  - Extendemos también la construcción de la función de squeezing $l$ para 
    el conjunto $T$ de la siguiente forma
$
l_T (x) = 
cases(
  h(x_1) + m_j (x - x_1) quad "si" x in [x_1 , x_(2)] ,
  h(x_2) + m_j (x - x_2) quad "si" x in [x_2 , x_(3)] ,
  dots.v,
  h(x_(k-1)) + m_(k-1) (x - x_(k-1)) quad "si" x in [x_(k-1), x_k],
  - oo quad "si " x in.not [x_1, x_k]
)\
m_j = ( h(x_(j+1)) - h(x_j))/(x_(j+1) - x_j) quad "para" j = 1,2,..., k-1
$

#pagebreak()
  - Para calcular el área debajo de cada una de las componentes de $g_T =
    exp(u_T)$ usamos la misma 
    formula que con el caso de dos abscisas cuando $ll'(x) != 0$ y agregamos el caso de cuando $ll'(x) = 0$ #footnote[Este 
    caso solo se da cuando $x = m$, en la construcción de dos puntos se supuso
    que ninguna de las abscisas era la ubicación de la moda]
    de la siguiente forma
#let dy = $dd(y)$
$
A_j = integral_(z_j)^(z_(j+1)) g_T(y) dy = cases(
  (z_(j+1) - z_j) exp(B_j) quad "si" ll'(x_j) = 0 ,
  (exp(ll'(x_j) z_(j+1)) - exp(ll(x_j) z_(j+1)) )/(ll'(x_j) * exp(B_j)) quad "si" ll'(x_j) != 0
)\
B_j = ll(x_j) - ll'(x_j) x_j \
A = sum_j A_j
$
    los pesos de cada componente de $g_T$ estarán dados por:
$
p_j = (integral_(z_j)^(z_(j+1)) g_T(y) dy)/A
$
#pagebreak()
    La CDF  de cada componente de $g_T$ esta dada por 
$
G_j(x) = (integral_(z_j)^(x) g_T(y) dy  )/(integral_(z_j)^(z_(j+1)) g_T(y) dy) =
cases(
  (exp(x ll'(x_j)) - exp( z_(j) ll'(x_j) ) )/exp(z_(j+1) ll'(x_j)) - exp( z_(j) ll'(x_j) ) quad "si " ll'(x_j) != 0,
  (x- z_j)/(z_(j+1) - z(j)) quad "si" ll'(z_j) = 0,
)
$
    Invirtiéndola tenemos:
$
Ginv_j (u) = cases(
  1/(ll'(x_j)) log( u exp( ll'(x_j) z_(j+1)) + (1-u) exp( ll'(x_j) z_j)  ) quad "si" ll'(z_j) != 0,
  z_j + u (z_(j+z) - z_j) quad "si" ll'(z_j) = 0,
)
$

#pagebreak()
#figure(
  image("figures/ars-envelope.svg", width: 60%),
  caption: [
Ejemplo de función envolvente $u_T$ y función de squezzing $l_T$ para una densidad $f prop normal(0,1)$
    usando $T= (-1, 0.1, 1)$.
  ],
)

#pagebreak()
    Definimos el siguiente algoritmo para muestrear de $(g_T)/A$

#pseudocode-list(booktabs: true, title:[Muestreo de $g_T$ ])[
  + Tomar $j$ tal que $PP(j = i) = p_i$
  + Tomar $U tilde uniform(0,1)$
  + *if* $ll'(x_j) == 0$
    + *return* $1/(ll'(x_j)) log( U exp( ll'(x_j) z_(j+1)) + (1-U) exp( ll'(x_j) z_j)  ) $
  + *else*
    + *return* $z_j + U (z_(j+1) - z(j)) $
  + *end*
]

    - Cuando tomamos un candidato $x' ~ g/A$ se realiza una prueba de squeezing entre $s_T$ y $g_T$
      antes de hacer la prueba de rechazo entre $f$ y $g_T$.

#pseudocode-list(booktabs: true, name:"Algoritmo ARS")[
  + *function* ARS($N,T_K$)
  + Definir $g_T, s_T$
  + Incializar $X = {}$
  + *while* $|X| < N$
    + Tomar $x tilde g/A$
    + $U tilde uniform(0,1)$
    + *if*  $U <= (s_T (x))/(g_T (x))$ *then*
      + Agregar $x$ a $X$
    + *else then*
      + *if*  $U <= (s_T (x))/(g_T (x))$ *then*
      +  Agregar $x$ a $X$
      + *end* 
      + Actualizar el conjunto de abscisas $T$ añadiendo $x$
      + Re definir $g_T, s_T$
    + *end*
  + *end*
  + *return* $X$
]


]

#new-section[Implementaciones]

#slide[
= Benchmarking

- Se realizo una implementación de cada uno de los algoritmos descritos en
  @devroyeInequalitiesSimulationMethods y de ARS en Julia.
- Para probar los algoritmos de caja negra utilizamos una distribución normal $f$
  con parámetros $mu = 2, sigma2 = 3$.
- Usando las funciones de benchmarking de Julia
  calculamos el costo en tiempo y memoria (alocaciones de memoria) de generar una sola muestra de 
  $f$. Para comparar se uso el 
  método `rand()`
  del paquete de *Distributions.jl*@dahualinJuliaStatsDistributionsjlV0251192025.
  Se registro el tiempo promedio y uso de memoria.

- En la @table:speed podemos ver que todos los algoritmos de caja negra requieren mucho mas tiempo que 
  el algoritmo base, esta diferencia es menor en los algoritmos *1* y *2*, pero aun así no deja de ser 
  significativamente mas lento que el algoritmo original de muestreo.

#figure(
  table(
  columns: 3,
  [], [*Alocaciones*], [*Tiempo promedio*],
  [Algoritmo base], [0 ], [4.401 ns ±  0.518 ns ],
  [Algoritmo 1], [21 ], [543.539 ns ± 957.936 ns],
  [Algoritmo 2], [33 ], [1.155 μs ±  1.253 μs ],
  [Algoritmo 3], [45 ], [2.217 μs ±   5.941 μs ],
  [Algoritmo 4], [28 ], [768.602 ns ± 654.707 ns ],
  [Algoritmo 5], [45 ], [1.098 μs ± 812.810 ns ],
  [Algoritmo 6], [39 ], [1.077 μs ±  1.162 μs ],
  )
)<table:speed>

- Las diferencias tan marcadas entre el desempeño de los algoritmos implementados en este proyecto contra el 
  algoritmo base usado en Julia también es importante considerar que el algoritmo de muestreo que se utiliza es 
  *Xorshiro256* @marsagliaXorshiftRNGs2003 el cual destaca por su eficiencia tanto en velocidad como en memoria.


- Para probar la implementación del algoritmo ARS se considero la misma densidad $f$, solo que 
  en este caso para la prueba se generaron $N = 10000$ muestras, inicializando el conjunto 
  de abscisas como $T_k = {-1, 0.5, 1}$. Los resultados en la tabla @table:speed2
  muestran que el tiempo y las alocaciones son altas para (ARS), al igual que con el Algoritmo 1

#figure(
  table(
  columns: 3,
  [], [*Alocaciones*], [*Tiempo promedio*],
  [Algoritmo base], [0 ], [3.245 μs ±  19.783 μs  ],
  [ARS], [21], [9.973 ms ±  2.190 ms  ],
  [Algoritmo 1], [21 ], [5.775 ms ±  1.117 ms  ],
  )
)<table:speed2>

= Comparación del numero de rechazos en ARS
- Otra comparación que se realizo entre el algoritmo 1 y (ARS) fue el numero de muestras que se rechazan 
  cuando se simula una muestra de tamaño $N=10000$.
- En la @table:rejections se puede ver que 
  (ARS) en promedio rechaza muy pocas muestras candidatas. 
- Esto se debe a que el conjunto final de abscisas $T$
  generado por (ARS) tiene (en promedio) un tamaño grande, lo cual resulta en una función envolvente $g_T$ extremadamente 
  ajustada a $f$, lo cual parece ser también la causa de que requiera tiempo extra, ya que al agregar múltiples puntos a 
  $T$ durante el _loop_ de muestreo de (ARS) entonces aumenta la complejidad de tomar muestras de $(g_T)/A$.


#figure(
  table(
  columns: 3,
  [], [*\# Rechazos (promedio)*], [*\# Abscisas (promedio)*],
  [Algoritmo base], [NA], [NA],
  [ARS], [5], [100],
  [Algoritmo 1], [5964], [NA],
  )
)<table:rejections>

]

#new-section[Puntos finales]
#slide[
- ARS *NO* es optimo ya que en la mayoría de los casos es posible diseñar un algoritmo especializado.
- La ventaja principal de ARS es su universalidad. Solo requerimos la forma cerrada de $f$ y su derivada 
  para poder simular   muestras 
- En un trabajo posterior Gilks y Walters propusieron una versión de ARS que no necesita la derivada de $f$  
  @gilks1992derivative. También se ha propuesto una generalización para densidades que no 
  son log cóncavas @martino2011generalization.
- También se ha diseñado un método para optimizar
  la construcción del conjunto de abscisas inicial para 
  aumentar la eficiencia del muestreo @james2024automated.
]
#focus[
Muchas gracias!!

#figure(
qr-code("https://typst.app", width: 10cm),
    caption:[Repositorio con el código del proyecto.],
    supplement: none
  )


]

#slide[

  #bibliography("references.bib")
]
