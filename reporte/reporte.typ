//#import "@preview/charged-ieee:0.1.3": ieee

#import "@preview/rubber-article:0.4.0": *
#import "@preview/lovelace:0.3.0": *
#import "@preview/physica:0.9.5": *


#show: article.with()

#maketitle(
  title: "Métodos de muestreo de densidades log-concavas",
  authors: ("Javier Roberto Rubalcava Cortes",),
  date: "September 2024",
)


#set text(lang:"es")
#set math.equation(numbering: "1.")

#let sigma2 = $sigma^2$
#let Mplus = $M^(+)$
#let Mminus = $M^(-)$

#let uniform = $op("Uniform")$
#let ll = $cal(l)$
= Introducción

El método de aceptación-rechazo (AR) es uno de los procedimientos fundamentales
para el muestreo de variables aleatorias con distribuciones complejas. A
diferencia del método de inversión, que requiere que la función de distribución
acumulativa (CDF) sea analíticamente invertible o al menos factible de invertir
numéricamente—una condición rara vez cumplida en situaciones prácticas—el
método AR no impone esta restricción. En cambio, AR se basa en la construcción
de una función envolvente $g$, que satisface la condición $f(x) <= K g(x)$ para
toda $x$, siendo $f$ la densidad objetivo y $c$ una constante
positiva conocida como constante de aceptación @devroyeNonUniformRandomVariate1986.

La eficiencia del método AR depende directamente de qué tan ajustada sea esta
función envolvente respecto a la densidad objetivo; cuanto más cercana sea la 
probabilidad de aceptacion 
 $1/(K)$ a $1$ más
eficiente será el muestreo.

En este proyecto analizaremos diversas desigualdades útiles para densidades
log-cóncavas y exploraremos cómo dichas desigualdades pueden emplearse para
construir funciones envolventes adecuadas para el método AR. Adicionalmente,
estudiaremos métodos específicos para construir funciones envolventes mediante
combinaciones de componentes exponenciales, derivadas directamente de la
función logarítmica de la densidad, $log(f)$. Finalmente, presentaremos
un algoritmo que utiliza estas funciones envolventes adaptativas para
realizar muestreo eficiente mediante aceptación-rechazo.

== Desigualdades para densidades log-cóncavas

Suponemos que la densidad $f$ es log-cóncava con dominio $D = RR$. La primera 
desigualdad que se expone en @devroyeInequalitiesSimulationMethods solo requiere de la 
localización de la moda $m$, entonces:
$
f(x) <= g_1(x) := M min(1,exp(1- abs(x-m) M)) quad "donde" f(m) = M
$<eq:desigualdad1>
Si ademas tenemos acceso a la desviación estándar y la media $sigma2, mu$ entonces 
se propone la siguiente desigualdad:
$
f(x) <= g_2(x) := cases(
  1/sigma quad "si " abs(x-mu) <= (1 + sqrt(3)) sigma,
  1/(abs(x-mu) - sigma sqrt(3)) quad "si " (sqrt(3)+ sqrt(12)) sigma >= abs(x-mu) >= (1 + sqrt(3)) sigma,
  1/(sigma sqrt(12)) exp(3/2 - (abs(x-mu))/(sigma sqrt(12)) ) quad "si " abs(x-mu) >= (sqrt(3)+ sqrt(12)) sigma,
)
$ 
La función envolvente $g_2$ esta compuesta por por tres segmentos, un segmento uniforme cuando $abs(x-mu) <= (1 + sqrt(3)) sigma$
y dos colas decrecientes. El área bajo $g_2$ es:
$
integral_RR g_2 = 2(1 + sqrt(3)) + 2(log( sqrt(12))) + sqrt(12)/sqrt(3) approx 9.94
$
El peso de cada una de las componentes es de:
$
p_1 = 2(1+sqrt(3))/q quad p_2 = 2 log(sqrt(12))/q quad p_3 / q = 2/q quad q = p_1 + p_2 + p_3
$
Por otro lado si solo se tiene la media $mu$ se puede proponer 
una envolvente similar usando la siguiente desigualdad:
$
f(x) <= g_3(x) := cases(
  f(mu) e sqrt(3) quad "si " abs(x-mu) <= (1 + 1/(e sqrt(3)))/(f(mu)),
  1/(abs(x-mu) - 1/(f(mu))) quad "si " 2/(f(mu)) >= abs(x-mu) >=  (1 + 1/(e sqrt(3)))/(f(mu)),
  f(mu) exp(2 - (abs(x-mu))/f(mu) ) quad "si " abs(x-mu) >= 2/(f(mu))
)
$
el área bajo $g_3$:
$
integral_RR g_3 = 6 + 2 e sqrt(3) + 2 log(sqrt(3))
$
de modo que los pesos de las componentes serán:
$
p_1 = (2 e sqrt(3) + 2)/q quad p_2 = (2 log(sqrt(12)) + 2)/q quad (p_3)/q = 2/q quad q = p_1 + p_2 + p_3
$
Si es costoso evaluar exactamente la densidad $f$ entonces se pueden considerar funciones 
proporcionales $h_1, h_2 prop f$ calibradas de la siguiente forma:
$
h_1(M) = 1 quad h_2(mu) =1
$
Está forma resulta útil por ejemplo en el caso de la distribución log-gamma donde para la densidad 
$f(x) = (exp(a x - e^x))/gamma(a)$, se puede definir una función proporcional:
$
h_1(x) = f(x)/f(m) = ((exp(a x - e^x))/gamma(a))/(1/gamma(a) (a/e)^a ) = 
exp(a (x-m) + a - e^x)
$
Usando la función $h_1$ se evita evaluar la función gamma en cada prueba del paso de rechazo.

Si se tiene una función proporcional $h_1$ y ademas se conoce una constante $Mminus in RR$ tal que 
$f(m)=M > Mminus$ entonces usado la desigualdad de la @eq:desigualdad1 se obtiene la 
siguiente función envolvente para $h_1$
$
h_1 <= g_4(x) = min(1, exp(1- abs(x-m)Mminus ) )
$
El área bajo la función envolvente $g_4$ es proporcional a $M/Mminus$, dependiendo de que tan 
ajustada sea la desigualdad $Mminus < M$ se podrá muestrear usando AR de forma mas eficiente. Si solo se 
tiene acceso $sigma2$ entonces se propone una envolvente la siguiente para $h_1$:
$
h_1 <= g_6 := min(1, exp( 1 - ( abs(x-m)/(sigma sqrt(12)))))
$
El área bajo esta curva es de $4 sigma sqrt(12)$.

Si se calibra una función $h_2 prop f$ y se tiene acceso a la $mu, sigma2$ de $f$ entonces usando la
desigualdad anterior junto a la 
desigualdad de Johnson-Rogers @johnsonMomentProblemUnimodal1951, la 
cual nos da la siguiente cota inferior: $abs(x-m)>= abs(z-mu) - sigma sqrt(3)$, y 
el hecho de que  $h_2(m) <= h_2(mu) e sqrt(3)$ entonces
se propone la siguiente envolvente para $h_2$
$
h_2 <= g_5 (x) := e sqrt(3) min(1,exp(3/2 - abs(x-mu)/(sigma sqrt(12))))
$



#figure(
  image("figures/cotas1.svg", width: 80%),
  caption: [
    Funciones envolventes definidas para para la densidad exacta $f$.
  A la izquierda se muestran las envolventes sobre $log(f)$, mientras que a la 
  derecha se muestran sobre $f$.
  ],
)

#figure(
  image("figures/cotas2.svg", width: 80%),
  caption: [
    Funciones envolventes definidas para la función proporcional $h_1 prop f$.
    A la izquierda se muestran las envolventes sobre $log(h_1)$, mientras que a la 
    derecha se muestran sobre $h_1$.
  ],
)

#figure(
  image("figures/cotas3.svg", width: 80%),
  caption: [
    Función envolvente definidas para $h_2$.
    A la izquierda se muestran las envolventes sobre $log(h_2)$, mientras que a la 
    derecha se muestran sobre $h_2$.
  ],
)
== Metodos blackbox

En @devroyeInequalitiesSimulationMethods se proponen varios algoritmos de muestreo 
a partir de las envolventes presentadas previamente. Para simplificar el proceso de muestreo
se hacen las siguientes observaciones:
- Si la densidad es trasladada de modo que la moda se encuentra en $m=0$ entonces cada una 
  de las funciones envolventes se vuelve simétrica. 
- Cuando se toma una muestra de la distribución dada por la envolvente normalizada, entonces basta
  con tomar la muestra del lado derecho y asignar el signo de forma aleatoria.
- Se puede muestrear de cada una de las componentes de las funciones envolventes usando 
  muestras de distribuciones exponenciales y uniformes, y aplicándoles el escalado apropiado

A partir de esto se proponen varios algoritmos basados en AR usando las 
funciones envolventes expuestas previamente, por ejemplo usando la desigualdad de la @eq:desigualdad1
se propone:

#let ber = $op("Bernoulli")$
#let expo = $op("Exponencial")$
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
)

El resto de funciones envolventes se pueden usar para definir algoritmos usando aceptación rechazo. Dependiendo 
de la función envolvente que se usa usamos la siguiente nomenclatura para los resultados de la siguiente sección:
#figure(
table(
  columns: 6,
  stroke: none,
  [Numero], [*Evaluación $f$?*], [*Función envolvente*],
  [Algoritmo 1], [exacta], [$g_1$],
  [Algoritmo 2], [exacta], [$g_2$],
  [Algoritmo 3], [exacta], [$g_3$],
  [Algoritmo 4], [$h_1$], [$g_4$],
  [Algoritmo 5], [$h_2$], [$g_5$],
  [Algoritmo 6], [$h_1$], [$g_6$],
),
  caption: [Nomenclatura de algoritmos.]
)<tab:nomen1>
== Envolvente construida con tangentes 
Las funciones envolventes presentadas previamente dependen de tener accesos a la media o la desviación estándar o la moda 
de la densidad $f$. En esta sección consideramos una función envolvente cuya construcción solo depende de la log concavidad de $f$.
Denotamos $ll (x) := log(f)$, como $f$ es log-cóncava entonces $ll$ es cóncava. Por la concavidad 
sabemos que toda recta tangente a $ll$ se mantiene encima de la curva. Si tenemos dos puntos $x_1 < x_2$
tales que $ll'(x_1) > 0 quad ll'(x_2) < 0 $ entonces podemos definir la siguiente función lineal a 
trozos:
$
u (x) &:= cases(
  ll(x_1) + (x- x_1) ll'(x_1) quad "si" x in (-oo, z ),
  ll(x_2) + (x- x_2) ll'(x_2) quad "si" x in ( z, oo ),
)\
z &:= (ll(x_(2)) - x_(2) ll'(x_(2))  - ll(x_(1)) + x_(1) ll'(x_(1))  )/(ll'(x_(2))  - ll(x_1)) 
$<eq:envelope1>
Y como $ll$ es cóncava se tiene que $ll <=  u $. Además como $ll$ es cóncava entonces
el segmento de recta secante que conecta a los puntos $x_1, x_2$ estará bajo la curva de $ll$.
Definimos la siguiente función lineal a trozos:
$
l (x) := cases(
  - oo &quad "Si" x in.not [x_1, x_2],
  ll(x_1) + ( ll(x_2) - ll(x_1))/(x_2 - x_1) (x - x_1) &quad "si" x in [x_1, x_2]
)
$<eq:squeezing1>
Al estar compuesta por el segmento de la tangente que conecta $x_1$ con $x_2$ y rectas constantes
iguales a $- oo$, tenemos la garantía de que $ll(x) > l(x)$. A partir de estas dos funciones 
lineales podemos acotar la densidad $f$ de la siguiente forma:
$
exp(l) <= f <= exp(u):= g(x)
$

#let normal = $op("Normal")$
#figure(
  image("figures/2point-envelope.svg", width: 80%),
  caption: [
    Ejemplo de función envolvente $u$ y función de squezzing $l$ para una densidad $f prop normal(0,1)$
    con la construcción dada en @eq:envelope1 y @eq:squeezing1 usando $x_1 = -2 , x_2 = 2$.
  ],
)
El área bajo cada uno de los segmentos de $g$ sera:

$
A_1 = integral_(-oo)^(z_1) exp(u(x)) d x = (exp(ll'(x_1) z_1 ))/(ll'(x_1))exp(b_1) \
A_2 = integral_(z_1)^(oo) exp(u(x)) d x = (exp(ll'(x_1) z_1 ))/(ll'(x_1))exp(b_1) \
A = A_1 + A_2
$
Definimos las CDF de los segmentos de $g/A$ de la siguiente forma:
$
G_1 (x) = (integral_(-oo)^(x) g(x) )/(A_1) = 1 - exp(ll'(x_1) (x - z_1)) \
G_2 (x) = (integral_(z_1)^(x) g(x) )/(A_2) = 1 - exp(ll'(x_2) (x - z_1)) \
$
#let Ginv = $G^(-1)$
Invertimos cada una de las funciones:
$
Ginv_1 (u) = z_1 + (log(1-u) )/(ll'(x_1)) \
Ginv_2 (u) = z_1 + (log(1-u) )/(ll'(x_2)) \
$
De modo que si queremos tomar una muestra $x tilde g/A$ basta con seleccionar una 
una de las componentes con probabilidad $p_1 = (A_1)/A , p_2 = 1 - p_1$ y tomar $U tilde uniform(0,1)$
entonces $Ginv_i (U) = x$.



Usando la construcción de la función envolvente a partir de las rectas tangentes y de las secantes
de $ll$ se propone el método de aceptación rechazo adaptativo (ARS) en @gilksAdaptiveRejectionSampling1992 
el cual consiste de lo siguiente

Se define un conjunto de absisas $T = x_1, x_2 ,... , x_k$ con la condicion 
de que $h'(x_1) > 0, quad h'(x_k) <0 $ y $x_1 < x_2 < ... < x_k$. A partir de este conjunto 
se extiende la construccion de la funcion envolvente $u$ con respecto al conjunto $T$.
Primero se calculan los puntos de intersección de las rectas tangentes 
a cada punto de $T$ : 
$
z_j = (ll(x_(j+1)) - x_(j+1) ll'(x_(j+1))  - ll(x_(j)) + x_(j) ll'(x_(j))  )/(ll'(x_(j+1))  - ll(x_j)) quad "para" j=1,2,...,k-1 \
$
Si el dominio de $f$ es $RR$ entonces asignamos $z_0 = -oo quad z_k = oo$, si el dominio esta acotado entonces 
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
Extendemos también la construcción de la función de squeezing $l$ para 
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

Para calcular el area debajo de cada una de las componentes de $g_T = exp(u_T)$ usamos la misma 
formula que con el caso de dos abscisas cuando $ll'(x) != 0$ y agregamos el caso de cuando $ll'(x) = 0$ #footnote[Este 
caso solo se da cuando $x = m$, en la construcción de dos puntos se supuso que ninguna de las abscisas era la ubicación de la moda]
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

#figure(
  image("figures/ars-envelope.svg", width: 80%),
  caption: [
Ejemplo de función envolvente $u_T$ y función de squezzing $l_T$ para una densidad $f prop normal(0,1)$
    usando $T= (-1, 0.1, 1)$.
  ],
)

De modo que proponemos el siguiente algoritmo para muestrear de $(g_T)/A$

#pseudocode-list(booktabs: true, title:[Muestreo de $g_T$ ])[
  + Tomar $j$ tal que $PP(j = i) = p_i$
  + Tomar $U tilde uniform(0,1)$
  + *if* $ll'(x_j) == 0$
    + *return* $1/(ll'(x_j)) log( U exp( ll'(x_j) z_(j+1)) + (1-U) exp( ll'(x_j) z_j)  ) $
  + *else*
    + *return* $z_j + U (z_(j+1) - z(j)) $
  + *end*
]

Otro elemento en el que difiere ARS de método de rechazo convencional es que 
se considera una prueba de _squeezing_ para esto se utiliza el cociente $l_K /g_K$, cuando 
las muestras son rechazadas por la prueba de _squeezing_ entonces se añaden a $T_k$,
de modo que $T_(k+1) = T_k union {x}$, este nuevo conjunto de abscisas induce a nuevas 
funciones $l_(K+1), u_(K+1)$, entre mas puntos tenga el conjunto de abscisas entonces 
$g_(K)$ estará mas ajustada sobre $f$ lo cual reduce las muestras rechazadas.


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






//
// == Métodos _blackbox_
//
// Para los metodos blackbox. Solo se asume que la densidad $f$ es log concava 
// con moda en $m$
// para construir las distintas funciones envolventes $g$ se considera la siguiente desigualdad
//
// Esta desigualdad se ajusta dependiendo de la información que se tiene de la densidad $f$ y de si 
// se esta usando una función proporcional $h prop x$ en lugar de $f$. 
//
// === Se puede evaluar la función $f$
//
// Cuando se tiene acceso a la densidad $f$ y es posible evaluarla 
// de forma *eficiente*. El primer caso que se considera es cuando
// tambien se tiene 
// la locación de la moda $m$, entonces se utiliza la
// desigualdad de @eq:desigualdad1 para hacer aceptación rechazo usando como función envolvente $g = M min(1,exp(1- abs(x-m) M))$,
// la función se puede expresar como:
// $
// g(x) = cases(
//   M quad "si" abs(x-m) <= 1/M,
// exp(1- abs(x-m) M) quad "si" abs(x-m) > 1/M
// )
// $
// Es decir si $x in [m-1/M, m + 1/M] $ entonces $g prop M$, mientras que 
// si $x not in [m-1/M, m + 1/M]$ se tiene $g prop exp (1- abs(x-m) M)$.
// que para $g(x) = $, de modo que ees una mezcla donde si $x in $ entonces 
// $g(x) prop U $ y si $g(x) prop E $, de modo que para tomar muestras de $X tilde G$
// basta con $B tilde B$ y con ello decidir . Ademas como la funcion $g$ es simetrica,
// se hace el muestreo sobre el semi eje derecho y se le asigna un signo de forma aleatoria.
//
// #let finv = $f^{-1}$
// Si la densidad $f$ es decreciente sobre $[0,oo]$ entonces se toma en cuenta que 
// la region $A = { x }$ puede ser descrita de forma usando la inversa $f$ como
// $A = {y }$ y se toman muestras $(X,Y)$ de $A$ de la siguiente forma: primero se 
// toma $Y tilde g(y) prop finv(x)$ y $X tilde U([0,Y])$ de modo que $(X,Y)$ es un punto 
// del rectagulo de altura $Y$ y se acepta cuando $Y <= f(X) $
//
// Si se conocen $mu, sigma2$ entonces se pueden definir  funciones $Mplus(sigma2), Mminus(sigma2)$ tales que 
// $
// Mminus <= f(m) <= Mplus
// $
// Entonces la desigualdad 1 puede ser ajustada
// usando la desigualdad de Johnson Rogers 
//
//
// Si solo se tiene acceso a $mu$ se toma la desiguald 1 y por la uni modealidad 
// de $f$ etnonces $abs(mu-m)f(mu)<= 1$
// entonces:
//
//
// El area de esta funcion evolvente es de $ 6 + 2 e sqrt(3) + log(3) approx 16$
// Estas dos envolventes son mezclas de tres componentes 
//
// Si no se puede evaluar directamente $f$ pero se tiene aceso a una funcion $h prop f$. 
// entonces $h(x) = f(x)/f(m)$, entonces $h(m) = 1$ y por la desigualdad 1 entonces
// $h(x)<= exp(1- abs(x-m)/(sigma sqrt(12)))$ esta nueva envolvente tiene area 
// $(4 sigma sqrt(12))/integral h = 8 sqrt(3) f(m) sigma <= 8 sqrt(3) approx 13$
//
// Si no se conoce la localizacion de la moda pero si $mu,sigma2$. Si consideramos $h prop f$ calibrada 
// de forma que $h(mu) = 1$. Usando la desigualdad 2 y $h(m) <= h(mu) e sqrt(3)$ se obtiene:
// $
// h(x) <= e sqrt(3) min( 1, exp(3/2 - abs(x-mu)/(sigma sqrt(12)) ))
// $
// El area de esta envolvente es de $30 e f(m) sigma <= 30 e approx 81$
//
// == Metodo adaptativo

= Implementaciones

Se realizo una implementación de cada uno de los algoritmos descritos en @tab:nomen1 y de ARS en Julia.
Para probar los algoritmos de la @tab:nomen1 utilizamos una distribución normal $f$
con parámetros $mu = 2, sigma2 = 3$. Usando las funciones de benchmarking de Julia
calculamos el costo en tiempo y memoria (alocaciones de memoria) de generar una sola muestra de 
$f$. Para comparar se implementaron funciones para tomar una muestra usando el 
método `rand()`
del paquete de *Distributions.jl*@dahualinJuliaStatsDistributionsjlV0251192025,
se registro el tiempo promedio y uso de memoria.
En la @table:speed podemos ver que todos los algoritmos de la @tab:nomen1 requieren mucho mas tiempo que 
el algoritmo base, esta diferencia es menor en los algoritmos *1* y *2*, pero aun así no deja de ser 
significativamente mas lento que el algoritmo original de muestreo

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

Para probar la implementación del algoritmo ARS se considero la misma densidad $f$, solo que 
en este caso para la prueba se generaron $N = 10000$ muestras, inicializando el conjunto 
de abscisas como $T_k = {-1, 0.5, 1}$. Los resultados en la tabla @table:speed2
muestran que el tiempo y las alocaciones son altas para (ARS), al igual que con el Algoritmo 1
Otra comparación que se realizo entre el algoritmo 1 y (ARS) fue el numero de muestras que se rechazan 
cuando se intenta tomar una muestra de tamaño $N=10000$, en la @table:rejections se puede ver que 
(ARS) en promedio rechaza muy pocas muestras candidatas. Esto se debe a que el conjunto final de abscisas $T$
generado por (ARS) tiene en promedio un tamaño grande, lo cual resulta en una función envolvente $g_T$ extremadamente 
ajustada a $f$, lo cual parece ser también la causa de que requiera tiempo extra, ya que al agregar múltiples puntos a 
$T$ durante el _loop_ de muestreo de (ARS) entonces aumenta la complejidad de tomar muestras de $(g_T)/A$.


#figure(
  table(
  columns: 3,
  [], [*Alocaciones*], [*Tiempo promedio*],
  [Algoritmo base], [0 ], [3.245 μs ±  19.783 μs  ],
  [ARS], [21], [9.973 ms ±  2.190 ms  ],
  [Algoritmo 1], [21 ], [5.775 ms ±  1.117 ms  ],
  )
)<table:speed2>


#figure(
  table(
  columns: 3,
  [], [*\# Rechazos (promedio)*], [*\# Abscisas (promedio)*],
  [Algoritmo base], [NA], [NA],
  [ARS], [5], [100],
  [Algoritmo 1], [5964], [NA],
  )
)<table:rejections>

Las diferencias tan marcadas entre el desempeño de los algoritmos implementados en este proyecto contra el 
algoritmo base usado en Julia también es importante considerar que el algoritmo de muestreo que se utiliza es 
*Xorshiro256* el cual destaca por su eficiencia tanto en velocidad como en memoria.


= Discusión

Los algoritmos de caja negra y adaptativo, se aprovechan de la log concavidad para poder realizar 
simulación sobre $f$ aun cuando no tenemos acceso a esta función directamente o a la locación de la 
moda. Esto los pone en ventaja sobre aceptación-rechazo estándar, ya que si se tratara de 
simular de una distribución $f$ para la cual se desconoce la locación de su moda $m$ entonces
se tendría que hacer un paso previo de optimización, lo cual es costoso computacional mente.

El algoritmo adaptativo ademas nos permite ir mejorando la función envolvente que se utiliza para 
tomar muestras candidatas  en el paso de aceptación rechazo, por lo que cada vez que simulamos una variable $X$ y la rechazamos estamos mejorando 
la capacidad de muestreo del algoritmo, un ejemplo de esto podemos ver en la @fig:area_ratios donde tomamos la densidad $f$, y vemos como 
el cociente del área de la envolvente entre la constante de normalización de $f$ reduce entre mas abscisas se usan para 
construir la envolvente.


#figure(
  image("figures/ratios.svg", width: 50%),
  caption: [
    Relación entre el cociente del área de la envolvente entre la constante de normalización de $f$. Se tomaron conjuntos de abscisas 
  con tamaño $k$, se construyo la función envolvente a partir de estas abscisas y se calculo su área. Se calculo el cociente entre el área 
  de la envolvente y la constante de normalización de $f$. El comportamiento decreciente indica que entre mas grande sea el conjunto de 
  el conjunto de abscisas entonces menor sera la probabilidad de generar una muestra de $g$ normalizada que tenga que ser rechazada.
  ],
)<fig:area_ratios>

La versatilidad de estos algoritmos, especialmente en el caso de ARS donde solo se 
requiere la su posición de que $f$ es log cóncava, da una clara ventaja que compensa su alto costo computacional 
que tienen.

Unos anos despues del articulo original donde se propuso (ARS) Gilks y Walters 
propusieron una version que no requiere la derivada de $f$  
@gilks1992derivative. Tambien se ha propuesto una generalizacion para densidades que no 
son log convacavs @martino2011generalization, la cual requiere de consideraciones extra.
Tambien se ha diseñado un metodo para optimizar
la construcción del conjunto de abscisas inicial para 
aumentar la eficiencia del muestreo @james2024automated.

#bibliography("references.bib")
