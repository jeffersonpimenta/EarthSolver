# Memorial de Metodologia — Cálculo de Sistemas de Aterramento (EarthSolver)

> **Documento técnico — base conceitual para Memorial de Metodologia**
> Descreve a metodologia físico-matemática implementada na engine de cálculo de
> aterramento *EarthSolver*. O texto é fiel ao código: cada equação apresentada
> corresponde a uma rotina efetivamente codificada nos módulos
> `solo.py`, `estratificacao.py`, `filtros.py`, `green.py`, `numerico.py`,
> `malha.py` e `seguranca.py`.

---

## Sumário

- [0. Visão geral: os dois solvers](#0-visão-geral-os-dois-solvers)
- [A. Modelo de solo e estratificação](#a-modelo-de-solo-e-estratificação)
- [B. Discretização da malha (segmentação)](#b-discretização-da-malha-segmentação)
- [C. Parâmetros globais: Rg e GPR](#c-parâmetros-globais-rg-e-gpr)
- [D. Potenciais de superfície: tensão absoluta, de toque e de passo](#d-potenciais-de-superfície-tensão-absoluta-de-toque-e-de-passo)
- [E. Critérios de segurança (IEEE Std 80 / NBR 15751)](#e-critérios-de-segurança-ieee-std-80--nbr-15751)
- [Apêndice: convenções e constantes numéricas](#apêndice-convenções-e-constantes-numéricas)

---

## 0. Visão geral: os dois solvers

A engine oferece **dois caminhos de cálculo independentes**, que compartilham o
mesmo modelo de solo e os mesmos critérios de segurança:

1. **Método simplificado IEEE Std 80** (`malha.py`, classe `EstudoAterramento`).
   Usa fórmulas analíticas fechadas (Sverak para a resistência; fatores
   geométricos $K_m$, $K_s$, $K_i$ para as tensões). O solo é reduzido a uma
   **resistividade uniforme equivalente**. Rápido, adequado a malhas retangulares
   regulares.

2. **Método numérico de segmentação de condutores** (`numerico.py`, classe
   `EstudoNumerico`). Discretiza os condutores em segmentos e resolve a
   distribuição de corrente por um sistema de impedâncias (resistências) mútuas,
   montado a partir da **função de Green do solo estratificado em N camadas**
   (`green.py`). Trata geometria arbitrária e solo multicamada.

Ambos terminam no mesmo conjunto de grandezas de saída — resistência de
aterramento $R_g$, elevação de potencial $V_{GPR}$, tensões de toque $E_m$ e de
passo $E_s$ — e na mesma verificação de segurança (`seguranca.py`).

---

## A. Modelo de solo e estratificação

### A.1 Estrutura do modelo de camadas

O solo é representado por um modelo **horizontalmente estratificado** em $N$
camadas planas e paralelas (`ModeloSolo` em `solo.py`). As $N-1$ camadas
superiores têm espessura finita; a camada de fundo é um **semi-espaço infinito**.

$$
\text{Solo} = \big\{\, (\rho_1, h_1),\ (\rho_2, h_2),\ \dots,\ (\rho_{N-1}, h_{N-1}),\ \rho_N \,\big\}
$$

**Onde:**

- $\rho_i$ — resistividade elétrica da camada $i$, contada do topo ($i=1$) ao
  fundo ($i=N$) $\;[\Omega \cdot \text{m}]$
- $h_i$ — espessura da camada $i$ (apenas para $i = 1, \dots, N-1$) $\;[\text{m}]$
- $N$ — número de camadas (adimensional). A camada $N$ não tem espessura: é
  infinita.

A descontinuidade de resistividade entre duas camadas adjacentes é medida pelo
**coeficiente de reflexão**:

$$
k_i = \frac{\rho_{i+1} - \rho_i}{\rho_{i+1} + \rho_i}, \qquad i = 1, \dots, N-1
$$

**Onde:**

- $k_i$ — coeficiente de reflexão na interface entre as camadas $i$ e $i+1$
  $\;[\text{adimensional}, \; -1 < k_i < 1]$. $k_i > 0$ indica camada inferior
  mais resistiva; $k_i < 0$, mais condutiva.

### A.2 Reduções para os métodos analíticos

O método simplificado IEEE 80 exige um solo de **uma camada**. A engine reduz o
modelo multicamada a uma **resistividade uniforme equivalente** por média
ponderada pelas espessuras (`uniforme_equivalente`):

$$
\rho_{eq} =
\frac{\displaystyle \sum_{i=1}^{N} \rho_i\, p_i}{\displaystyle \sum_{i=1}^{N} p_i},
\qquad
p_i =
\begin{cases}
h_i, & i = 1, \dots, N-1 \\[4pt]
h_{N-1}, & i = N
\end{cases}
$$

**Onde:**

- $\rho_{eq}$ — resistividade uniforme equivalente $\;[\Omega \cdot \text{m}]$
- $p_i$ — peso atribuído à camada $i$ $\;[\text{m}]$. A camada de fundo
  (semi-espaço) recebe peso igual à espessura da última camada finita, $h_{N-1}$,
  limitando a influência de contrastes profundos elevados.

Para um solo de uma camada, $\rho_{eq} = \rho_1$. Há ainda uma redução opcional a
**duas camadas** $(\rho_1, \rho_2^{*}, h_1)$, em que $\rho_2^{*}$ é a média
ponderada das camadas $2 \dots N$.

> Essas reduções são apropriadas ao **método simplificado**. Para contraste de
> resistividade forte, recomenda-se usar o **solver numérico** (seção B), que
> trabalha com o perfil multicamada completo via função de Green.

### A.3 Resistividade aparente medida em campo (Wenner)

A estratificação parte das medições do **método de Wenner** (NBR 7117), em que
quatro eletrodos igualmente espaçados de $a$ são cravados em linha. Para cada
espaçamento $a$ mede-se uma resistência $R$. A **resistividade aparente** é:

$$
\rho_a(a) = 2\pi\, a\, R
$$

**Onde:**

- $\rho_a(a)$ — resistividade aparente para o espaçamento $a$ $\;[\Omega \cdot \text{m}]$
- $a$ — espaçamento entre eletrodos adjacentes do arranjo de Wenner $\;[\text{m}]$
- $R$ — resistência medida pelo terrômetro para esse espaçamento $\;[\Omega]$

A varredura de espaçamentos $a_1 < a_2 < \dots < a_m$ produz a **curva de campo**
$\rho_a(a)$: pequenos $a$ "enxergam" as camadas rasas; grandes $a$ penetram as
profundas.

### A.4 Modelo direto: resistividade aparente teórica (transformada de Hankel)

Para ajustar um modelo às medições é preciso o **problema direto** — prever
$\rho_a(a)$ a partir de $\{\rho_i, h_i\}$. O potencial de uma fonte de corrente
pontual num solo estratificado se resolve no domínio de **Hankel** (simetria
cilíndrica). Para o arranjo de Wenner, a engine (`filtros.py`) usa a forma
numericamente estável:

$$
\rho_a(a) = \rho_1 + 2a \int_{0}^{\infty}
\big[\, T(\lambda) - \rho_1 \,\big]\,
\big[\, J_0(\lambda a) - J_0(2\lambda a) \,\big]\, d\lambda
$$

**Onde:**

- $\lambda$ — número de onda (variável de integração da transformada de Hankel)
  $\;[\text{m}^{-1}]$
- $T(\lambda)$ — **kernel** (função de resistividade transformada) do solo
  estratificado $\;[\Omega \cdot \text{m}]$
- $J_0(\cdot)$ — função de Bessel de primeira espécie, ordem zero
  $\;[\text{adimensional}]$
- $\rho_1$ — resistividade da primeira camada $\;[\Omega \cdot \text{m}]$

A subtração de $\rho_1$ e a diferença $J_0(\lambda a) - J_0(2\lambda a)$ fazem o
integrando decair exponencialmente (pois $T(\lambda) \to \rho_1$ quando
$\lambda h_1 \to \infty$), eliminando a cauda oscilatória lenta e permitindo
integração direta por regra do trapézio.

#### Kernel pela recursão de Pekeris

O kernel é avaliado **de baixo para cima** (do semi-espaço para a superfície)
pela recursão de Pekeris (`kernel` em `filtros.py`):

$$
T_N = \rho_N, \qquad
T_i = \frac{T_{i+1} + \rho_i \tanh(\lambda h_i)}
            {1 + \dfrac{T_{i+1}}{\rho_i}\tanh(\lambda h_i)},
\quad i = N-1, \dots, 1
$$

O valor procurado é $T(\lambda) \equiv T_1$.

**Onde:**

- $T_i$ — kernel acumulado a partir da camada $i$ até o fundo
  $\;[\Omega \cdot \text{m}]$
- $\tanh(\cdot)$ — tangente hiperbólica $\;[\text{adimensional}]$

#### Função de Bessel $J_0$ por aproximação polinomial

Para manter a dependência única de `numpy`, $J_0$ é avaliada pela aproximação de
Abramowitz & Stegun (cap. 9, erro $< 10^{-7}$):

$$
J_0(x) \approx
\begin{cases}
\displaystyle \sum_{n=0}^{6} c_n\, t^n, & |x| \le 3, \quad t = (x/3)^2 \\[10pt]
\displaystyle \frac{f_0(z)}{\sqrt{|x|}}\,\cos\!\big(\theta_0(z)\big),
   & |x| \ge 3, \quad z = 3/|x|
\end{cases}
$$

**Onde:**

- $x = \lambda a$ — argumento da função de Bessel $\;[\text{adimensional}]$
- $c_n$ — coeficientes da série de potências (região $|x|\le 3$)
- $f_0(z)$ — amplitude assintótica; $\theta_0(z)$ — fase assintótica (região
  $|x|\ge 3$) $\;[\text{rad}]$

### A.5 Ajuste de curva por Levenberg-Marquardt (curve fitting)

A estratificação propriamente dita (`estratificacao.py`, classe
`Estratificador`) é um **problema inverso**: encontrar $\{\rho_i, h_i\}$ que faça
a curva teórica $\rho_a^{model}(a)$ casar com a curva medida $\rho_a^{meas}(a)$.
Minimiza-se a soma dos quadrados dos resíduos **em escala logarítmica** (que dá
peso relativo igual a todas as décadas de resistividade):

$$
\min_{\mathbf{p}}\ \Phi(\mathbf{p}) = \sum_{j=1}^{m} r_j(\mathbf{p})^2,
\qquad
r_j(\mathbf{p}) = \ln \rho_a^{meas}(a_j) - \ln \rho_a^{model}(a_j;\ \mathbf{p})
$$

**Onde:**

- $\Phi$ — função-custo (soma dos quadrados dos resíduos) $\;[\text{adimensional}]$
- $r_j$ — resíduo logarítmico no espaçamento $a_j$ $\;[\text{adimensional}]$
- $\mathbf{p}$ — vetor de parâmetros do modelo (ver abaixo)
- $m$ — número de pontos medidos (espaçamentos)

#### Parametrização logarítmica

Para garantir positividade ($\rho_i > 0$, $h_i > 0$) sem restrições explícitas, a
otimização ocorre no **logaritmo** dos parâmetros:

$$
\mathbf{p} = \big[\, \ln\rho_1,\ \dots,\ \ln\rho_N,\ \ln h_1,\ \dots,\ \ln h_{N-1} \,\big]^\top
\in \mathbb{R}^{2N-1}
$$

**Onde:**

- $2N-1$ — número de parâmetros livres: $N$ resistividades $+$ $(N-1)$ espessuras.

#### Iteração de Levenberg-Marquardt

A cada iteração resolve-se o sistema linear que combina Gauss-Newton (rápido
perto do mínimo) com gradiente descendente (robusto longe dele), via parâmetro de
amortecimento $\mu$:

$$
\big(\, \mathbf{J}^\top \mathbf{J} + \mu\, \mathrm{diag}(\mathbf{J}^\top \mathbf{J}) \,\big)\, \Delta\mathbf{p}
= -\, \mathbf{J}^\top \mathbf{r}
$$

$$
\mathbf{p}^{(k+1)} = \mathbf{p}^{(k)} + \Delta\mathbf{p}
$$

**Onde:**

- $\mathbf{J}$ — matriz Jacobiana, $J_{j\ell} = \partial r_j / \partial p_\ell$,
  calculada por **diferenças finitas** $\;[\text{adimensional}]$
- $\mathbf{r}$ — vetor de resíduos $\;[\text{adimensional}]$
- $\mu$ — parâmetro de amortecimento de Marquardt $\;[\text{adimensional}]$. Cresce
  (×3) quando o passo piora o custo e diminui (÷3) quando melhora.
- $\Delta\mathbf{p}$ — incremento dos parâmetros na iteração $\;[\text{adimensional}]$
- $\mathrm{diag}(\mathbf{J}^\top \mathbf{J})$ — matriz diagonal com os elementos
  diagonais de $\mathbf{J}^\top \mathbf{J}$ (escala o amortecimento por parâmetro).

O laço interno de busca de amortecimento aceita o passo se
$\Phi(\mathbf{p}+\Delta\mathbf{p}) < \Phi(\mathbf{p})$; converge quando
$\lVert \Delta\mathbf{p} \rVert < \text{tol}$ (padrão $10^{-8}$).

#### Multi-start (mínimos locais)

A inversão de resistividade tem **mínimos locais**. A engine roda o ajuste a
partir de **5 pontos de partida**: o chute inicial principal (amostragem da
própria curva medida) mais 4 versões perturbadas aleatoriamente. Mantém-se a
solução de menor custo:

$$
\mathbf{p}^{*} = \arg\min_{\text{partidas}}\ \Phi\!\left(\mathbf{p}^{(\infty)}_{\text{partida}}\right)
$$

#### Qualidade do ajuste (RMS) e seleção automática de N

O erro do ajuste é reportado como **RMS relativo percentual**:

$$
\mathrm{RMS} = 100 \times
\sqrt{\frac{1}{m}\sum_{j=1}^{m}
\left(\frac{\rho_a^{model}(a_j) - \rho_a^{meas}(a_j)}{\rho_a^{meas}(a_j)}\right)^{2}}
$$

**Onde:**

- $\mathrm{RMS}$ — desvio quadrático médio relativo do ajuste $\;[\%]$

A rotina `auto_estratificar` testa $N = 1, 2, \dots, N_{max}$ e escolhe o modelo
por um **critério penalizado** (princípio da parcimônia — evita sobreajuste):

$$
\text{score}(N) = \mathrm{RMS}(N) + 0{,}5\,(2N-1)
$$

Vence o $N$ de menor *score*: só se aceita uma camada extra se o ganho de RMS
compensar a penalização de $0{,}5$ por parâmetro adicional.

---

## B. Discretização da malha (segmentação)

### B.1 Conceito: do contínuo ao discreto

O solver numérico (`numerico.py`) modela o sistema de aterramento como uma
coleção de **condutores retos** (`Condutor`): cada condutor é um segmento
$\mathbf{p}_1 \to \mathbf{p}_2$ com raio $a$. Condutores horizontais (cabos do
reticulado) e verticais (hastes) são tratados de forma unificada — uma haste é
apenas um condutor cuja direção tem componente $z$.

A física do problema (equação de Laplace para o potencial no solo, com a corrente
dispersando pela superfície metálica) **não tem solução fechada** para geometria
arbitrária. A estratégia é **discretizar**: dividir cada condutor em $M$ pequenos
**segmentos**, sobre cada qual se admite que a densidade de corrente drenada para
o solo é **constante** (função de base tipo "pulso"). Esta é a essência do
**Método dos Momentos (MoM)**.

### B.2 Regras de segmentação

Cada condutor de comprimento $L_c$ é subdividido em $n_{seg}$ segmentos de
comprimento $\le \ell_{alvo}$ (`segmentar`):

$$
n_{seg} = \left\lceil \frac{L_c}{\ell_{alvo}} \right\rceil
$$

**Onde:**

- $\ell_{alvo}$ — comprimento-alvo de segmento (`comp_alvo`) $\;[\text{m}]$
- $L_c$ — comprimento do condutor $\;[\text{m}]$
- $\lceil \cdot \rceil$ — função teto (arredondamento para cima)

Regra adicional importante para solo estratificado: os pontos de quebra incluem
os cruzamentos das **interfaces de camada** (profundidades acumuladas
$d_i = \sum_{k\le i} h_k$). Assim **cada segmento fica inteiramente dentro de uma
única camada de solo**, e a ele se associa a resistividade local $\rho_{(c)}$ da
camada que o contém. Cada segmento $i$ é representado por:

$$
\text{segmento}_i = \big\{\ \mathbf{m}_i,\ \hat{\mathbf{u}}_i,\ \ell_i/2,\ a_i,\ c_i\ \big\}
$$

**Onde:**

- $\mathbf{m}_i$ — ponto médio do segmento $\;[\text{m}]$ (vetor em $\mathbb{R}^3$)
- $\hat{\mathbf{u}}_i$ — versor (direção unitária) do segmento $\;[\text{adimensional}]$
- $\ell_i/2$ — meio-comprimento do segmento $\;[\text{m}]$
- $a_i$ — raio do condutor do segmento $\;[\text{m}]$
- $c_i$ — índice (0-based) da camada de solo que contém o segmento

### B.3 Convergência numérica (GPR / $R_g$ vs. número de segmentos)

A solução discreta **converge** para a solução exata à medida que
$\ell_{alvo} \to 0$ (ou seja, $M \to \infty$). Como em todo método numérico, há
um compromisso entre exatidão e custo: a montagem da matriz custa
$\mathcal{O}(M^2)$ e a solução do sistema $\mathcal{O}(M^3)$.

A engine quantifica isso com um **estudo de convergência**
(`estudo_convergencia`): resolve o mesmo eletrodo para uma série de
$\ell_{alvo}$ decrescentes e registra $R_g$ e $V_{GPR}$ em função de $M$:

$$
\big\{\, (M_1, R_{g,1}, V_{GPR,1}),\ (M_2, R_{g,2}, V_{GPR,2}),\ \dots \,\big\}
$$

O critério prático de convergência é a estabilização relativa:

$$
\left| \frac{R_g(M_{k+1}) - R_g(M_k)}{R_g(M_k)} \right| < \varepsilon
$$

**Onde:**

- $M_k$ — número de segmentos na discretização $k$ $\;[\text{adimensional}]$
- $\varepsilon$ — tolerância relativa adotada (tipicamente $1\%$–$2\%$)

O gráfico $R_g$ e $V_{GPR}$ vs. $M$ (rotina `plot_convergencia`) mostra a curva
achatando num platô: o engenheiro adota o $M$ a partir do qual o resultado deixa
de variar significativamente. **Validação registrada:** para o exemplo do IEEE
Std 80, o solver converge para $R_g \approx 2{,}65\ \Omega$ (Sverak analítico:
$2{,}78\ \Omega$; diferença $\approx 5\%$).

---

## C. Parâmetros globais: $R_g$ e GPR

### C.1 Método numérico: matriz de resistências e sistema equipotencial

#### Potencial mútuo entre segmentos (montagem da matriz)

O elemento $R_{ij}$ da matriz de resistências é o **potencial médio induzido no
segmento $i$ pela corrente unitária drenada pelo segmento $j$**. Para um segmento
fonte $j$ (extremos $\mathbf{A}_j \to \mathbf{B}_j$, comprimento $L_j$) que drena
corrente uniformemente, o potencial num ponto de campo $\mathbf{P}$ vale, no solo
homogêneo de referência (mais a imagem de superfície — ver C.1.3):

$$
\varphi_j(\mathbf{P}) = \frac{\rho}{4\pi L_j}
\int_{0}^{L_j} \frac{d\ell}{\big|\,\mathbf{P} - (\mathbf{A}_j + \ell\,\hat{\mathbf{u}}_j)\,\big|}
$$

Essa integral de linha tem **forma fechada** (implementada em `_phi_seg`):

$$
\int_{0}^{L} \frac{d\ell}{\sqrt{p^2 + (\ell - s_0)^2}}
= \operatorname{arcsinh}\!\left(\frac{L - s_0}{p}\right)
+ \operatorname{arcsinh}\!\left(\frac{s_0}{p}\right)
$$

**Onde:**

- $\varphi_j$ — potencial gerado pelo segmento $j$ $\;[\text{V}]$ (por corrente
  unitária, é uma resistência $[\Omega]$)
- $\rho$ — resistividade local do segmento fonte $\;[\Omega \cdot \text{m}]$
- $L_j$ — comprimento do segmento fonte $\;[\text{m}]$
- $s_0 = (\mathbf{P} - \mathbf{A}_j)\cdot\hat{\mathbf{u}}_j$ — projeção do ponto de
  campo sobre o eixo do segmento $\;[\text{m}]$
- $p$ — distância perpendicular do ponto de campo ao eixo do segmento, **regularizada
  pelo raio**: $p = \sqrt{p_\perp^2 + a^2}$ $\;[\text{m}]$
- $\operatorname{arcsinh}$ — seno hiperbólico inverso $\;[\text{adimensional}]$

A regularização $p = \sqrt{p_\perp^2 + a^2}$ é o que torna a expressão válida
inclusive para segmentos colineares vizinhos ($p_\perp \to 0$) e dá origem ao
termo próprio logarítmico.

O potencial médio sobre o segmento de campo $i$ (necessário para $R_{ij}$) é
obtido por **quadratura de Gauss-Legendre** de ordem $n_g$ ao longo do segmento:

$$
R_{ij} = \frac{1}{I_j}\,\langle \varphi_j \rangle_{\text{seg } i}
\approx \sum_{g=1}^{n_g} w_g\, \varphi_j(\mathbf{P}_{i,g})
$$

**Onde:**

- $R_{ij}$ — resistência mútua entre os segmentos $i$ e $j$ $\;[\Omega]$
- $w_g$ — pesos de Gauss-Legendre normalizados ($\sum_g w_g = 1$) $\;[\text{adimensional}]$
- $\mathbf{P}_{i,g}$ — nós de Gauss ao longo do segmento de campo $i$ $\;[\text{m}]$
- $n_g$ — ordem da quadratura (`n_gauss`, padrão 4)

#### Termo próprio (diagonal)

Na diagonal ($i = j$) a integral é singular; usa-se a **forma fechada do
auto-potencial** de um segmento reto fino:

$$
R_{ii}^{\,dir} = \frac{\rho}{2\pi L_i}
\left[\, \ln\!\left(\frac{2 L_i}{a_i}\right) - 1 \,\right]
$$

**Onde:**

- $R_{ii}^{\,dir}$ — auto-resistência (parte direta) do segmento $i$ $\;[\Omega]$
- $a_i$ — raio do condutor $\;[\text{m}]$

À diagonal soma-se ainda a contribuição da **imagem** (avaliada por quadratura,
não singular).

#### Imagem de superfície (interface ar/solo)

A superfície do solo ($z=0$) é uma fronteira com o ar isolante: a componente
normal da densidade de corrente é nula ($\partial V/\partial z = 0$). Isso é
imposto exatamente pelo **método das imagens** — uma fonte-espelho de mesmo sinal
refletida em $z \to -z$. Cada distância fonte-campo passa a ter um par:

$$
r = \sqrt{r_h^2 + (z - z')^2}, \qquad
r_{img} = \sqrt{r_h^2 + (z + z')^2}
$$

**Onde:**

- $r$ — distância direta fonte-campo $\;[\text{m}]$
- $r_{img}$ — distância à fonte-imagem $\;[\text{m}]$
- $r_h$ — distância **horizontal** entre fonte e campo $\;[\text{m}]$
- $z, z'$ — profundidades do ponto de campo e da fonte $\;[\text{m}]$ ($z \ge 0$,
  crescente para baixo)

#### Correção de solo estratificado (função de Green N-camadas)

Para solo de **mais de uma camada** ($N>1$), a parte direta+imagem é apenas a
aproximação de solo homogêneo. A engine soma a essa base a **função de Green
completa do solo estratificado** (`green.py`), decomposta como:

$$
G(r_h, z, z') = \underbrace{\frac{\rho_m}{4\pi}\left(\frac{1}{r} + \frac{1}{r_{img}}\right)}_{\text{forma fechada (direta + imagem)}}
\;+\; \underbrace{G_{rem}(r_h, z, z')}_{\text{resto de camadas}}
$$

**Onde:**

- $G$ — potencial completo de uma fonte pontual unitária no solo N-camadas
  $\;[\text{V/A} = \Omega]$
- $\rho_m$ — resistividade da camada que contém a **fonte** $\;[\Omega \cdot \text{m}]$
- $G_{rem}$ — **resto** de camadas: termo limitado (finito mesmo na coincidência
  fonte = campo); é **nulo** para solo uniforme $\;[\Omega]$

O resto é obtido por **transformada de Hankel** de um problema de contorno 1-D
resolvido por número de onda $\lambda$. Em cada camada $j$ o potencial
transformado tem a forma:

$$
\psi_j(\lambda, z) = A_j\, e^{-\lambda (z - \text{topo}_j)} + B_j\, e^{-\lambda (\text{base}_j - z)}
$$

(expoentes $\le 0$ dentro da camada — sem *overflow* numérico). Os coeficientes
$A_j, B_j$ (total $2N-1$ incógnitas) saem de um **sistema linear** montado pelas
condições de contorno:

- continuidade de $\psi$ nas interfaces;
- continuidade de $\tfrac{1}{\rho}\,\partial\psi/\partial z$ (densidade de corrente
  normal) nas interfaces;
- $\partial\psi/\partial z = 0$ na superfície (ar isolante);
- decaimento no semi-espaço de fundo.

E o potencial físico do resto é a integral de Hankel:

$$
G_{rem}(r_h, z, z') = \int_{0}^{\infty}
\big[\, \psi_{obs}(\lambda; z, z') - \psi_{fechada}(\lambda; z, z') \,\big]\,
J_0(\lambda r_h)\, d\lambda
$$

**Onde:**

- $\psi_{obs}$ — potencial transformado na camada de observação (solução do
  problema de contorno)
- $\psi_{fechada}$ — transformada exata da parte direta+imagem (subtraída para
  condicionar a integral; o resultado é **exato**, a subtração só melhora o
  condicionamento)
- $A_j, B_j$ — coeficientes do potencial transformado na camada $j$

A integral é discretizada por regra do trapézio numa grade de $\lambda$ adaptada à
escala de profundidade (decaimento $\sim e^{-\lambda \cdot \text{profundidade}}$).

#### Sistema equipotencial e $R_g$

Reunidos todos os termos, monta-se a matriz simétrica $\mathbf{R}\ (M\times M)$.
O condutor metálico, sendo bom condutor e equipotencializado, força **todos os
segmentos ao mesmo potencial** $V$ (o GPR). A relação potencial-corrente é:

$$
\mathbf{R}\, \mathbf{I} = V\, \mathbf{1}
\quad\Longrightarrow\quad
\mathbf{I} = V\, \mathbf{R}^{-1}\, \mathbf{1}
$$

**Onde:**

- $\mathbf{R}$ — matriz de resistências mútuas/próprias $\;[\Omega]$
- $\mathbf{I}$ — vetor das correntes drenadas para o solo por cada segmento
  $\;[\text{A}]$
- $V$ — potencial comum de todos os segmentos = GPR $\;[\text{V}]$
- $\mathbf{1}$ — vetor unitário (todos os elementos iguais a 1)

A corrente total injetada é $I_g = \mathbf{1}^\top \mathbf{I}
= V\,(\mathbf{1}^\top \mathbf{R}^{-1} \mathbf{1})$. Logo a **resistência de
aterramento** é:

$$
\boxed{\;
R_g = \frac{V}{I_g} = \frac{1}{\mathbf{1}^\top \mathbf{R}^{-1} \mathbf{1}}
\;}
$$

**Onde:**

- $R_g$ — resistência de aterramento da malha $\;[\Omega]$
- $I_g$ — corrente total que escoa da malha para o solo $\;[\text{A}]$

Numericamente: resolve-se $\mathbf{y} = \mathbf{R}^{-1}\mathbf{1}$ (sem inverter a
matriz), e então $R_g = 1/\sum_i y_i$. A **elevação de potencial de terra**
segue de:

$$
\boxed{\; V_{GPR} = I_g\, R_g \;}
$$

**Onde:**

- $V_{GPR}$ — Ground Potential Rise: potencial absoluto da malha em relação à terra
  remota $\;[\text{V}]$
- $I_g$ — corrente de malha injetada (dado de entrada `Ig`) $\;[\text{A}]$

Por fim, a distribuição de corrente por segmento $\mathbf{I} = V\,\mathbf{y}$
revela onde a malha mais dispersa corrente (tipicamente **picos nos cantos**) —
visualizado em `plot_corrente`.

### C.2 Método simplificado: resistência de Sverak (IEEE Std 80)

No caminho analítico (`malha.py`), $R_g$ vem da fórmula de **Sverak**:

$$
R_g = \rho \left[\, \frac{1}{L_T}
+ \frac{1}{\sqrt{20\,A}}
\left(1 + \frac{1}{1 + h\sqrt{20/A}}\right) \right]
$$

**Onde:**

- $R_g$ — resistência de aterramento $\;[\Omega]$
- $\rho$ — resistividade do solo (uniforme equivalente) $\;[\Omega \cdot \text{m}]$
- $L_T$ — comprimento total enterrado (condutores + hastes), $L_T = L_c + L_r$
  $\;[\text{m}]$
- $A$ — área coberta pela malha $\;[\text{m}^2]$
- $h$ — profundidade de enterramento da malha $\;[\text{m}]$

E o GPR analiticamente: $V_{GPR} = I_g\, R_g$, idêntico em forma ao caso numérico.

---

## D. Potenciais de superfície: tensão absoluta, de toque e de passo

### D.1 Potencial absoluto no solo $V(x, y, z)$

Resolvida a distribuição de corrente $\mathbf{I}$, o **potencial absoluto** em
qualquer ponto $\mathbf{P}=(x,y,z)$ é a **superposição** das contribuições de
todos os segmentos (cada um já incluindo imagem e correção de camadas):

$$
V(\mathbf{P}) = \sum_{j=1}^{M} I_j\, G(\mathbf{P}, \text{segmento}_j)
$$

**Onde:**

- $V(\mathbf{P})$ — potencial absoluto no ponto $\mathbf{P}$, referido à terra
  remota $\;[\text{V}]$
- $I_j$ — corrente drenada pelo segmento $j$ $\;[\text{A}]$
- $G(\mathbf{P}, \text{segmento}_j)$ — função de Green: potencial em $\mathbf{P}$
  por corrente unitária do segmento $j$ $\;[\Omega]$

Para os mapas de superfície avalia-se $z = 0$. A engine
(`_calcular_superficie`) monta um **raster** regular $(X, Y)$ cobrindo a malha
mais uma margem, com passo $\Delta$ (`passo_raster`) e folga `margem_raster`:

$$
V(x, y, 0) = \sum_{j=1}^{M} I_j\, G\big((x,y,0),\ \text{segmento}_j\big)
$$

Esse campo é a base dos gráficos de potencial 2D (`plot_potencial`) e 3D
(`plot_potencial_3d`).

### D.2 Tensão de toque ($E_m$)

A **tensão de toque** é a diferença entre o potencial do metal (que a pessoa
segura — o GPR) e o potencial do solo onde ela pisa, a $1\ \text{m}$ de distância.
O pior caso (`Em`) ocorre onde o potencial de superfície **dentro da projeção da
malha** é mínimo:

$$
E_m = V_{GPR} - \min_{(x,y)\in \text{malha}} V(x, y, 0)
$$

**Onde:**

- $E_m$ — tensão de toque (pior caso na área da malha) $\;[\text{V}]$
- $V_{GPR}$ — potencial da malha $\;[\text{V}]$
- $\min V(x,y,0)$ — menor potencial de superfície dentro da projeção do reticulado
  $\;[\text{V}]$

O campo completo $V_{GPR} - V(x,y,0)$ vira o **mapa de tensão de toque**
(`plot_tensao_toque`), com hachura nas regiões que excedem o limite tolerável.

### D.3 Tensão de passo ($E_s$)

A **tensão de passo** é a diferença de potencial entre os dois pés de uma pessoa,
separados por $1\ \text{m}$. A engine avalia, em cada ponto do raster, a diferença
nas direções $x$ e $y$ e toma a pior:

$$
\Delta V_x = V(x{+}1, y, 0) - V(x, y, 0), \qquad
\Delta V_y = V(x, y{+}1, 0) - V(x, y, 0)
$$

$$
E_s = \max_{(x,y)} \ \max\big(\, |\Delta V_x|,\ |\Delta V_y| \,\big)
$$

**Onde:**

- $E_s$ — tensão de passo (pior caso de toda a área) $\;[\text{V}]$
- $\Delta V_x, \Delta V_y$ — variação do potencial de superfície sobre $1\ \text{m}$
  nas direções $x$ e $y$ $\;[\text{V}]$

Diferente do toque, o passo é avaliado em **toda a área** (a pessoa pode pisar em
qualquer ponto, inclusive na periferia onde o gradiente é maior). O campo vira o
**mapa de tensão de passo** (`plot_tensao_passo`).

### D.4 Varreduras paramétricas (espaçamento e profundidade)

A influência do **espaçamento entre condutores** $D$ e da **profundidade** $h$
sobre as curvas de toque e passo é obtida por **varredura paramétrica**:
re-resolve-se o estudo (numérico ou IEEE 80) para uma série de valores de $D$ ou
$h$, montando geometrias com `malha_retangular`. Tendência física capturada pelo
modelo:

- **Diminuir $D$** (malha mais densa) reduz $E_m$ e $E_s$ — superfície mais
  equipotencializada;
- **Aumentar $h$** reduz o gradiente raso, atenuando principalmente $E_s$ na
  periferia.

Os **perfis em corte** pelas linhas centrais (`plot_perfis`) traçam, lado a lado,
o potencial de superfície, a tensão de toque e a de passo ao longo de um eixo,
com as linhas de GPR e dos limites toleráveis — a apresentação clássica do
IEEE 80 para inspeção visual da margem.

### D.5 Método simplificado: tensões de malha e de passo (IEEE 80)

No caminho analítico (`malha.py`), $E_m$ e $E_s$ não vêm de um campo, mas de
**fórmulas com fatores geométricos**:

$$
E_m = \frac{\rho\, I_g\, K_m\, K_i}{L_M}, \qquad
E_s = \frac{\rho\, I_g\, K_s\, K_i}{L_S}
$$

**Onde:**

- $E_m, E_s$ — tensões de malha (toque) e de passo $\;[\text{V}]$
- $\rho$ — resistividade do solo $\;[\Omega \cdot \text{m}]$
- $I_g$ — corrente de malha $\;[\text{A}]$
- $K_m$ — fator geométrico de malha (toque) $\;[\text{adimensional}]$
- $K_s$ — fator geométrico de passo $\;[\text{adimensional}]$
- $K_i$ — fator de irregularidade $\;[\text{adimensional}]$
- $L_M$ — comprimento efetivo para tensão de malha $\;[\text{m}]$
- $L_S$ — comprimento efetivo para tensão de passo $\;[\text{m}]$

Os fatores são (todos codificados em `EstudoAterramento`):

$$
n = n_a\, n_b\, n_c\, n_d, \qquad
n_a = \frac{2 L_c}{L_p}, \quad
n_b = \sqrt{\frac{L_p}{4\sqrt{A}}}, \quad
n_c = n_d = 1 \ (\text{retangular})
$$

$$
K_i = 0{,}644 + 0{,}148\, n
$$

$$
K_m = \frac{1}{2\pi}\left[
\ln\!\left(\frac{D^2}{16\,h\,d} + \frac{(D+2h)^2}{8\,D\,d} - \frac{h}{4d}\right)
+ \frac{K_{ii}}{K_h}\,\ln\!\left(\frac{8}{\pi(2n-1)}\right)
\right]
$$

$$
K_s = \frac{1}{\pi}\left[
\frac{1}{2h} + \frac{1}{D+h} + \frac{1}{D}\left(1 - 0{,}5^{\,n-2}\right)
\right]
$$

$$
K_h = \sqrt{1 + \frac{h}{h_0}}\ (h_0 = 1\,\text{m}), \qquad
K_{ii} = \begin{cases} 1, & \text{com hastes} \\ (2n)^{-2/n}, & \text{sem hastes} \end{cases}
$$

**Onde:**

- $L_p$ — perímetro da malha $\;[\text{m}]$
- $L_c$ — comprimento total de condutores horizontais $\;[\text{m}]$
- $D$ — espaçamento entre condutores paralelos $\;[\text{m}]$
- $h$ — profundidade da malha $\;[\text{m}]$
- $d$ — diâmetro do condutor $\;[\text{m}]$
- $n$ — fator geométrico de número de condutores $\;[\text{adimensional}]$
- $K_h$ — fator de profundidade; $K_{ii}$ — fator de condutores internos;
  $h_0$ — profundidade de referência $\;[\text{m}]$

Os comprimentos efetivos:

$$
L_M = \begin{cases}
L_c + L_r, & \text{sem hastes} \\[4pt]
L_c + \left(1{,}55 + 1{,}22\,\dfrac{L_{haste}}{\sqrt{L_x^2 + L_y^2}}\right) L_r, & \text{com hastes}
\end{cases}
\qquad
L_S = 0{,}75\, L_c + 0{,}85\, L_r
$$

**Onde:**

- $L_r$ — comprimento total de hastes $\;[\text{m}]$
- $L_{haste}$ — comprimento de cada haste $\;[\text{m}]$
- $L_x, L_y$ — dimensões da malha; $\sqrt{L_x^2+L_y^2}$ é a diagonal $\;[\text{m}]$

### D.6 Sobre o "Método dos Elementos Finitos" — esclarecimento

O solver numérico do EarthSolver **não usa FEM**. Ele usa o **Método dos
Momentos (MoM)**, que pertence à família das **equações integrais de contorno**
(BEM). A distinção é essencial num memorial técnico:

| Aspecto | MoM / BEM (**implementado**) | FEM (**não implementado**) |
|---|---|---|
| O que se discretiza | apenas a **superfície dos condutores** (segmentos 1-D) | todo o **volume do solo** (malha 3-D de elementos) |
| Incógnita | densidade de corrente drenada por segmento | potencial $V$ nos nós da malha de volume |
| Domínio infinito do solo | tratado **exatamente** pela função de Green (imagem + camadas) | exige truncar o domínio e impor condições de contorno artificiais |
| Equação resolvida | equação integral de Fredholm (potencial = superposição de Green) | forma fraca de $\nabla\!\cdot(\sigma\nabla V) = 0$ |
| Tamanho do sistema | $M$ = nº de segmentos (pequeno, denso) | nº de nós do volume (grande, esparso) |

Para o problema de aterramento — **solo semi-infinito** com eletrodos finos — o
MoM é a escolha natural e mais eficiente: a condição de radiação ao infinito já
está embutida na função de Green, sem necessidade de malhar e truncar o solo. Uma
formulação FEM resolveria, na forma fraca, a equação de continuidade da corrente
estacionária:

$$
\nabla \cdot \big(\sigma(\mathbf{r})\,\nabla V(\mathbf{r})\big) = 0,
\qquad \sigma = \frac{1}{\rho}
$$

**Onde:**

- $V(\mathbf{r})$ — potencial elétrico no ponto $\mathbf{r}$ do solo $\;[\text{V}]$
- $\sigma(\mathbf{r})$ — condutividade local do solo $\;[\text{S/m}]$
- $\rho$ — resistividade local $\;[\Omega \cdot \text{m}]$

discretizando o **volume** do solo em elementos (tetraedros/hexaedros) com funções
de forma nodais — abordagem distinta da adotada. **Conclusão:** sempre que este
memorial menciona o "solver numérico", trata-se de MoM com função de Green
N-camadas, equivalente em propósito ao FEM mas mais adequado ao domínio aberto do
aterramento.

---

## E. Critérios de segurança (IEEE Std 80 / NBR 15751)

A verificação final compara as tensões **calculadas** ($E_m$, $E_s$) com os
**limites suportáveis pelo corpo humano** (`seguranca.py`), comuns aos dois
solvers.

### E.1 Fator de redução da camada superficial $C_s$

Uma camada superficial de alta resistividade (brita, asfalto) aumenta a
resistência de contato dos pés e, portanto, a tensão suportável. Seu efeito é o
**fator de redução** $C_s$ (ajuste empírico do IEEE 80):

$$
C_s = 1 - \frac{0{,}09\left(1 - \dfrac{\rho}{\rho_s}\right)}{2\,h_s + 0{,}09}
$$

**Onde:**

- $C_s$ — fator de redução da camada superficial $\;[\text{adimensional}, \le 1]$
- $\rho$ — resistividade do solo natural $\;[\Omega \cdot \text{m}]$
- $\rho_s$ — resistividade da camada superficial (brita/asfalto) $\;[\Omega \cdot \text{m}]$
- $h_s$ — espessura da camada superficial $\;[\text{m}]$

Quando não há camada superficial ($\rho_s = \rho$), resulta $C_s = 1$.

### E.2 Corrente tolerável pelo corpo e tensões limite

A corrente máxima tolerável pelo corpo segue o critério de fibrilação de
Dalziel, função da duração do choque:

$$
I_B = \frac{c}{\sqrt{t}}
$$

**Onde:**

- $I_B$ — corrente tolerável pelo corpo $\;[\text{A}]$
- $t$ — duração do choque / da falta $\;[\text{s}]$
- $c$ — constante corporal de energia $\;[\text{A}\cdot\sqrt{\text{s}}]$:
  $c = 0{,}116$ para corpo de $50\ \text{kg}$ e $c = 0{,}157$ para $70\ \text{kg}$.

Modelando o corpo como resistência $R_B = 1000\ \Omega$ e os pés como resistências
de contato com o solo ($R_f = 3\,C_s\,\rho_s$ por pé), obtêm-se as **tensões
limite**. Na **tensão de toque** os dois pés ficam **em paralelo**
($1{,}5\,C_s\,\rho_s$); na **tensão de passo**, **em série** ($6{,}0\,C_s\,\rho_s$):

$$
\boxed{\;
E_{toque}^{lim} = \big(1000 + 1{,}5\, C_s\, \rho_s\big)\,\frac{c}{\sqrt{t}}
\;}
$$

$$
\boxed{\;
E_{passo}^{lim} = \big(1000 + 6{,}0\, C_s\, \rho_s\big)\,\frac{c}{\sqrt{t}}
\;}
$$

**Onde:**

- $E_{toque}^{lim}$ — tensão de toque máxima tolerável $\;[\text{V}]$
- $E_{passo}^{lim}$ — tensão de passo máxima tolerável $\;[\text{V}]$
- $1000$ — resistência do corpo humano $R_B$ $\;[\Omega]$
- $1{,}5\,C_s\,\rho_s$ — resistência dos dois pés em paralelo (toque) $\;[\Omega]$
- $6{,}0\,C_s\,\rho_s$ — resistência dos dois pés em série (passo) $\;[\Omega]$
- $c$ — constante corporal ($0{,}116$ para $50\ \text{kg}$; $0{,}157$ para
  $70\ \text{kg}$)

### E.3 Veredito de aprovação

O sistema é considerado **seguro** quando ambas as tensões calculadas ficam
abaixo dos respectivos limites:

$$
\text{APROVADO} \iff
\big(E_m \le E_{toque}^{lim}\big)
\;\wedge\;
\big(E_s \le E_{passo}^{lim}\big)
$$

**Onde:**

- $E_m$ — tensão de toque calculada (seção D) $\;[\text{V}]$
- $E_s$ — tensão de passo calculada (seção D) $\;[\text{V}]$

A engine reporta, além do veredito booleano, a **margem de utilização** ponto a
ponto (`plot_margem`):

$$
U(x,y) = 100 \times \max\!\left(
\frac{E_{toque}(x,y)}{E_{toque}^{lim}},\
\frac{E_{passo}(x,y)}{E_{passo}^{lim}}
\right)\ [\%]
$$

Valores $U > 100\%$ indicam pontos reprovados, destacados por hachura no mapa de
segurança. **Validação registrada:** para o exemplo do IEEE Std 80, o solver
numérico produz tensão de toque $\approx 732\ \text{V}$ contra $749\ \text{V}$ do
cálculo analítico de referência — diferença inferior a $3\%$.

---

## Apêndice: convenções e constantes numéricas

| Símbolo / parâmetro | Significado | Valor / faixa | Onde no código |
|---|---|---|---|
| $z \ge 0$ | profundidade (cresce para baixo, superfície em $z=0$) | — | `numerico.py`, `green.py` |
| $\ell_{alvo}$ | comprimento-alvo de segmento (`comp_alvo`) | padrão $2{,}0\ \text{m}$ | `EstudoNumerico` |
| $n_g$ | ordem da quadratura de Gauss-Legendre (`n_gauss`) | padrão 4 | `_gauss` |
| $\Delta$ | passo do raster de superfície (`passo_raster`) | padrão $2{,}0\ \text{m}$ | `_calcular_superficie` |
| $N_\lambda$ | nº de pontos da quadratura de Hankel | 4001 | `green.py`, `filtros.py` |
| $L_{lim}$ | limite adimensional de integração ($e^{-12}$) | 12,0 | `green.py`, `filtros.py` |
| $c$ | constante corporal ($I_B = c/\sqrt{t}$) | 0,116 (50 kg) / 0,157 (70 kg) | `seguranca.py` |
| $R_B$ | resistência do corpo humano | $1000\ \Omega$ | `seguranca.py` |
| $h_0$ | profundidade de referência do fator $K_h$ | $1{,}0\ \text{m}$ | `malha.py` |
| tol | tolerância de convergência do LM | $10^{-8}$ | `Estratificador._lm` |

### Mapa módulo → responsabilidade

| Módulo | Responsabilidade |
|---|---|
| `solo.py` | modelo de solo N-camadas, reduções (uniforme/2 camadas), $k_i$ |
| `filtros.py` | modelo direto $\rho_a$ (Pekeris + Hankel), $J_0$ (Abramowitz & Stegun) |
| `estratificacao.py` | inversão Wenner → modelo de solo (Levenberg-Marquardt, multi-start, auto-N) |
| `green.py` | função de Green do solo N-camadas (problema de contorno por $\lambda$ + Hankel) |
| `numerico.py` | solver MoM: segmentação, matriz $\mathbf{R}$, $R_g$, GPR, superfície, toque/passo |
| `malha.py` | método simplificado IEEE 80 (Sverak, $K_m$, $K_s$, $K_i$) |
| `seguranca.py` | $C_s$ e tensões toleráveis (IEEE 80 / NBR 15751) |
| `plot.py` | gráficos: planta, potencial 2D/3D, toque, passo, margem, perfis, convergência |

---

*Documento gerado a partir da engine EarthSolver. As equações refletem a
implementação corrente dos módulos citados; constantes e padrões numéricos podem
ser ajustados nos parâmetros das respectivas classes.*
