# earthsolver

Suite de analise de aterramento em Python (apenas `numpy`), no estilo do projeto
**EletroSolver**. Permite:

1. **Estratificar o solo** a partir de medicoes de resistividade pelo metodo de
   Wenner (NBR 7117), em modelos de N camadas.
2. **Simular uma malha de aterramento** pelo metodo simplificado do **IEEE Std 80**
   (validado tambem pela NBR 15751), extraindo:
   - resistencia de aterramento `Rg` (formula de Sverak);
   - elevacao de potencial de terra `GPR`;
   - tensao de toque `Em` e tensao de passo `Es`, comparadas aos limites
     toleraveis (com fator de camada superficial `Cs`, pesos de 50/70 kg).

## Instalacao

```bash
pip install -e .
```

## Uso como biblioteca

```python
from earthsolver import Estratificador, ModeloSolo, Malha, EstudoAterramento

# 1. Estratificacao a partir de uma sondagem de Wenner (a em m, R em Ohm).
estrat = Estratificador(
    espacamentos=[0.5, 1, 2, 4, 8, 16, 32],
    resistencias=[95.49, 52.84, 28.94, 18.06, 15.92, 18.84, 22.74],
)
solo = estrat.auto_estratificar(max_camadas=4)
estrat.imprimir_modelo()

# 2. Simulacao da malha.
malha = Malha(area=4900, Lc=1540, comprimento_x=70, comprimento_y=70,
              espac_D=7, prof_h=0.5, d=0.01, n_hastes=20, comp_haste=7.5,
              rho_s=2500, h_s=0.102)
estudo = EstudoAterramento(solo, malha, Ig=1908, t=0.5, peso=70)
estudo.imprimir_resultado()
```

## Uso pela linha de comando

```bash
# Estratificar uma sondagem CSV (colunas: spacing,resistance ou spacing,rho_a)
earthsolver estratificar --entrada exemplos/solo_wenner.csv --camadas auto --saida solo.json

# Simular uma malha a partir de solo + geometria
earthsolver malha --solo solo.json --malha malha.json --ig 1908 --t 0.5 --peso 70

# Pipeline completo a partir de um projeto JSON
earthsolver analisar exemplos/projeto.json --exportar saida.json
```

## Metodo

- **Modelo direto do solo:** transformada de resistividade pela recursao de
  Pekeris e integral de Hankel de Wenner avaliada numericamente (Bessel J0 por
  aproximacao de Abramowitz & Stegun). So depende de `numpy`.
- **Estratificacao (inversao):** Levenberg-Marquardt escrito a mao, com residuos
  em `log(rho_a)` e parametros em log (positividade garantida).
- **Malha:** equacoes do IEEE Std 80 (Sverak para `Rg`; fatores `n, Ki, Km, Ks,
  Kii, Kh`; comprimentos efetivos `L_M`, `L_S`).

## Metodo numerico de segmentacao

Alem do metodo simplificado do IEEE 80, ha o **solver numerico de segmentacao**
(`EstudoNumerico`), que dispensa fatores empiricos:

1. divide cada condutor (e haste) em segmentos curtos;
2. monta a matriz de resistencias `R` pela **funcao de Green do solo N-camadas**
   (problema de contorno por transformada de Hankel, recursao tipo Pekeris, com a
   imagem da superficie ar/solo) - `green.py`;
3. impoe o eletrodo equipotencial e resolve `R I = V 1`, `sum(I) = Ig`, obtendo
   `Rg = 1/(1^T R^{-1} 1)` e a distribuicao real de corrente;
4. dela calcula o campo de potencial na superficie e as tensoes de toque/passo
   reais, comparadas aos limites toleraveis (`seguranca.py`).

Vantagens sobre o IEEE 80 simplificado: corrente nao uniforme (cantos drenam mais),
**geometria arbitraria** de condutores e **solo estratificado** na propria Green.

```python
from earthsolver import ModeloSolo, Eletrodo, EstudoNumerico

solo = ModeloSolo(rho=[400.0], espessura=[])           # ou N camadas
eletrodo = Eletrodo.malha_retangular(
    comprimento_x=70, comprimento_y=70, espac_x=7, espac_y=7,
    prof_h=0.5, d=0.01, n_hastes=20, comp_haste=7.5)    # ou condutores explicitos
estudo = EstudoNumerico(solo, eletrodo, Ig=1908, t=0.5, peso=70,
                        comp_alvo=3.5, rho_s=2500, h_s=0.102)
estudo.resolver()
estudo.imprimir_resultado()
estudo.exportar_raster("potencial.json")               # mapa de potencial p/ plotar
```

Geometria explicita de condutores (segmentos `p1->p2`, em metros, z = profundidade):

```python
from earthsolver import Condutor, Eletrodo
eletrodo = Eletrodo([Condutor((0, 0, 0.5), (70, 0, 0.5), raio=0.005), ...])
```

Pela linha de comando:

```bash
earthsolver numerico --eletrodo cond.json --solo solo.json \
    --ig 1908 --t 0.5 --peso 70 --comp-alvo 3.5 \
    --rho-s 2500 --h-s 0.102 --exportar saida.json --raster potencial.json
```

`cond.json` aceita `{"condutores": [{"p1":[..],"p2":[..],"raio":..}, ...]}` ou
`{"malha_retangular": {"comprimento_x":70, ...}}`.

## Arquitetura

`EstudoAterramento` (metodo IEEE 80) recebe a `Malha` agregada. O metodo numerico
usa a classe `EstudoNumerico` com geometria explicita de condutores; chamar
`EstudoAterramento(metodo="numerico")` levanta um erro apontando para ela.

## Normas de referencia

- IEEE Std 80 - Guide for Safety in AC Substation Grounding.
- ABNT NBR 7117 - Medicao da resistividade e determinacao da estratificacao do solo.
- ABNT NBR 15751 - Sistemas de aterramento de subestacoes.

## Limitacoes (fora do escopo desta versao)

- Impedancia interna dos condutores (assume eletrodo equipotencial perfeito).
- Resposta em frequencia / acoplamento indutivo (analise resistiva / DC-equivalente).
- Variacao lateral do solo (apenas estratificacao horizontal em camadas).
- Interface grafica (saida via `imprimir_*`, `exportar` e raster em JSON).
