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

## Arquitetura (extensibilidade)

`EstudoAterramento` recebe `metodo="ieee80"` (atual). O parametro reserva espaco
para um futuro `metodo="numerico"` (segmentacao de condutores / metodo matricial)
sem alterar a interface publica.

## Normas de referencia

- IEEE Std 80 - Guide for Safety in AC Substation Grounding.
- ABNT NBR 7117 - Medicao da resistividade e determinacao da estratificacao do solo.
- ABNT NBR 15751 - Sistemas de aterramento de subestacoes.

## Limitacoes (fora do escopo desta versao)

- Solver numerico de segmentacao (apenas a porta de extensao esta pronta).
- Geometrias de malha irregulares (o IEEE 80 simplificado assume malha retangular).
- Interface grafica (saida via `imprimir_*` e `exportar` em JSON).
