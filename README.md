# earthsolver

Suíte de análise de aterramento em Python, no estilo do projeto **EletroSolver**.
O núcleo de cálculo depende só de `numpy`; a importação de DXF e os gráficos usam
`ezdxf` e `matplotlib` (instalados junto). Permite:

1. **Estratificar o solo** a partir de medições de resistividade pelo método de
   Wenner (NBR 7117), em modelos de N camadas.
2. **Simular uma malha de aterramento** pelo método simplificado do **IEEE Std 80**
   (validado também pela NBR 15751), extraindo:
   - resistência de aterramento `Rg` (fórmula de Sverak);
   - elevação de potencial de terra `GPR`;
   - tensão de toque `Em` e tensão de passo `Es`, comparadas aos limites
     toleráveis (com fator de camada superficial `Cs`, pesos de 50/70 kg).

## Instalação

```bash
pip install -e .
```

## Uso como biblioteca

```python
from earthsolver import Estratificador, ModeloSolo, Malha, EstudoAterramento

# 1. Estratificação a partir de uma sondagem de Wenner (a em m, R em Ohm).
estrat = Estratificador(
    espacamentos=[0.5, 1, 2, 4, 8, 16, 32],
    resistencias=[95.49, 52.84, 28.94, 18.06, 15.92, 18.84, 22.74],
)
solo = estrat.auto_estratificar(max_camadas=4)
estrat.imprimir_modelo()

# 2. Simulação da malha.
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

## Método

- **Modelo direto do solo:** transformada de resistividade pela recursão de
  Pekeris e integral de Hankel de Wenner avaliada numericamente (Bessel J0 por
  aproximação de Abramowitz & Stegun). Só depende de `numpy`.
- **Estratificação (inversão):** Levenberg-Marquardt escrito à mão, com resíduos
  em `log(rho_a)` e parâmetros em log (positividade garantida).
- **Malha:** equações do IEEE Std 80 (Sverak para `Rg`; fatores `n, Ki, Km, Ks,
  Kii, Kh`; comprimentos efetivos `L_M`, `L_S`).

## Método numérico de segmentação

Além do método simplificado do IEEE 80, há o **solver numérico de segmentação**
(`EstudoNumerico`), que dispensa fatores empíricos:

1. divide cada condutor (e haste) em segmentos curtos;
2. monta a matriz de resistências `R` pela **função de Green do solo N-camadas**
   (problema de contorno por transformada de Hankel, recursão tipo Pekeris, com a
   imagem da superfície ar/solo) - `green.py`;
3. impõe o eletrodo equipotencial e resolve `R I = V 1`, `sum(I) = Ig`, obtendo
   `Rg = 1/(1^T R^{-1} 1)` e a distribuição real de corrente;
4. dela calcula o campo de potencial na superfície e as tensões de toque/passo
   reais, comparadas aos limites toleráveis (`seguranca.py`).

Vantagens sobre o IEEE 80 simplificado: corrente não uniforme (cantos drenam mais),
**geometria arbitrária** de condutores e **solo estratificado** na própria Green.

```python
from earthsolver import ModeloSolo, Eletrodo, EstudoNumerico

solo = ModeloSolo(rho=[400.0], espessura=[])           # ou N camadas
eletrodo = Eletrodo.malha_retangular(
    comprimento_x=70, comprimento_y=70, espac_x=7, espac_y=7,
    prof_h=0.5, d=0.01, n_hastes=20, comp_haste=7.5)    # ou condutores explícitos
estudo = EstudoNumerico(solo, eletrodo, Ig=1908, t=0.5, peso=70,
                        comp_alvo=3.5, rho_s=2500, h_s=0.102)
estudo.resolver()
estudo.imprimir_resultado()
estudo.exportar_raster("potencial.json")               # mapa de potencial p/ plotar
```

Geometria explícita de condutores (segmentos `p1->p2`, em metros, z = profundidade):

```python
from earthsolver import Condutor, Eletrodo
eletrodo = Eletrodo([Condutor((0, 0, 0.5), (70, 0, 0.5), raio=0.005), ...])
```

Pela linha de comando (com `--plot` gera PNGs da malha e do mapa de potencial;
`--plot-3d`/`--plot-malha-3d` geram as vistas em perspectiva 3D):

```bash
earthsolver numerico --eletrodo cond.json --solo solo.json \
    --ig 1908 --t 0.5 --peso 70 --comp-alvo 3.5 \
    --rho-s 2500 --h-s 0.102 --exportar saida.json --raster potencial.json \
    --plot potencial.png --plot-malha malha.png \
    --plot-3d potencial_3d.png --plot-malha-3d malha_3d.png
```

`cond.json` aceita `{"condutores": [{"p1":[..],"p2":[..],"raio":..}, ...]}` ou
`{"malha_retangular": {"comprimento_x":70, ...}}`.

## Importar uma malha de um DXF

Em vez de escrever a geometria à mão, importe o desenho da malha (AutoCAD etc.).
O DXF é plano (só `x, y`); a profundidade de enterramento, o raio do condutor e
quais layers são hastes verticais vêm de um **mapa de layers**:

```bash
# 1. Converter o DXF em geometria. Sem --mapa, um assistente escaneia as layers
#    e pergunta o que é cada uma (condutor/haste, profundidade, raio, comprimento).
earthsolver dxf exemplos/malha.dxf --salvar-mapa layers.json -o cond.json

# 2. Reaproveitando um mapa já salvo (não interativo):
earthsolver dxf exemplos/malha.dxf --mapa layers.json -o cond.json --plot-malha malha.png

# 3. Simular (o cond.json gerado entra direto no comando numerico):
earthsolver numerico --eletrodo cond.json --solo exemplos/solo_uniforme.json \
    --ig 1908 --t 0.5 --comp-alvo 7 --plot potencial.png
```

Mapa de layers (`exemplos/layers.json`):

```json
{
  "escala": 1.0,
  "padrao": { "prof": 0.5, "raio": 0.005 },
  "layers": {
    "MALHA":  { "prof": 0.5, "raio": 0.005 },
    "HASTES": { "rod": true, "prof": 0.5, "comp": 7.5, "raio": 0.005 }
  }
}
```

- `escala`: fator para metros (desenho em mm -> `0.001`; também aceito via `--escala`).
- Condutores vêm de `LINE` / `LWPOLYLINE` / `POLYLINE`; em layer com `rod: true`,
  `POINT` / `INSERT` / `CIRCLE` viram hastes verticais (de `prof` a `prof+comp`).
- `ARC` / `SPLINE` / `ELLIPSE` ficam fora do escopo desta versão (são ignorados com aviso).

## Formatos de arquivo (entradas e saídas)

Fluxo típico - importar a malha de um DXF, simular e gerar os gráficos:

```bash
earthsolver dxf exemplos/malha.dxf --mapa exemplos/layers.json -o cond.json
earthsolver numerico --eletrodo cond.json --solo exemplos/solo_uniforme.json \
    --ig 1908 --t 0.5 --comp-alvo 7 --plot potencial.png --plot-malha malha.png
```

### Entradas (JSON)

**Solo** (`--solo`):
```json
{ "rho": [400.0], "espessura": [] }
```
`rho`: N resistividades em Ohm.m (topo -> fundo). `espessura`: N-1 espessuras em m
(a última camada é um semi-espaço infinito). Solo uniforme -> `espessura: []`.

**Geometria do eletrodo** (`--eletrodo`, comando `numerico`) - duas formas. Explícita
(é o que o comando `dxf` gera), com `p1`/`p2` = `[x, y, z]` em m (`z` = profundidade
>= 0) e `raio` em m:
```json
{ "condutores": [ { "p1": [0, 0, 0.5], "p2": [70, 0, 0.5], "raio": 0.005 } ] }
```
Ou malha retangular gerada automaticamente:
```json
{ "malha_retangular": { "comprimento_x": 70, "comprimento_y": 70, "espac_x": 7,
  "espac_y": 7, "prof_h": 0.5, "d": 0.01, "n_hastes": 20, "comp_haste": 7.5 } }
```

**Malha IEEE 80** (`--malha`, comando `malha`):
```json
{ "area": 4900, "Lc": 1540, "comprimento_x": 70, "comprimento_y": 70, "espac_D": 7,
  "prof_h": 0.5, "d": 0.01, "n_hastes": 20, "comp_haste": 7.5,
  "rho_s": 2500, "h_s": 0.102 }
```
`area` (m2), `Lc` (comprimento total de condutores horizontais, m), `espac_D`
(espaçamento entre condutores paralelos, m), `d` (diâmetro, m). Opcionais:
`n_hastes` (0), `comp_haste` (0), `rho_s` (brita; ausente -> Cs = 1), `h_s` (0.1).

**Projeto completo** (comando `analisar`):
```json
{
  "sondagem": { "espacamentos": [0.5, 1, 2, 4, 8, 16, 32],
                "resistencias": [95.5, 52.8, 28.9, 18.1, 15.9, 18.8, 22.7],
                "camadas": "auto", "max_camadas": 4 },
  "malha":    { "...": "campos da Malha IEEE 80 acima" },
  "falta":    { "Ig": 1908, "t": 0.5, "peso": 70 }
}
```
Use `"solo": {...}` no lugar de `"sondagem"` para informar o solo já pronto; dentro
de `sondagem`, `resistividades` substitui `resistencias` (resistividade aparente direta).

**Mapa de layers do DXF** (`--mapa`, comando `dxf`): ver
[Importar uma malha de um DXF](#importar-uma-malha-de-um-dxf).

> Sondagem por CSV (comando `estratificar --entrada`): não é JSON - colunas
> `spacing,resistance` ou `spacing,rho_a`, com cabeçalho opcional.

### Saídas

**Console** (todos os estudos imprimem um relatório):
```
Estudo de Aterramento (metodo numerico de segmentacao):
--------------------------------------------------------
  Segmentos                 : 260
  Resistencia de malha  Rg  : 2.50 Ohm
  Elevacao de potencial GPR : 4774.46 V
--------------------------------------------------------
  Tensao de toque  Em       : 734.33 V (limite 840.55 V) -> OK
  Tensao de passo  Es       : 382.12 V (limite 2696.10 V) -> OK
--------------------------------------------------------
  Veredito: APROVADO (peso 70 kg, t = 0.5 s)
```

**Resultado** (`--exportar saida.json`):
```json
{
  "entrada": { "Ig": 1908, "t": 0.5, "peso": 70, "metodo": "numerico",
               "rho_s": 2500, "h_s": 0.102, "n_condutores": 42,
               "solo": { "rho": [400.0], "espessura": [] } },
  "resultado": { "Rg": 2.50, "GPR": 4774.46, "V": 4774.46,
                 "Em": 734.33, "Es": 382.12,
                 "E_toque": 840.55, "E_passo": 2696.10, "Cs": 0.74,
                 "toque_ok": true, "passo_ok": true, "aprovado": true,
                 "n_segmentos": 260, "rho_eq": 400.0 }
}
```
No comando `malha` (IEEE 80), `resultado` traz ainda os fatores empíricos
`n, Ki, Km, Ks, Kh, Kii, L_M, L_S`.

**Raster do potencial** (`--raster potencial.json`): matrizes 2D (grade de pontos)
do potencial de superfície em V, para plotar por fora.
```json
{ "x": [[...]], "y": [[...]], "phi": [[...]], "GPR": 4774.46 }
```

**Gráficos** (`--plot potencial.png`, `--plot-malha malha.png`): PNG do mapa de
potencial de superfície (com o contorno da malha sobreposto) e da vista em planta.
**Gráficos 3D** (`--plot-3d potencial_3d.png`, `--plot-malha-3d malha_3d.png`):
perspectiva da malha (condutores enterrados, hastes verticais) e a elevação do
potencial — superfície 3D onde a altura é o potencial em V, com a malha desenhada
no plano da base.

**Modelo de solo** (`estratificar --saida solo.json`): o solo mais um bloco
`ajuste` com a qualidade da inversão:
```json
{ "rho": ["..."], "espessura": ["..."],
  "ajuste": { "rms_percent": 1.2, "n_iter": 8, "convergiu": true,
              "espacamentos": ["..."], "rho_aparente": ["..."] } }
```

## Arquitetura

`EstudoAterramento` (método IEEE 80) recebe a `Malha` agregada. O método numérico
usa a classe `EstudoNumerico` com geometria explícita de condutores; chamar
`EstudoAterramento(metodo="numerico")` levanta um erro apontando para ela.

## Normas de referência

- IEEE Std 80 - Guide for Safety in AC Substation Grounding.
- ABNT NBR 7117 - Medição da resistividade e determinação da estratificação do solo.
- ABNT NBR 15751 - Sistemas de aterramento de subestações.

## Limitações (fora do escopo desta versão)

- Impedância interna dos condutores (assume eletrodo equipotencial perfeito).
- Resposta em frequência / acoplamento indutivo (análise resistiva / DC-equivalente).
- Variação lateral do solo (apenas estratificação horizontal em camadas).
- App/GUI interativo (há gráficos PNG via `--plot` e importação de DXF pela CLI,
  mas a saída continua por `imprimir_*`, `exportar`, raster JSON e PNG).
