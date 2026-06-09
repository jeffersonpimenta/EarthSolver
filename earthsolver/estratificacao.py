"""Estratificacao do solo a partir de medicoes de resistividade (Wenner).

O Estratificador recebe as medicoes de campo (espacamentos e resistencias do
metodo de Wenner, NBR 7117), gera a curva de resistividade aparente e ajusta um
ModeloSolo de N camadas por Levenberg-Marquardt (escrito a mao, somente numpy),
no mesmo espirito do Newton-Raphson do EletroSolver.
"""

import json

import numpy as np

from .filtros import resistividade_aparente
from .solo import ModeloSolo


class Estratificador:
    """Ajusta um modelo de solo estratificado a uma sondagem de Wenner.

    Parametros (1-based na apresentacao, arrays internos 0-based):
      espacamentos: lista dos espacamentos `a` (m) do arranjo de Wenner.
      resistencias: resistencias medidas R (Ohm); se informadas, a resistividade
        aparente e rho_a = 2*pi*a*R (Wenner igualmente espacado, NBR 7117).
      resistividades: alternativa a `resistencias` - rho_a ja calculada.
    """

    def __init__(self, espacamentos, resistencias=None, resistividades=None):
        a = np.asarray(espacamentos, dtype=float)
        if a.ndim != 1 or a.size < 2:
            raise ValueError("informe ao menos 2 espacamentos")
        if np.any(a <= 0):
            raise ValueError("espacamentos devem ser > 0")
        if np.any(np.diff(a) <= 0):
            raise ValueError("espacamentos devem ser crescentes e sem repeticoes")

        if resistencias is not None and resistividades is not None:
            raise ValueError("informe resistencias OU resistividades, nao ambas")
        if resistencias is not None:
            R = np.asarray(resistencias, dtype=float)
            if R.shape != a.shape:
                raise ValueError("resistencias e espacamentos com tamanhos diferentes")
            if np.any(R <= 0):
                raise ValueError("resistencias devem ser > 0")
            rho_a = 2.0 * np.pi * a * R
        elif resistividades is not None:
            rho_a = np.asarray(resistividades, dtype=float)
            if rho_a.shape != a.shape:
                raise ValueError("resistividades e espacamentos com tamanhos diferentes")
            if np.any(rho_a <= 0):
                raise ValueError("resistividades aparentes devem ser > 0")
        else:
            raise ValueError("informe resistencias ou resistividades")

        self.a = a
        self.rho_a = rho_a
        self.modelo = None      # ModeloSolo ajustado (apos estratificar)
        self.rms = None         # erro RMS relativo do ajuste (%)
        self.n_iter = None
        self.convergiu = False

    # ----------------------------------------------------------- modelo direto
    def modelo_direto(self, modelo: ModeloSolo, a=None):
        """Resistividade aparente prevista pelo modelo para os espacamentos `a`."""
        if a is None:
            a = self.a
        return resistividade_aparente(modelo.rho, modelo.espessura, a)

    # --------------------------------------------------------- parametrizacao
    @staticmethod
    def _empacotar(rho, espessura):
        """Vetor de parametros em log (garante positividade no ajuste)."""
        return np.log(np.concatenate([rho, espessura])) if len(espessura) \
            else np.log(np.asarray(rho, dtype=float))

    @staticmethod
    def _desempacotar(p, n_camadas):
        # Limita os parametros em log para evitar overflow em chutes extremos
        # (faixa fisica ampla: ~1e-3 a ~1e7 Ohm.m / m).
        vals = np.exp(np.clip(p, -7.0, 16.0))
        rho = vals[:n_camadas]
        espessura = vals[n_camadas:]
        return rho, espessura

    def _residuo(self, p, n_camadas):
        rho, espessura = self._desempacotar(p, n_camadas)
        prev = resistividade_aparente(rho, espessura, self.a)
        # Protege contra valores nao-positivos do modelo direto em chutes
        # extremos do multi-start (evita NaN no logaritmo).
        prev = np.maximum(prev, 1e-9)
        return np.log(self.rho_a) - np.log(prev)

    def _jacobiano(self, p, n_camadas, r0, eps=1e-5):
        """Jacobiano numerico (diferencas finitas) do residuo."""
        J = np.zeros((r0.size, p.size))
        for j in range(p.size):
            dp = np.zeros_like(p)
            dp[j] = eps
            J[:, j] = (self._residuo(p + dp, n_camadas) - r0) / eps
        return J

    def _chute_inicial(self, n_camadas):
        """Estimativa inicial das resistividades e espessuras.

        rho0 amostra a propria curva medida em `n_camadas` espacamentos
        (em escala log), capturando minimos/maximos intermediarios (ex.: uma
        camada condutiva enterrada). Espessuras espalhadas geometricamente na
        faixa de espacamentos.
        """
        loga = np.log(self.a)
        alvos = np.linspace(loga[0], loga[-1], n_camadas)
        rho0 = np.interp(alvos, loga, self.rho_a)
        if n_camadas == 1:
            return rho0, np.array([])
        a_min, a_max = self.a[0], self.a[-1]
        bordas = np.geomspace(a_min, a_max, n_camadas)
        espessura0 = np.diff(bordas)
        espessura0 = np.clip(espessura0, a_min * 0.1, None)
        return rho0, espessura0

    # -------------------------------------------------------------- ajuste LM
    def estratificar(self, n_camadas: int, x0=None, max_iter: int = 200,
                     tol: float = 1e-8):
        """Ajusta um ModeloSolo de `n_camadas` por Levenberg-Marquardt.

        x0: ModeloSolo opcional como chute inicial. Retorna o ModeloSolo
        ajustado e preenche self.modelo, self.rms, self.n_iter, self.convergiu.
        """
        if n_camadas < 1:
            raise ValueError("n_camadas deve ser >= 1")
        if n_camadas == 1:
            return self._ajuste_uniforme()

        if x0 is not None:
            rho0, espessura0 = np.asarray(x0.rho), np.asarray(x0.espessura)
        else:
            rho0, espessura0 = self._chute_inicial(n_camadas)

        # Multi-start: o chute principal mais alguns perturbados (a inversao de
        # resistividade tem minimos locais; mantem-se a melhor solucao).
        partidas = [self._empacotar(rho0, espessura0)]
        rng = np.random.default_rng(0)
        for _ in range(4):
            jit = rng.uniform(-0.5, 0.5, size=2 * n_camadas - 1)
            partidas.append(self._empacotar(rho0, espessura0) + jit)

        melhor_p, melhor_custo, melhor_conv, melhor_it = None, np.inf, False, 0
        for p0 in partidas:
            p, custo, convergiu, it = self._lm(p0, n_camadas, max_iter, tol)
            if custo < melhor_custo:
                melhor_p, melhor_custo, melhor_conv, melhor_it = p, custo, convergiu, it

        rho, espessura = self._desempacotar(melhor_p, n_camadas)
        self.modelo = ModeloSolo(rho=rho.tolist(), espessura=espessura.tolist())
        self.n_iter = melhor_it
        self.convergiu = melhor_conv
        self.rms = self._rms(self.modelo)
        return self.modelo

    def _lm(self, p, n_camadas, max_iter, tol):
        """Nucleo Levenberg-Marquardt. Retorna (p, custo, convergiu, n_iter)."""
        lam = 1e-2
        r = self._residuo(p, n_camadas)
        custo = float(r @ r)
        convergiu = False
        it = 0
        for it in range(1, max_iter + 1):
            J = self._jacobiano(p, n_camadas, r)
            JTJ = J.T @ J
            JTr = J.T @ r
            diag = np.diag(np.diag(JTJ))
            melhorou = False
            for _ in range(12):  # busca de amortecimento
                try:
                    dp = np.linalg.solve(JTJ + lam * diag, -JTr)
                except np.linalg.LinAlgError:
                    lam *= 10.0
                    continue
                p_novo = p + dp
                r_novo = self._residuo(p_novo, n_camadas)
                custo_novo = float(r_novo @ r_novo)
                if custo_novo < custo:
                    p, r, custo = p_novo, r_novo, custo_novo
                    lam = max(lam / 3.0, 1e-12)
                    melhorou = True
                    break
                lam *= 3.0
            if not melhorou:
                break
            if np.linalg.norm(dp) < tol:
                convergiu = True
                break
        return p, custo, convergiu, it

    def _ajuste_uniforme(self):
        rho_med = float(np.mean(self.rho_a))
        self.modelo = ModeloSolo(rho=[rho_med], espessura=[])
        self.n_iter = 0
        self.convergiu = True
        self.rms = self._rms(self.modelo)
        return self.modelo

    def _rms(self, modelo: ModeloSolo) -> float:
        prev = self.modelo_direto(modelo)
        erro_rel = (prev - self.rho_a) / self.rho_a
        return float(np.sqrt(np.mean(erro_rel ** 2)) * 100.0)

    def auto_estratificar(self, max_camadas: int = 4):
        """Testa N=1..max_camadas e escolhe o melhor por RMS penalizado.

        A penalizacao por numero de parametros evita sobreajuste (preferir o
        modelo mais simples quando o ganho de RMS e marginal).
        """
        if max_camadas < 1:
            raise ValueError("max_camadas deve ser >= 1")
        melhor = None
        for n in range(1, max_camadas + 1):
            estrat_n = Estratificador(self.a, resistividades=self.rho_a)
            estrat_n.estratificar(n)
            n_param = 2 * n - 1
            score = estrat_n.rms + 0.5 * n_param  # penalizacao leve
            if melhor is None or score < melhor[0]:
                melhor = (score, estrat_n.modelo, estrat_n.rms,
                          estrat_n.n_iter, estrat_n.convergiu)
        _, self.modelo, self.rms, self.n_iter, self.convergiu = melhor
        return self.modelo

    # --------------------------------------------------------------- relatorio
    def imprimir_modelo(self, cd: int = 4) -> None:
        if self.modelo is None:
            raise ValueError("rode estratificar() ou auto_estratificar() antes")
        self.modelo.imprimir_modelo(cd)
        print(f"  RMS do ajuste: {self.rms:.{cd}f} %  "
              f"(iteracoes: {self.n_iter}, convergiu: {self.convergiu})")

    def exportar(self, arquivo: str = "solo.json") -> None:
        if self.modelo is None:
            raise ValueError("rode estratificar() ou auto_estratificar() antes")
        dados = self.modelo.to_dict()
        dados["ajuste"] = {
            "rms_percent": self.rms,
            "n_iter": self.n_iter,
            "convergiu": self.convergiu,
            "espacamentos": self.a.tolist(),
            "rho_aparente": self.rho_a.tolist(),
        }
        with open(arquivo, "w") as f:
            json.dump(dados, f, indent=4)
        print(f"Estratificacao exportada para {arquivo}.")
