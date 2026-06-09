"""Modelo de solo estratificado em camadas horizontais.

ModeloSolo guarda as resistividades e espessuras das camadas e oferece reducoes
(uniforme equivalente, duas camadas) usadas para alimentar as formulas de malha
do IEEE Std 80.
"""

import json
from dataclasses import dataclass, field


@dataclass
class ModeloSolo:
    """Solo horizontalmente estratificado.

    rho: resistividades das camadas em Ohm.m (rho_1 = topo .. rho_N = fundo).
    espessura: espessuras em metros das N-1 camadas superiores; a ultima camada
        e um semi-espaco infinito (por isso len(espessura) == len(rho) - 1).
    """

    rho: list = field(default_factory=list)
    espessura: list = field(default_factory=list)

    def __post_init__(self):
        self.rho = [float(r) for r in self.rho]
        self.espessura = [float(h) for h in self.espessura]
        if len(self.rho) < 1:
            raise ValueError("e necessaria ao menos uma camada (rho nao vazio)")
        if len(self.espessura) != len(self.rho) - 1:
            raise ValueError(
                f"espessura deve ter {len(self.rho) - 1} valores "
                f"(len(rho)-1); recebidos {len(self.espessura)}"
            )
        if any(r <= 0 for r in self.rho):
            raise ValueError("todas as resistividades devem ser > 0")
        if any(h <= 0 for h in self.espessura):
            raise ValueError("todas as espessuras devem ser > 0")

    @property
    def n_camadas(self) -> int:
        return len(self.rho)

    def coef_reflexao(self):
        """Coeficientes de reflexao k_i entre camadas adjacentes."""
        return [
            (self.rho[i + 1] - self.rho[i]) / (self.rho[i + 1] + self.rho[i])
            for i in range(self.n_camadas - 1)
        ]

    def uniforme_equivalente(self) -> float:
        """Resistividade uniforme equivalente (Ohm.m).

        Para solo de uma camada retorna rho_1. Caso contrario usa a media das
        camadas ponderada pela espessura; a camada de fundo (semi-espaco) recebe
        um peso igual a espessura total das camadas superiores, representando a
        profundidade de influencia tipica de uma malha. Simplificacao adequada
        ao metodo simplificado do IEEE 80; para contraste forte prefira informar
        o solo ja como uniforme ou de duas camadas.
        """
        if self.n_camadas == 1:
            return self.rho[0]
        # A camada de fundo (semi-espaco) recebe peso igual a ultima espessura
        # finita, limitando a influencia de contrastes profundos elevados.
        pesos = list(self.espessura) + [self.espessura[-1]]
        return sum(r * p for r, p in zip(self.rho, pesos)) / sum(pesos)

    def duas_camadas(self):
        """Reduz o modelo a (rho_1, rho_2, h) para metodos de duas camadas.

        Mantem a primeira camada e agrega as demais numa resistividade de fundo
        ponderada pela espessura (a camada infinita recebe peso igual a soma das
        espessuras restantes).
        """
        if self.n_camadas == 1:
            return self.rho[0], self.rho[0], None
        if self.n_camadas == 2:
            return self.rho[0], self.rho[1], self.espessura[0]
        h = self.espessura[0]
        sup = self.rho[1:]
        esp = self.espessura[1:]
        pesos = list(esp) + [esp[-1]]
        rho2 = sum(r * p for r, p in zip(sup, pesos)) / sum(pesos)
        return self.rho[0], rho2, h

    def imprimir_modelo(self, cd: int = 4) -> None:
        print("Modelo de Solo Estratificado:")
        print("-" * 40)
        for i in range(self.n_camadas):
            if i < self.n_camadas - 1:
                print(f"  Camada {i + 1}: rho = {self.rho[i]:.{cd}f} Ohm.m, "
                      f"espessura = {self.espessura[i]:.{cd}f} m")
            else:
                print(f"  Camada {i + 1}: rho = {self.rho[i]:.{cd}f} Ohm.m, "
                      f"espessura = infinita")
        print("-" * 40)

    def to_dict(self) -> dict:
        return {"rho": list(self.rho), "espessura": list(self.espessura)}

    @classmethod
    def from_dict(cls, d: dict) -> "ModeloSolo":
        return cls(rho=d["rho"], espessura=d.get("espessura", []))

    def exportar(self, arquivo: str = "solo.json") -> None:
        with open(arquivo, "w") as f:
            json.dump(self.to_dict(), f, indent=4)
        print(f"Modelo de solo exportado para {arquivo}.")
