"""earthsolver - Suite de analise de aterramento.

Estratificacao de solo (metodo de Wenner / NBR 7117) e simulacao de malhas de
aterramento pelo metodo simplificado do IEEE Std 80 (resistencia de aterramento,
tensoes de toque e de passo), no estilo do projeto EletroSolver.
"""

from .solo import ModeloSolo
from .estratificacao import Estratificador
from .malha import Malha, EstudoAterramento
from .numerico import Condutor, Eletrodo, EstudoNumerico, estudo_convergencia
from .dxf import from_dxf

__all__ = [
    "ModeloSolo", "Estratificador", "Malha", "EstudoAterramento",
    "Condutor", "Eletrodo", "EstudoNumerico", "estudo_convergencia", "from_dxf",
]

__version__ = "0.1.0"
