__version__ = "0.0.1"
__copyright__ = "Copyright 2025"
__license__ = "GPL-3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__ = "jianfeng.sunmt@gmail.com"

from typing import List

import random


class Sample:

    def dirichlet(self, k: int, alpha: float, rng: random.Random) -> List[float]:
        return [rng.gammavariate(alpha, 1.0) for _ in range(k)]

    def lognormal(self, k: int, mu: float, sigma: float, rng: random.Random) -> List[float]:
        import math as _m
        return [_m.exp(rng.normalvariate(mu, sigma)) for _ in range(k)]

    def zipf_like(self, k: int, s: float) -> List[float]:
        return [1.0 / (i ** s) for i in range(1, k + 1)]