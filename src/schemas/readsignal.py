from dataclasses import dataclass

import numpy as np


@dataclass
class ReadSignal:
    name: str
    reverse: bool
    signal: np.ndarray
