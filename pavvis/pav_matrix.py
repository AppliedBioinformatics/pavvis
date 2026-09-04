from __future__ import annotations

from pathlib import Path

import pandas as pd


class PavMatrix:
    def __init__(self, pav_path: str | Path, metadata_path: str | Path) -> None:

        # Paths.
        pav_path = Path(pav_path)
        metadata_path = Path(metadata_path)
        
        # Load the PAV data from the user-specified .csv file.
        self.pav: pd.DataFrame = pd.read_csv(pav_path, index_col=0)
        
        # Load the metadata from the user-specified .csv file.
        self.metadata: pd.DataFrame = pd.read_csv(metadata_path, index_col=0)