# Extracted from qAOPweaver/app/upstream/profile.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from typing import Literal

from pydantic import BaseModel, Field


class ReplicateConfig(BaseModel):
    """Pairing parameters used by the case and its concentration-held-out folds."""

    mode: Literal["repeated_one_to_one"] = "repeated_one_to_one"
    n_bootstrap: int = Field(default=3, ge=1)
    bootstrap_seed: int = 42
