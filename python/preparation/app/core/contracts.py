# Extracted from qAOPweaver/app/core/contracts.py at 0143dd56dd21a5268b17ecbe8383647b1c6323c9.

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, Field

class ActionLogEntry(BaseModel):
    process: str
    action: str
    rows_before: int
    rows_after: int
    rows_affected: int
    details: str = ""

class HandoffCheck(BaseModel):
    check_id: str
    status: bool
    notes: str

class PreprocessSummary(BaseModel):
    n_rows_before: int
    n_rows_after: int
    n_cols: int
    missing_handled: int
    missing_response_rows: int = 0
    outliers_removed: int
    normalization: str
    column_stats: dict[str, Any]    # per-column describe() keyed by column name
