"""
ClimData Narrator
=================

Phase 4. Turns a machine-readable ``ClimateSummary`` (the flat dict produced by
``analytics.py``) into a human-readable climate assessment written in the voice
of a climate scientist.

Hard constraints (enforced, not just requested):
  * Use ONLY the supplied statistics — never invent a value.
  * Do NOT perform calculations (no unit conversions, no per-decade rescaling,
    no deriving new numbers). The narrator interprets; it does not compute.
  * Explain trends, uncertainty, model agreement/disagreement, and future
    projections WHEN those statistics are present.

Enforcement: after generation, every numeric token in the prose is checked
against the values actually supplied (``grounding_check``). Any ungrounded
number is fed back for repair; if it cannot be grounded, narration fails closed
rather than emitting hallucinated figures.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional

from pydantic import BaseModel, Field


# ===========================================================================
# Schemas
# ===========================================================================

class ClimateSummary(BaseModel):
    """Flexible container — accepts whatever keys analytics produced."""
    model_config = {"extra": "allow"}

    def as_dict(self) -> Dict:
        return self.model_dump()


class NarrativeAssessment(BaseModel):
    text: str
    grounded: bool
    ungrounded_numbers: List[str] = Field(default_factory=list)
    used_statistics: Dict = Field(default_factory=dict)


# Human-meaning of known statistic keys (steers interpretation, adds no numbers).
_GLOSSARY = {
    "annual_trend": "linear trend of the annual mean, in variable units per year (positive = increase/warming)",
    "summer_trend": "linear trend of the JJA summer mean, in units per year",
    "winter_trend": "linear trend of the DJF winter mean, in units per year",
    "hot_days_change": "change in the annual count of hot days (later period minus earlier period)",
    "overall_mean": "long-term mean over the record",
    "model_bias": "mean model-minus-reference difference (positive = model reads high)",
    "rmse": "root-mean-square error against the reference (lower = better)",
    "mae": "mean absolute error against the reference (lower = better)",
    "correlation": "temporal correlation with the reference (1 = perfect agreement)",
}


def _describe_key(key: str) -> str:
    if key in _GLOSSARY:
        return _GLOSSARY[key]
    if key.endswith("_trend"):
        return "linear trend in units per year (positive = increase)"
    if key.endswith("_change"):
        return "change between later and earlier periods"
    if key.startswith(("ssp", "rcp")) or "projection" in key or "future" in key:
        return "future projection statistic"
    return "supplied statistic (interpret literally)"


# ===========================================================================
# Numeric grounding  (the anti-hallucination guard)
# ===========================================================================

_NUM_RE = re.compile(r"-?\d+(?:\.\d+)?")


def _collect_values(obj) -> List[float]:
    """Recursively gather every numeric value from the summary."""
    out: List[float] = []
    if isinstance(obj, dict):
        for v in obj.values():
            out += _collect_values(v)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            out += _collect_values(v)
    elif isinstance(obj, bool):
        pass
    elif isinstance(obj, (int, float)):
        out.append(float(obj))
    return out


def grounding_check(text: str, summary: Dict, tol: float = 0.05) -> List[str]:
    """
    Return numeric tokens in `text` that do NOT match any supplied value.

    A token matches if it is within `tol` (absolute) or 2% (relative) of a
    supplied value or its magnitude. This deliberately fails any number the
    narrator derived by calculation (e.g. per-decade rescaling).
    """
    allowed = _collect_values(summary)
    allowed = allowed + [abs(v) for v in allowed]
    ungrounded: List[str] = []
    for tok in _NUM_RE.findall(text):
        x = float(tok)
        if not any(abs(x - v) <= max(tol, 0.02 * abs(v)) for v in allowed):
            ungrounded.append(tok)
    return ungrounded


# ===========================================================================
# Prompt
# ===========================================================================

def build_narrator_prompt(summary: Dict) -> str:
    lines = []
    for k, v in summary.items():
        lines.append(f"  - {k} = {v}   # {_describe_key(k)}")
    stats_block = "\n".join(lines) if lines else "  (no statistics supplied)"

    return f"""You are a climate scientist writing the findings section of a
short technical assessment. Interpret the statistics below and explain what they
mean for the climate of the location.

ABSOLUTE RULES
- Use ONLY the statistics provided. NEVER state a number that is not in the list.
- Do NOT calculate, convert units, rescale (e.g. per-year → per-decade), or
  derive new figures. Quote the supplied numbers as given.
- Do NOT mention specific calendar years or dates — none are provided.
- If a category of statistic is absent, simply do not discuss it. Do not
  speculate to fill gaps.

WHAT TO COVER (only for statistics that are present)
- Trends: state direction (increase/decrease) and magnitude using the supplied
  per-year trend values; distinguish annual vs seasonal trends if both given.
- Uncertainty: use rmse / mae to characterise error magnitude; be appropriately
  cautious in wording.
- Model agreement / disagreement: use correlation (closeness of agreement) and
  model_bias (systematic over/under-estimation, with its sign) to judge how well
  the model matches the reference.
- Future projections: discuss only if projection statistics are present.

STYLE
- Objective, measured, scientific. 1–3 short paragraphs. No headings, no bullet
  lists, no recommendations, no narrative embellishment beyond the data.

SUPPLIED STATISTICS
{stats_block}

Write the assessment now."""


# ===========================================================================
# Narration
# ===========================================================================

def narrate(
    summary,
    client=None,
    model: str = "qwen3-30b-a3b-instruct-2507",
    max_retries: int = 2,
) -> NarrativeAssessment:
    """
    ClimateSummary → grounded human-readable assessment.

    `client` is an OpenAI-compatible client (e.g. the SAIA client). Ungrounded
    numbers trigger a repair request; if they persist, returns the assessment
    flagged ``grounded=False`` with the offending tokens (fails closed — the
    caller can reject it) rather than silently publishing invented figures.
    """
    data = summary.as_dict() if isinstance(summary, ClimateSummary) else dict(summary)
    if client is None:
        raise ValueError("Provide an OpenAI-compatible `client` (e.g. the SAIA client).")

    messages = [
        {"role": "system", "content": build_narrator_prompt(data)},
        {"role": "user", "content": "Produce the assessment."},
    ]
    text, ungrounded = "", []
    for _ in range(max_retries + 1):
        resp = client.chat.completions.create(
            model=model, messages=messages, temperature=0.2,
        )
        text = resp.choices[0].message.content.strip()
        ungrounded = grounding_check(text, data)
        if not ungrounded:
            return NarrativeAssessment(text=text, grounded=True, used_statistics=data)
        messages += [
            {"role": "assistant", "content": text},
            {"role": "user", "content":
                f"These numbers are NOT in the supplied statistics: {ungrounded}. "
                f"Rewrite using ONLY the supplied values; remove any derived or "
                f"converted figures. Do not calculate."},
        ]
    return NarrativeAssessment(
        text=text, grounded=False, ungrounded_numbers=ungrounded, used_statistics=data
    )


if __name__ == "__main__":
    # Offline test of the grounding guard (no LLM needed).
    summary = {
        "annual_trend": 0.042, "summer_trend": 0.0432,
        "hot_days_change": 28, "model_bias": 0.2,
        "rmse": 0.2, "mae": 0.2, "correlation": 1.0,
    }
    good = ("The annual mean shows a warming trend of 0.042 units per year, with "
            "a comparable summer trend of 0.0432. Hot days increased by 28. The "
            "model agrees closely with the reference (correlation 1.0) with a "
            "small high bias of 0.2 and an rmse of 0.2.")
    bad = ("Warming reaches 0.42 per decade and temperatures rose by 1.5 degrees "
           "since 1980.")  # 0.42, 1.5, 1980 are all invented / calculated
    print("good ungrounded:", grounding_check(good, summary))
    print("bad  ungrounded:", grounding_check(bad, summary))
