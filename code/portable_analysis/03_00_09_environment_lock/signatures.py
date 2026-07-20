"""Frozen, performance-blind host-RNA score functions for Stage 2.

These functions implement only primary-source arithmetic. They do not fit
weights, choose probes by outcome, optimize cutoffs, or impute absent genes.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from math import exp, isfinite, log
from statistics import fmean
from typing import Any


class SignatureInputError(ValueError):
    """Base class for invalid signature input."""


class MissingGeneError(SignatureInputError):
    """A required gene is absent."""


class DuplicateGeneError(SignatureInputError):
    """A long-form input contains a duplicate gene without an explicit rule."""


class NonNumericExpressionError(SignatureInputError):
    """An expression value is nonnumeric, missing, or nonfinite."""


HISTORICAL_TO_CURRENT = {
    "C9ORF95": "NMRK1",
    "KIAA1370": "FAM214A",
    "FCMR": "FAIM3",
    "C9ORF103": "IDNK",
}


def _canonical(symbol: str) -> str:
    text = str(symbol).strip().upper()
    return HISTORICAL_TO_CURRENT.get(text, text)


def _numeric(value: Any, gene: str) -> float:
    if isinstance(value, bool):
        raise NonNumericExpressionError(f"{gene}: boolean is not an expression value")
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise NonNumericExpressionError(f"{gene}: nonnumeric expression value") from exc
    if not isfinite(out):
        raise NonNumericExpressionError(f"{gene}: missing or nonfinite expression value")
    return out


def prepare_expression(
    expression: Mapping[str, Any] | Iterable[tuple[str, Any]],
    *,
    multi_probe_rule: str = "mean",
) -> dict[str, float]:
    """Canonicalize gene symbols and validate expression.

    Mapping values may be a numeric scalar or an explicit sequence of probe
    values. Probe sequences are summarized only under the declared mean rule.
    Long-form repeated gene rows are rejected to avoid silent aggregation.
    """

    if isinstance(expression, Mapping):
        items = list(expression.items())
    else:
        items = list(expression)

    out: dict[str, float] = {}
    for raw_gene, raw_value in items:
        gene = _canonical(raw_gene)
        if gene in out:
            raise DuplicateGeneError(f"duplicate gene after alias mapping: {gene}")
        if isinstance(raw_value, Sequence) and not isinstance(raw_value, (str, bytes, bytearray)):
            if multi_probe_rule != "mean":
                raise SignatureInputError(f"unsupported multi-probe rule: {multi_probe_rule}")
            values = [_numeric(v, gene) for v in raw_value]
            if not values:
                raise NonNumericExpressionError(f"{gene}: empty probe list")
            out[gene] = fmean(values)
        else:
            out[gene] = _numeric(raw_value, gene)
    return out


def _required(data: Mapping[str, float], genes: Sequence[str]) -> list[float]:
    missing = [_canonical(g) for g in genes if _canonical(g) not in data]
    if missing:
        raise MissingGeneError("missing required gene(s): " + ", ".join(missing))
    return [data[_canonical(g)] for g in genes]


def _geometric_mean(values: Sequence[float], label: str) -> float:
    if any(v <= 0 for v in values):
        raise SignatureInputError(f"{label}: geometric mean requires positive values")
    return exp(fmean([log(v) for v in values]))


SMS_UP = ("CEACAM1", "ZDHHC19", "NMRK1", "GNA15", "BATF", "C3AR1")
SMS_DOWN = ("FAM214A", "TGFBI", "MTCH1", "RPGRIP1", "HLA-DPB1")


def sepsis_metascore_dgm(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """Published Sepsis MetaScore: GM(up six) - (5/6) * GM(down five).

    The original within-dataset z standardization is intentionally not done
    here because no immutable single-sample reference mean/SD was published.
    Higher values indicate infection/sepsis.
    """

    data = prepare_expression(expression)
    return _geometric_mean(_required(data, SMS_UP), "SMS up") - (5.0 / 6.0) * _geometric_mean(
        _required(data, SMS_DOWN), "SMS down"
    )


def septicyte_lab_microarray(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """SeptiCyte LAB microarray arithmetic on log intensity.

    (PLAC8 + LAMP1) - (PLA2G7 + CEACAM4); higher indicates sepsis.
    """

    data = prepare_expression(expression)
    plac8, lamp1, pla2g7, ceacam4 = _required(data, ("PLAC8", "LAMP1", "PLA2G7", "CEACAM4"))
    return plac8 + lamp1 - pla2g7 - ceacam4


def faim3_plac8_ratio(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """Original FAIM3:PLAC8 ratio; higher favors no-CAP/noninfected."""

    data = prepare_expression(expression)
    faim3, plac8 = _required(data, ("FAIM3", "PLAC8"))
    if plac8 == 0:
        raise SignatureInputError("PLAC8 denominator is zero")
    return faim3 / plac8


def faim3_plac8_sepsis_positive(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """Explicit harmonized orientation: negative original FAIM3:PLAC8 ratio."""

    return -faim3_plac8_ratio(expression)


def snip_score(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """sNIP = (NLRP1 - IDNK) / PLAC8; higher indicates abdominal sepsis."""

    data = prepare_expression(expression)
    nlrp1, idnk, plac8 = _required(data, ("NLRP1", "IDNK", "PLAC8"))
    if plac8 == 0:
        raise SignatureInputError("PLAC8 denominator is zero")
    return (nlrp1 - idnk) / plac8


def herberg_two_gene_drs(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """FAM89A - IFI44L; higher favors bacterial rather than viral infection."""

    data = prepare_expression(expression)
    fam89a, ifi44l = _required(data, ("FAM89A", "IFI44L"))
    return fam89a - ifi44l


SWEENEY_BV_VIRAL = ("IFI27", "JUP", "LAX1")
SWEENEY_BV_BACTERIAL = ("HK3", "TNIP1", "GPAA1", "CTSB")


def sweeney_bv_ordering_score(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """GM(viral three) - GM(bacterial four); positive favors viral infection.

    The publication's unresolved positive gene-count scalar is omitted because
    it does not change ordering. Therefore this is A2 discrimination-only.
    """

    data = prepare_expression(expression)
    return _geometric_mean(_required(data, SWEENEY_BV_VIRAL), "BV viral") - _geometric_mean(
        _required(data, SWEENEY_BV_BACTERIAL), "BV bacterial"
    )


SEPSIGLNC_PAIRS = (
    ("ECRP", "CTD-2012K14.6"), ("ECRP", "LOC101926943"), ("ECRP", "MCM3AP-AS1"),
    ("ECRP", "MGC27345"), ("ECRP", "STARD7-AS1"), ("C5ORF66", "AC090627.1"),
    ("GUSBP11", "AC090627.1"), ("LOC284112", "AC090627.1"), ("SDCBP2-AS1", "AC090627.1"),
    ("LOC101927974", "ITFG2"), ("LOC101927974", "LOC100506990"),
    ("LOC101927974", "MIRLET7DHG"), ("AC008753.4", "LOC101928817"),
    ("RP11-533E19.7", "AC100830.4"),
)


def sepsiglnc_pair_count(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """Count 14 sepsis-direction lncRNA orderings; ties count as a point."""

    data = prepare_expression(expression)
    required = tuple(dict.fromkeys(g for pair in SEPSIGLNC_PAIRS for g in pair))
    _required(data, required)
    return float(sum(data[_canonical(left)] >= data[_canonical(right)] for left, right in SEPSIGLNC_PAIRS))


SEPSIG28_WEIGHTS = {
    "BOLA3-AS1": 0.254, "LINC00354": 0.1996, "C5ORF27": 0.1537, "RP1-187B23.1": -0.1427,
    "MBNL1-AS1": -0.1419, "LINC01420": -0.1140, "RP13-436F16.1": 0.1060,
    "CTB-31O20.2": 0.1023, "LINC01425": 0.0949, "C10ORF25": -0.0763,
    "RP11-111M22.3": 0.0743, "LAMTOR5-AS1": 0.0739, "FLJ37453": 0.0713,
    "AX746755": -0.0690, "TTTY12": 0.0678, "ASMTL-AS1": -0.0535,
    "LOC101928491": 0.0461, "RBM26-AS1": -0.0438, "ANP32A-IT1": 0.0437,
    "LOC101060691": 0.0319, "MSH5": -0.0311, "LOC100507221": 0.0289,
    "RP11-1137G4.3": -0.0245, "LOC100506457": 0.0237, "MIR612": -0.0189,
    "AC114730.11": 0.0079, "LOC101927526": 0.0026, "LINC01019": -0.0020,
}


def sepsig28_weighted_score(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """Published 28-lncRNA weighted sum; higher favors sepsis."""

    data = prepare_expression(expression)
    _required(data, tuple(SEPSIG28_WEIGHTS))
    return sum(weight * data[_canonical(gene)] for gene, weight in SEPSIG28_WEIGHTS.items())


LIN7_WEIGHTS = {
    "ADRB2": -0.4102, "CTSG": 0.1825, "CX3CR1": -0.1810, "CXCR6": 0.8549,
    "IL4R": -0.4270, "LTB": -0.5605, "TMSB10": -0.6836,
}


def lin7_mortality_score(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """Seven-immune-gene mortality risk score; higher indicates higher risk."""

    data = prepare_expression(expression)
    _required(data, tuple(LIN7_WEIGHTS))
    return sum(weight * data[_canonical(gene)] for gene, weight in LIN7_WEIGHTS.items())


SOM_M1 = ("NQO2", "SLPI", "ORM1", "KLHL2", "ANXA3", "TXN", "AQP9", "BCL6", "DOK3", "PFKFB4", "TYK2")
SOM_M2 = ("BCL2L11", "BCAT1", "BTBD7", "CEP55", "HMMR", "PRC1", "KIF15", "CAMP", "CEACAM8", "DEFA4", "LCN2", "CTSG", "AZU1")
SOM_M3 = ("MAFB", "OASL", "UBE2L6", "VAMP5", "CCL2", "NAPA", "ATG3", "VRK2", "TMEM123", "CASP7")
SOM_M4 = ("DOK2", "HLA-DPB1", "BUB3", "SMYD2", "SIDT1", "EXOC2", "TRIB2", "KLRB1")


def severe_or_mild_score(expression: Mapping[str, Any] | Iterable[tuple[str, Any]]) -> float:
    """[GM(M1)+GM(M2)]/[GM(M3)+GM(M4)]; higher indicates severe response."""

    data = prepare_expression(expression)
    numerator = _geometric_mean(_required(data, SOM_M1), "SoM module 1") + _geometric_mean(_required(data, SOM_M2), "SoM module 2")
    denominator = _geometric_mean(_required(data, SOM_M3), "SoM module 3") + _geometric_mean(_required(data, SOM_M4), "SoM module 4")
    if denominator == 0:
        raise SignatureInputError("SoM denominator is zero")
    return numerator / denominator


@dataclass(frozen=True)
class SignatureSpec:
    signature_id: str
    function_name: str
    output_direction: str
    grade: str


SPECS = (
    SignatureSpec("SIG001", "sepsis_metascore_dgm", "higher=infection/sepsis", "A2"),
    SignatureSpec("SIG002", "septicyte_lab_microarray", "higher=sepsis", "A2"),
    SignatureSpec("SIG003", "faim3_plac8_ratio", "higher=no-CAP/noninfected", "A2"),
    SignatureSpec("SIG004", "snip_score", "higher=abdominal sepsis", "A2"),
    SignatureSpec("SIG023", "herberg_two_gene_drs", "higher=bacterial", "A2"),
    SignatureSpec("SIG022", "sweeney_bv_ordering_score", "higher=viral", "A2"),
    SignatureSpec("SIG010", "sepsiglnc_pair_count", "higher=sepsis", "A2"),
    SignatureSpec("SIG032", "sepsig28_weighted_score", "higher=sepsis", "A2"),
    SignatureSpec("SIG033", "lin7_mortality_score", "higher=mortality risk", "A2"),
    SignatureSpec("SIG034", "severe_or_mild_score", "higher=severe host response", "A2"),
)


# Stage 3 cross-platform set. The function outputs preserve the primary-source
# orientation; SCORE_DIRECTION below is applied exactly once to obtain Q.
STAGE3_FUNCTIONS = {
    "SIG001": sepsis_metascore_dgm,
    "SIG002": septicyte_lab_microarray,
    "SIG003": faim3_plac8_ratio,
    "SIG004": snip_score,
    "SIG022": sweeney_bv_ordering_score,
    "SIG023": herberg_two_gene_drs,
    "SIG033": lin7_mortality_score,
    "SIG034": severe_or_mild_score,
}

STAGE3_REQUIRED_GENES = {
    "SIG001": SMS_UP + SMS_DOWN,
    "SIG002": ("PLAC8", "LAMP1", "PLA2G7", "CEACAM4"),
    "SIG003": ("FAIM3", "PLAC8"),
    "SIG004": ("NLRP1", "IDNK", "PLAC8"),
    "SIG022": SWEENEY_BV_VIRAL + SWEENEY_BV_BACTERIAL,
    "SIG023": ("FAM89A", "IFI44L"),
    "SIG033": tuple(LIN7_WEIGHTS),
    "SIG034": SOM_M1 + SOM_M2 + SOM_M3 + SOM_M4,
}

# d_j in Q=d_j*S. SIG003 is inverted because the original FAIM3:PLAC8 ratio
# increases toward noninfection; all other functions already increase toward
# their explicitly frozen target-positive state.
SCORE_DIRECTION = {
    "SIG001": 1,
    "SIG002": 1,
    "SIG003": -1,
    "SIG004": 1,
    "SIG022": 1,
    "SIG023": 1,
    "SIG033": 1,
    "SIG034": 1,
}

TARGET_POSITIVE_STATE = {
    "SIG001": "sepsis/infection",
    "SIG002": "sepsis",
    "SIG003": "sepsis/infection",
    "SIG004": "abdominal sepsis",
    "SIG022": "viral infection",
    "SIG023": "bacterial infection",
    "SIG033": "higher mortality risk",
    "SIG034": "severe host response",
}
