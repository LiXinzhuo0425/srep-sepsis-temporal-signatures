from __future__ import annotations

import math
import unittest

from signatures import (
    DuplicateGeneError,
    MissingGeneError,
    NonNumericExpressionError,
    SignatureInputError,
    faim3_plac8_ratio,
    herberg_two_gene_drs,
    prepare_expression,
    sepsis_metascore_dgm,
    septicyte_lab_microarray,
    snip_score,
    sepsig28_weighted_score,
    sepsiglnc_pair_count,
    lin7_mortality_score,
    severe_or_mild_score,
    sweeney_bv_ordering_score,
    SCORE_DIRECTION,
    STAGE3_FUNCTIONS,
    STAGE3_REQUIRED_GENES,
)


class TestInputGuardrails(unittest.TestCase):
    def test_gene_order_does_not_change_score(self):
        a = {"PLAC8": 8, "LAMP1": 6, "PLA2G7": 2, "CEACAM4": 1}
        b = {"CEACAM4": 1, "PLA2G7": 2, "LAMP1": 6, "PLAC8": 8}
        self.assertEqual(septicyte_lab_microarray(a), septicyte_lab_microarray(b))

    def test_missing_gene_raises(self):
        with self.assertRaises(MissingGeneError):
            septicyte_lab_microarray({"PLAC8": 1, "LAMP1": 2, "PLA2G7": 3})

    def test_duplicate_gene_raises(self):
        with self.assertRaises(DuplicateGeneError):
            prepare_expression([("PLAC8", 1), ("PLAC8", 2)])

    def test_historical_aliases_map(self):
        data = {
            "CEACAM1": 4, "ZDHHC19": 4, "C9orf95": 4, "GNA15": 4, "BATF": 4, "C3AR1": 4,
            "KIAA1370": 1, "TGFBI": 1, "MTCH1": 1, "RPGRIP1": 1, "HLA-DPB1": 1,
        }
        self.assertAlmostEqual(sepsis_metascore_dgm(data), 3.0)

    def test_multi_probe_mean_is_explicit(self):
        score = septicyte_lab_microarray({"PLAC8": [7, 9], "LAMP1": 6, "PLA2G7": 2, "CEACAM4": 1})
        self.assertEqual(score, 11.0)

    def test_nonnumeric_raises(self):
        with self.assertRaises(NonNumericExpressionError):
            herberg_two_gene_drs({"FAM89A": "bad", "IFI44L": 2})

    def test_nan_raises(self):
        with self.assertRaises(NonNumericExpressionError):
            herberg_two_gene_drs({"FAM89A": math.nan, "IFI44L": 2})


class TestManualArithmetic(unittest.TestCase):
    def test_stage3_set_is_frozen(self):
        self.assertEqual(
            tuple(STAGE3_FUNCTIONS),
            ("SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"),
        )
        self.assertEqual(set(STAGE3_FUNCTIONS), set(STAGE3_REQUIRED_GENES))
        self.assertEqual(SCORE_DIRECTION["SIG003"], -1)
        self.assertTrue(all(value in (-1, 1) for value in SCORE_DIRECTION.values()))

    def test_sepsis_metascore(self):
        data = {g: 9 for g in ("CEACAM1", "ZDHHC19", "NMRK1", "GNA15", "BATF", "C3AR1")}
        data.update({g: 4 for g in ("FAM214A", "TGFBI", "MTCH1", "RPGRIP1", "HLA-DPB1")})
        self.assertAlmostEqual(sepsis_metascore_dgm(data), 5.0)

    def test_septicyte(self):
        self.assertEqual(septicyte_lab_microarray({"PLAC8": 8, "LAMP1": 6, "PLA2G7": 2, "CEACAM4": 1}), 11)

    def test_faim3_plac8(self):
        self.assertEqual(faim3_plac8_ratio({"FAIM3": 6, "PLAC8": 3}), 2)

    def test_snip(self):
        self.assertEqual(snip_score({"NLRP1": 9, "IDNK": 3, "PLAC8": 2}), 3)

    def test_herberg(self):
        self.assertEqual(herberg_two_gene_drs({"FAM89A": 9, "IFI44L": 2}), 7)

    def test_sweeney_bv(self):
        data = {"IFI27": 9, "JUP": 9, "LAX1": 9, "HK3": 4, "TNIP1": 4, "GPAA1": 4, "CTSB": 4}
        self.assertAlmostEqual(sweeney_bv_ordering_score(data), 5)

    def test_direction_checks(self):
        self.assertGreater(herberg_two_gene_drs({"FAM89A": 10, "IFI44L": 1}), 0)
        self.assertGreater(sweeney_bv_ordering_score({"IFI27": 10, "JUP": 10, "LAX1": 10, "HK3": 1, "TNIP1": 1, "GPAA1": 1, "CTSB": 1}), 0)

    def test_extreme_values_remain_finite(self):
        score = septicyte_lab_microarray({"PLAC8": 1e12, "LAMP1": 1e12, "PLA2G7": -1e12, "CEACAM4": -1e12})
        self.assertTrue(math.isfinite(score))

    def test_zero_denominator_raises(self):
        with self.assertRaises(SignatureInputError):
            snip_score({"NLRP1": 1, "IDNK": 2, "PLAC8": 0})

    def test_sepsiglnc_ties_count(self):
        genes = {
            "ECRP", "CTD-2012K14.6", "LOC101926943", "MCM3AP-AS1", "MGC27345", "STARD7-AS1",
            "C5orf66", "AC090627.1", "GUSBP11", "LOC284112", "SDCBP2-AS1", "LOC101927974",
            "ITFG2", "LOC100506990", "MIRLET7DHG", "AC008753.4", "LOC101928817",
            "RP11-533E19.7", "AC100830.4",
        }
        self.assertEqual(sepsiglnc_pair_count({g: 5 for g in genes}), 14)

    def test_sepsig28_manual_weight(self):
        from signatures import SEPSIG28_WEIGHTS
        data = {g: 0 for g in SEPSIG28_WEIGHTS}
        data["BOLA3-AS1"] = 2
        self.assertAlmostEqual(sepsig28_weighted_score(data), 0.508)

    def test_lin7_manual(self):
        from signatures import LIN7_WEIGHTS
        data = {g: 1 for g in LIN7_WEIGHTS}
        self.assertAlmostEqual(lin7_mortality_score(data), sum(LIN7_WEIGHTS.values()))

    def test_severe_or_mild_manual(self):
        from signatures import SOM_M1, SOM_M2, SOM_M3, SOM_M4
        data = {g: 4 for g in SOM_M1 + SOM_M2}
        data.update({g: 2 for g in SOM_M3 + SOM_M4})
        self.assertAlmostEqual(severe_or_mild_score(data), 2.0)


if __name__ == "__main__":
    unittest.main()
