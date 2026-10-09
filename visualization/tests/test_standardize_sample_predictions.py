import tempfile
import unittest
from pathlib import Path

import pandas as pd

from analysis_shared_functions import load_sample_matching_results
from standardize_heatmap_data import export_heatmap_matrices
from standardize_sample_predictions import PREDICTION_THRESHOLDS, export_sample_predictions


FIXTURE_DATA = Path(__file__).parent / "fixtures" / "standardize_heatmap"
SAMPLE_1 = "GSM0000001_pseudobulk_n_barcodes_1000_1"
SAMPLE_2 = "GSM0000002_pseudobulk_n_barcodes_1000_1"
SAMPLE_1_2 = "GSM0000001_pseudobulk_n_barcodes_1000_2"
SAMPLE_2_2 = "GSM0000002_pseudobulk_n_barcodes_1000_2"


class StandardizeSamplePredictionParsingTests(unittest.TestCase):
    def test_shared_prediction_loader_dispatches_to_each_prediction_parser(self):
        tools = [
            "Conpair",
            "CrosscheckFingerprints",
            "HYSYS",
            "NGSCheckmate",
            "ntsm",
            "OmicsPrint",
            "Peddy",
            "Somalier",
            "TimeAttackGenComp",
            "Vireo",
        ]
        expected_row = SAMPLE_1
        expected_column = SAMPLE_1_2

        with tempfile.TemporaryDirectory() as heatmap_dir:
            heatmap_paths = export_heatmap_matrices(
                data_path=FIXTURE_DATA,
                output_dir=heatmap_dir,
                dataset="demo",
                tools=tools,
                read_depth="0",
                mod1="_1",
                mod2="_2",
                pseudobulk=True,
                ncells=1000,
            )
            self.assertEqual(len(heatmap_paths), len(tools))

            for tool in tools:
                with self.subTest(tool=tool):
                    matrix = load_sample_matching_results(
                        str(FIXTURE_DATA),
                        True,
                        tool,
                        "demo",
                        1000,
                        "0",
                        mod1="_1",
                        mod2="_2",
                        standardized_heatmaps_path=heatmap_dir,
                    )
                    self.assertIsInstance(matrix, pd.DataFrame)
                    self.assertIn(expected_row, matrix.index)
                    self.assertIn(expected_column, matrix.columns)
                    self.assertEqual(int(matrix.loc[expected_row, expected_column]), 1)

    def test_raw_results_and_standardized_scores_export_binary_calls(self):
        tools = [
            "Conpair",
            "CrosscheckFingerprints",
            "HYSYS",
            "NGSCheckmate",
            "ntsm",
            "OmicsPrint",
            "Peddy",
            "Somalier",
            "TimeAttackGenComp",
            "Vireo",
        ]
        with tempfile.TemporaryDirectory() as heatmap_dir, tempfile.TemporaryDirectory() as prediction_dir:
            heatmap_paths = export_heatmap_matrices(
                data_path=FIXTURE_DATA,
                output_dir=heatmap_dir,
                dataset="demo",
                tools=tools,
                read_depth="0",
                mod1="_1",
                mod2="_2",
                pseudobulk=True,
                ncells=1000,
            )
            self.assertEqual(len(heatmap_paths), len(tools))

            prediction_paths = export_sample_predictions(
                data_path=FIXTURE_DATA,
                output_dir=prediction_dir,
                dataset="demo",
                tools=tools,
                read_depth="0",
                mod1="_1",
                mod2="_2",
                pseudobulk=True,
                ncells=1000,
                standardized_heatmaps_path=heatmap_dir,
            )
            self.assertEqual(len(prediction_paths), len(tools))

            results = {}
            for path in prediction_paths:
                tool = path.name.removesuffix(".csv").split(
                    "_threshold_", maxsplit=1
                )[0]
                matrix = pd.read_csv(path, index_col=0)
                values = {
                    value
                    for value in matrix.to_numpy().ravel()
                    if pd.notna(value)
                }
                self.assertLessEqual(values, {0, 1}, msg=f"unexpected value in {path.name}")
                results[tool] = matrix

        self.assertEqual(results["Conpair"].loc[SAMPLE_1, SAMPLE_1_2], 1)

        crosscheck = results["CrosscheckFingerprints"]
        self.assertEqual(crosscheck.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(crosscheck.loc[SAMPLE_2, SAMPLE_1_2], 0)
        self.assertEqual(crosscheck.loc[SAMPLE_2, SAMPLE_2_2], 1)

        hysys = results["HYSYS"]
        self.assertEqual(hysys.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertTrue(pd.isna(hysys.loc[SAMPLE_2, SAMPLE_1_2]))

        ngscheckmate = results["NGSCheckmate"]
        self.assertEqual(ngscheckmate.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(ngscheckmate.loc[SAMPLE_2, SAMPLE_1_2], 0)
        self.assertEqual(ngscheckmate.loc[SAMPLE_2, SAMPLE_2_2], 1)

        ntsm = results["ntsm"]
        self.assertEqual(ntsm.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(ntsm.loc[SAMPLE_2, SAMPLE_2_2], 1)
        self.assertEqual(ntsm.loc[SAMPLE_2, SAMPLE_1_2], 0)

        omicsprint = results["OmicsPrint"]
        self.assertEqual(omicsprint.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(omicsprint.loc[SAMPLE_2, SAMPLE_1_2], 0)
        self.assertEqual(omicsprint.loc[SAMPLE_2, SAMPLE_2_2], 0)
        self.assertTrue(pd.isna(omicsprint.loc[SAMPLE_1, SAMPLE_2_2]))

        peddy = results["Peddy"]
        self.assertEqual(peddy.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(peddy.loc[SAMPLE_2, SAMPLE_1_2], 0)
        self.assertEqual(peddy.loc[SAMPLE_2, SAMPLE_2_2], 1)

        somalier = results["Somalier"]
        self.assertEqual(somalier.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(somalier.loc[SAMPLE_2, SAMPLE_1_2], 0)
        self.assertEqual(somalier.loc[SAMPLE_2, SAMPLE_2_2], 1)

        timeattack = results["TimeAttackGenComp"]
        self.assertEqual(timeattack.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(timeattack.loc[SAMPLE_1, SAMPLE_2_2], 0)
        self.assertEqual(timeattack.loc[SAMPLE_2, SAMPLE_2_2], 1)

        vireo = results["Vireo"]
        self.assertEqual(vireo.loc[SAMPLE_1, SAMPLE_1_2], 1)
        self.assertEqual(vireo.loc[SAMPLE_2, SAMPLE_2_2], 1)
        self.assertEqual(vireo.loc[SAMPLE_1, SAMPLE_2_2], 0)


if __name__ == "__main__":
    unittest.main()
