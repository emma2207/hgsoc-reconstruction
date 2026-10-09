import tempfile
import unittest
from pathlib import Path

import pandas as pd

from analysis_shared_functions import load_heatmap_data
from standardize_heatmap_data import export_heatmap_matrices


FIXTURE_DATA = Path(__file__).parent / "fixtures" / "standardize_heatmap"


class StandardizeHeatmapParsingTests(unittest.TestCase):
    def _export_tool(self, tool, *, read_depth="0"):
        output_dir = tempfile.TemporaryDirectory()
        self.addCleanup(output_dir.cleanup)
        paths = export_heatmap_matrices(
            data_path=FIXTURE_DATA,
            output_dir=output_dir.name,
            dataset="demo",
            tools=[tool],
            read_depth=read_depth,
            mod1="_1",
            mod2="_2",
            pseudobulk=True,
            ncells=1000,
        )
        self.assertEqual(len(paths), 1)
        return pd.read_csv(paths[0], index_col=0)

    def test_shared_heatmap_loader_dispatches_to_each_raw_parser(self):
        sample_1 = "GSM0000001_pseudobulk_n_barcodes_1000_1"
        sample_2 = "GSM0000002_pseudobulk_n_barcodes_1000_1"
        sample_1_2 = "GSM0000001_pseudobulk_n_barcodes_1000_2"
        sample_2_2 = "GSM0000002_pseudobulk_n_barcodes_1000_2"
        expected = {
            "Conpair": (sample_1, sample_1_2, 97.25),
            "CrosscheckFingerprints": (sample_1, sample_1_2, 12.5),
            "HYSYS": (sample_1, sample_1_2, 0.987),
            "NGSCheckmate": (sample_1, sample_1_2, 0.87),
            "ntsm": (sample_1, sample_1_2, 0.42),
            "OmicsPrint": (sample_1, sample_1_2, 1.75),
            "Peddy": (sample_2, sample_1_2, 0.125),
            "Somalier": (sample_1, sample_1_2, 1.0),
            "TimeAttackGenComp": (sample_1, sample_1_2, 0.01),
            "Vireo": (sample_1, sample_1_2, 0.01),
        }
        for tool, (row, column, expected_value) in expected.items():
            with self.subTest(tool=tool):
                matrix = load_heatmap_data(
                    str(FIXTURE_DATA),
                    True,
                    tool,
                    "demo",
                    1000,
                    "0",
                    "_1",
                    "_2",
                )
                self.assertIsInstance(matrix, pd.DataFrame)
                self.assertIn(row, matrix.index)
                self.assertIn(column, matrix.columns)
                self.assertAlmostEqual(
                    float(str(matrix.loc[row, column])), expected_value
                )

    def test_conpair_raw_pair_file_is_parsed_and_exported(self):
        with tempfile.TemporaryDirectory() as output_dir:
            paths = export_heatmap_matrices(
                data_path=FIXTURE_DATA,
                output_dir=output_dir,
                dataset="demo",
                tools=["Conpair"],
                read_depth="0",
                mod1="_1",
                mod2="_2",
                pseudobulk=True,
                ncells=1000,
            )

            self.assertEqual(len(paths), 1)
            result = pd.read_csv(paths[0], index_col=0)
            self.assertEqual(result.shape, (1, 1))
            self.assertEqual(
                result.index.tolist(),
                ["GSM0000001_pseudobulk_n_barcodes_1000_1"],
            )
            self.assertEqual(
                result.columns.tolist(),
                ["GSM0000001_pseudobulk_n_barcodes_1000_2"],
            )
            self.assertAlmostEqual(float(str(result.iat[0, 0])), 97.25)

    def test_hysys_raw_tsv_is_parsed_filtered_and_exported(self):
        with tempfile.TemporaryDirectory() as output_dir:
            paths = export_heatmap_matrices(
                data_path=FIXTURE_DATA,
                output_dir=output_dir,
                dataset="demo",
                tools=["HYSYS"],
                read_depth="0",
                mod1="_1",
                mod2="_2",
                pseudobulk=True,
                ncells=1000,
            )

            self.assertEqual(len(paths), 1)
            result = pd.read_csv(paths[0], index_col=0)
            self.assertEqual(result.shape, (2, 1))
            self.assertEqual(
                result.columns.tolist(),
                ["GSM0000001_pseudobulk_n_barcodes_1000_2"],
            )
            self.assertEqual(
                set(result.index),
                {
                    "GSM0000001_pseudobulk_n_barcodes_1000_1",
                    "GSM0000002_pseudobulk_n_barcodes_1000_1",
                },
            )
            self.assertAlmostEqual(
                float(
                    str(
                        result.loc[
                            "GSM0000001_pseudobulk_n_barcodes_1000_1",
                            "GSM0000001_pseudobulk_n_barcodes_1000_2",
                        ]
                    )
                ),
                0.987,
            )
            self.assertAlmostEqual(
                float(
                    str(
                        result.loc[
                            "GSM0000002_pseudobulk_n_barcodes_1000_1",
                            "GSM0000001_pseudobulk_n_barcodes_1000_2",
                        ]
                    )
                ),
                0.125,
            )

    def test_crosscheck_fingerprints_raw_metrics_are_parsed(self):
        result = self._export_tool("CrosscheckFingerprints")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            12.5,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            -3.25,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            9.0,
        )

    def test_ngscheckmate_full_precision_matrix_is_parsed(self):
        result = self._export_tool("NGSCheckmate")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.87,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            0.88,
        )

    def test_ntsm_pairwise_table_is_parsed_with_crosscheck_sample_universe(self):
        result = self._export_tool("ntsm")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.42,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            0.33,
        )
        self.assertTrue(pd.isna(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"]))

    def test_omicsprint_allele_sharing_table_is_parsed(self):
        result = self._export_tool("OmicsPrint")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            1.75,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.25,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            1.2,
        )
        self.assertTrue(pd.isna(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"]))

    def test_peddy_csv_is_parsed(self):
        result = self._export_tool("Peddy")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.0,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.125,
        )

    def test_somalier_pairs_tsv_is_parsed(self):
        result = self._export_tool("Somalier")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            1.0,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.12,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            1.0,
        )

    def test_timeattackgencomp_matrix_is_parsed(self):
        result = self._export_tool("TimeAttackGenComp")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.01,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            0.02,
        )

    def test_vireo_similarity_matrix_is_parsed_with_crosscheck_sample_universe(self):
        result = self._export_tool("Vireo")

        self.assertEqual(result.shape, (2, 2))
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000001_pseudobulk_n_barcodes_1000_1", "GSM0000001_pseudobulk_n_barcodes_1000_2"])),
            0.01,
        )
        self.assertAlmostEqual(
            float(str(result.loc["GSM0000002_pseudobulk_n_barcodes_1000_1", "GSM0000002_pseudobulk_n_barcodes_1000_2"])),
            0.02,
        )


if __name__ == "__main__":
    unittest.main()
