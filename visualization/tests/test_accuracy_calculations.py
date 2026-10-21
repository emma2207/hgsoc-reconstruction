import unittest

import pandas as pd

from accuracy_calculations import calculate_accuracy_metrics_pseudobulk
from analysis_shared_functions import true_matches_pseudobulk


class CalculateAccuracyMetricsPseudobulkTests(unittest.TestCase):
    def setUp(self):
        self.index = ["patientA_1", "patientB_1"]
        self.columns = ["patientA_2", "patientB_2"]

    def assertMetricsAlmostEqual(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(float(str(actual_value)), expected_value)

    def test_true_matches_uses_sample_prefix_and_preserves_axes(self):
        rows = [
            "GSM7305260_pseudobulk_n_barcodes_10000_1",
            "GSM7305261_pseudobulk_n_barcodes_10000_1",
        ]
        columns = [
            "GSM7305260_pseudobulk_n_barcodes_10000_2",
            "GSM7305261_pseudobulk_n_barcodes_10000_2",
        ]
        matrix = pd.DataFrame([[0, 0], [0, 0]], index=rows, columns=columns)

        true_matches = true_matches_pseudobulk(matrix)

        self.assertEqual(true_matches.shape, matrix.shape)
        self.assertEqual(true_matches.index.tolist(), rows)
        self.assertEqual(true_matches.columns.tolist(), columns)
        self.assertEqual(
            true_matches.to_numpy().tolist(),
            [[1.0, 0.0], [0.0, 1.0]],
        )

    def test_perfect_matches_and_nonmatches(self):
        inferred = pd.DataFrame(
            [[1, 0], [0, 1]], index=self.index, columns=self.columns
        )

        metrics = calculate_accuracy_metrics_pseudobulk(inferred)

        self.assertMetricsAlmostEqual(metrics, (0.0, 1.0, 1.0, 1.0, 1.0, 1.0))

    def test_nan_is_counted_inconclusive_then_treated_as_nonmatch(self):
        inferred = pd.DataFrame(
            [[1, float("nan")], [0, 1]], index=self.index, columns=self.columns
        )

        metrics = calculate_accuracy_metrics_pseudobulk(inferred)

        self.assertMetricsAlmostEqual(metrics, (0.25, 1.0, 1.0, 1.0, 1.0, 1.0))

    def test_false_positive_changes_precision_and_accuracy(self):
        inferred = pd.DataFrame(
            [[1, 1], [0, 1]], index=self.index, columns=self.columns
        )

        metrics = calculate_accuracy_metrics_pseudobulk(inferred)

        self.assertMetricsAlmostEqual(
            metrics,
            (0.0, 0.75, 0.75, 2 / 3, 1.0, 0.8),
        )


if __name__ == "__main__":
    unittest.main()
