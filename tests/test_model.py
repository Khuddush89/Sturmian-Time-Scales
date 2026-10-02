"""Check analytical limits, domain conventions, and scientific failure cases."""

import numpy as np
import pytest

from sturmian.model import (
    count_zeros,
    dense_matrix,
    first_pulse_root,
    gap_matrix,
    gauge,
    propagate,
    pulse_profile,
    pulse_transfer,
    spectrum,
)


@pytest.mark.parametrize("alpha", [0.6, 0.8, 1.0])
def test_uniform_grid_matches_closed_form(alpha):
    n, h = 12, 0.5
    t = np.arange(n + 1) * h
    k = np.arange(1, n)
    expected = 4 * alpha * (alpha - h * (1 - alpha)) / h**2 * np.sin(k * np.pi / (2 * n)) ** 2
    np.testing.assert_allclose(spectrum(t, alpha), expected, rtol=1e-12, atol=1e-12)


def test_positive_gauge_preserves_signs_and_rejects_nonadmissible_order():
    t = np.array([0.0, 0.25, 0.8, 1.0])
    _, e = gauge(t, 0.75, 0.25)
    y = np.array([0.0, 1.0, -2.0, 0.0])
    np.testing.assert_array_equal(np.sign(e * y), np.sign(y))
    with pytest.raises(ValueError, match="positive-gauge"):
        gauge([0, 1], 0.5, 0.5)


@pytest.mark.parametrize("grid", [[0, 0, 1], [1, 0], [0, float("nan")], [0]])
def test_invalid_time_scales_raise_clear_errors(grid):
    with pytest.raises(ValueError):
        gauge(grid, 0.75, 0.25)


def test_spectral_mass_must_be_positive():
    with pytest.raises(ValueError, match="positive"):
        spectrum([0, 0.5, 1], 0.75, r=0)


def test_zero_counts_exclude_initial_zero_and_include_terminal_zero():
    crossings, actual = count_zeros([0, 1, 2], 1.0, [2.0, 0.0])
    assert crossings == []
    assert actual == [2]


def test_crossing_is_at_left_endpoint_and_requires_contained_gap():
    full, actual = count_zeros([0, 1, 2], 1.0, [3.0, 0.0])
    assert full == [1]
    assert actual == []
    assert count_zeros([0, 1], 1.0, [3.0]) == ([], [])


def test_tiny_amplitudes_do_not_hide_a_strict_sign_change():
    _, _, crossings, _ = propagate([0, 1], 1.0, initial=(1e-200, -2e-200))
    assert crossings == [0]


@pytest.mark.parametrize("s", [-0.5, 0.0, 0.45])
def test_dense_transfer_composes_and_preserves_wronskian(s):
    left, right = dense_matrix(s, 0.8, 0.3), dense_matrix(s, 0.8, 0.7)
    np.testing.assert_allclose(right @ left, dense_matrix(s, 0.8), atol=1e-13)
    assert np.linalg.det(gap_matrix(s, 0.8)) == pytest.approx(1.0, abs=1e-13)


def test_principal_pulse_root_satisfies_boundary_and_is_positive_inside():
    root = first_pulse_root(0.85)
    assert abs(pulse_transfer(root, 0.85)[0, 1]) < 1e-12
    assert all(np.min(y[1:-1]) > 0 for _, y, _ in pulse_profile(0.85, root))
