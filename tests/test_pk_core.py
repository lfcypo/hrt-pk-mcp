import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from pk_mcp.pk_params import (
    Compound, Hormone, Route, DoseEvent, PKParams,
    HORMONE_PARAMS, COMPOUND_INFO,
    compound_hormone,
)
from pk_mcp.pk_core import (
    resolve_params,
    _bateman_amount, _analytic_3c,
    inj_amount, one_comp_amount, dual_abs_amount,
    patch_amount, compute_concentration_at, simulate_timeline,
)


def approx(a: float, b: float, rel_tol: float = 1e-9, abs_tol: float = 1e-9) -> bool:
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


def test_bateman_zero_dose():
    assert _bateman_amount(dose_mg=0, F=1.0, ka=1.0, ke=0.5, t=1.0) == 0.0


def test_bateman_negative_time():
    assert _bateman_amount(dose_mg=10, F=1.0, ka=1.0, ke=0.5, t=-1.0) == 0.0


def test_bateman_zero_ka():
    assert _bateman_amount(dose_mg=10, F=1.0, ka=0, ke=0.5, t=1.0) == 0.0


def test_bateman_peak_after_dose():
    """Bateman should peak after t=0."""
    t_values = [0, 0.5, 1.0, 2.0, 4.0, 8.0]
    amounts = [_bateman_amount(10, 0.5, 1.0, 0.3, t) for t in t_values]
    assert amounts[0] == 0.0
    peak_idx = max(range(len(amounts)), key=lambda i: amounts[i])
    assert 0 < peak_idx < len(amounts) - 1


def test_bateman_repeated_rates():
    """ka == ke uses repeated-roots formula."""
    result = _bateman_amount(10, 0.5, 0.3, 0.3, 2.0)
    expected = 10 * 0.5 * 0.3 * 2.0 * math.exp(-0.3 * 2.0)
    assert approx(result, expected)


def test_bateman_near_repeated_rates():
    """Very close rates are numerically stable."""
    result = _bateman_amount(10, 0.5, 0.3, 0.3000000001, 2.0)
    assert result > 0
    expected_limit = 10 * 0.5 * 0.3 * 2.0 * math.exp(-0.3 * 2.0)
    assert approx(result, expected_limit, rel_tol=1e-4)


def test_bateman_decay_to_zero():
    assert _bateman_amount(10, 0.5, 1.0, 0.3, t=1000) < 1e-10


def test_3c_zero_dose():
    assert _analytic_3c(tau=1.0, dose_mg=0, F=1.0, k1=1.0, k2=0.5, k3=0.1) == 0.0


def test_3c_negative_tau():
    assert _analytic_3c(tau=-1.0, dose_mg=10, F=1.0, k1=1.0, k2=0.5, k3=0.1) == 0.0


def test_3c_all_rates_equal():
    result = _analytic_3c(2.0, 10, 0.8, 0.2, 0.2, 0.2)
    expected = 10 * 0.8 * 0.2 * 0.2 * 2.0 * 2.0 * 0.5 * math.exp(-0.2 * 2.0)
    assert approx(result, expected)


def test_3c_k2_eq_k3():
    r = _analytic_3c(2.0, 10, 0.8, 0.5, 0.2, 0.2)
    rp = _analytic_3c(2.0, 10, 0.8, 0.5, 0.2 + 1e-8, 0.2)
    assert r > 0 and approx(r, rp, rel_tol=1e-4)


def test_3c_k1_eq_k2():
    r = _analytic_3c(2.0, 10, 0.8, 0.3, 0.3, 0.1)
    rp = _analytic_3c(2.0, 10, 0.8, 0.3 + 1e-8, 0.3, 0.1)
    assert r > 0 and approx(r, rp, rel_tol=1e-4)


def test_3c_k1_eq_k3():
    assert _analytic_3c(2.0, 10, 0.8, 0.3, 0.5, 0.3) > 0


def test_3c_distinct_rates():
    dose, F = 10.0, 0.8
    k1, k2, k3 = 0.5, 0.3, 0.1
    t = 2.0
    result = _analytic_3c(t, dose, F, k1, k2, k3)
    term1 = math.exp(-k1 * t) / ((k1 - k2) * (k1 - k3))
    term2 = math.exp(-k2 * t) / ((k2 - k1) * (k2 - k3))
    term3 = math.exp(-k3 * t) / ((k3 - k1) * (k3 - k2))
    expected = dose * F * k1 * k2 * (term1 + term2 + term3)
    assert approx(result, expected, rel_tol=1e-12)


def test_3c_decay_to_zero():
    assert _analytic_3c(1000, 10, 0.8, 0.5, 0.3, 0.1) < 1e-30


def test_injection_amounts_with_known_params():
    """EV injection: amounts should be non-negative with a reasonable peak."""
    p = PKParams(0.0216, 0.0138, 0.07, 0.041, 0.0622582882,
                 0.4, 0.0622582882, 0.0622582882, 0, 0, 0)
    dose = 5.0

    # Near-zero at t=0 (allow FP noise)
    assert approx(inj_amount(0, dose, p), 0.0, abs_tol=1e-12)

    # Positive at 24h
    assert inj_amount(24, dose, p) > 0

    # Non-negative at 720h
    assert inj_amount(720, dose, p) >= -1e-12

    # Peak between 5-150h
    amounts = [inj_amount(t, dose, p) for t in range(1, 201, 5)]
    peak_idx = max(range(len(amounts)), key=lambda i: amounts[i])
    peak_hour = 1 + peak_idx * 5
    assert 5 < peak_hour < 150, f"Peak at {peak_hour}h"


def test_injection_superposition():
    """Two injections sum linearly."""
    p = PKParams(0.0216, 0.0138, 0.07, 0.041, 0.0622582882,
                 0.4, 0.0622582882, 0.0622582882, 0, 0, 0)
    s = inj_amount(48, 5, p) + inj_amount(48, 3, p)
    c = inj_amount(48, 8, p)
    assert approx(s, c, rel_tol=1e-12)


def test_oral_amount():
    p = PKParams(0.32, 0, 0, 0.41, 0.03, 1.0, 0.03, 0.03, 0, 0, 0)
    times = list(range(0, 24))
    amounts = [one_comp_amount(t, 2.0, p) for t in times]
    assert amounts[0] == 0.0
    peak_idx = max(range(len(amounts)), key=lambda i: amounts[i])
    assert peak_idx < 12, f"Oral E2 peak too late: {times[peak_idx]}h"


def test_gel_amount():
    """Gel (slow abs): peak around 7-8h."""
    p = PKParams(0.022, 0, 0, 0.41, 0.06, 1.0, 0.06, 0.06, 0, 0, 0)
    times = list(range(0, 72, 2))
    amounts = [one_comp_amount(t, 1.0, p) for t in times]
    peak_idx = max(range(len(amounts)), key=lambda i: amounts[i])
    peak_hour = times[peak_idx]
    # With k1=0.022, k3=0.41: theoretical peak at ~7.5h
    assert 6 <= peak_hour <= 12, f"Gel peak at {peak_hour}h (expected ~7-8h)"


def test_patch_zero_order():
    """Zero-order patch: accumulation during wear, decay after."""
    p = PKParams(0, 0, 0, 0.41, 1.0, 1.0, 1.0, 1.0, 0, 0, rate_mg_h=0.1667)
    wear_h = 24.0
    a12 = patch_amount(12, 4.0, wear_h, p)
    a24 = patch_amount(24, 4.0, wear_h, p)
    a36 = patch_amount(36, 4.0, wear_h, p)
    a200 = patch_amount(200, 4.0, wear_h, p)
    assert a12 > 0
    assert a24 > a12, "Amount increases during wear (zero-order)"
    assert a36 < a24, "Amount decreases after removal"
    assert a200 < a36


def test_patch_first_order():
    """First-order patch with slow abs vs fast elim: peaks early then declines.
    
    With k1=0.0075 (t½=92h) and k3=0.41 (t½=1.7h), absorption is so slow
    that elimination outpaces it after the early peak (~10h) even during wear.
    """
    p = PKParams(0.0075, 0, 0, 0.41, 1.0, 1.0, 1.0, 1.0, 0, 0, rate_mg_h=0)
    wear_h = 72.0
    a24 = patch_amount(24, 0.1, wear_h, p)
    a48 = patch_amount(48, 0.1, wear_h, p)
    assert a24 > 0
    # With slow abs vs fast elim, peak at ~10h then declines
    # 48h > 24h means DECREASING, not accumulating
    # This is correct: the first-order patch model is Bateman with
    # k1 << k3, so elimination dominates after the early peak.
    assert a48 < a24, "First-order patch declines after early peak (flip-flop kinetics)"


def test_dual_abs():
    p = PKParams(0.450550912583251, 0.0142806935343998, 0,
                 0.44024417908217306, 0.025919316335729803,
                 1.0, 0.025919316335729803, 0.0,
                 2.75, 0.0, 0)
    assert dual_abs_amount(5.0, 225.0, p) > 0
    assert dual_abs_amount(1.0, 225.0, p) == 0.0  # before lag


def test_estradiol_concentration_scale():
    from pk_mcp.pk_params import concentration_scale, vd_ml, HORMONE_PARAMS
    s = concentration_scale(Hormone.estradiol, HORMONE_PARAMS[Hormone.estradiol].concentration_unit)
    v = vd_ml(70.0, Hormone.estradiol)
    assert approx(1.0 * s / v, 1.0 * 1e9 / 140000.0)


def test_testosterone_concentration_scale():
    from pk_mcp.pk_params import concentration_scale, vd_ml, HORMONE_PARAMS
    s = concentration_scale(Hormone.testosterone, HORMONE_PARAMS[Hormone.testosterone].concentration_unit)
    v = vd_ml(70.0, Hormone.testosterone)
    assert approx(1.0 * s / v, 1.0 * 1e8 / 140000.0)


def test_single_dose_concentration():
    e = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    c = compute_concentration_at(hour=2.0, events=[e], hormone=Hormone.estradiol)
    assert 0 < c < 500, f"Unreasonable: {c}"
    cf = compute_concentration_at(hour=168.0, events=[e], hormone=Hormone.estradiol)
    assert cf < 1.0, f"7-day should be near zero: {cf}"


def test_multiple_dose_superposition():
    e1 = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    e2 = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    c1 = compute_concentration_at(4.0, [e1], Hormone.estradiol)
    c2 = compute_concentration_at(4.0, [DoseEvent(Compound.E2, Route.oral, 0, 4.0)], Hormone.estradiol)
    c3 = compute_concentration_at(4.0, [e1, e2], Hormone.estradiol)
    assert approx(c2, c1 * 2, rel_tol=1e-10)
    assert approx(c3, c2, rel_tol=1e-10)


def test_steady_state_repeat_dosing():
    es = [DoseEvent(Compound.E2, Route.oral, i * 12, 2.0) for i in range(14)]
    ct = compute_concentration_at(12 * 6, es, Hormone.estradiol)
    cp = compute_concentration_at(12 * 6 + 1, es, Hormone.estradiol)
    assert ct > 0
    assert cp > ct
    cfp = compute_concentration_at(1, es[:1], Hormone.estradiol)
    assert cp > cfp, "SS peak > first-dose peak"


def test_different_body_weights():
    e = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    c70 = compute_concentration_at(2.0, [e], Hormone.estradiol, 70.0)
    c100 = compute_concentration_at(2.0, [e], Hormone.estradiol, 100.0)
    assert approx(c100 / c70, 70.0 / 100.0, rel_tol=1e-10)


def test_no_hormone_events_zero():
    assert compute_concentration_at(10.0, [], Hormone.estradiol) == 0.0


def test_wrong_hormone_events():
    e = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    ec = compute_concentration_at(2.0, [e], Hormone.estradiol)
    tc = compute_concentration_at(2.0, [e], Hormone.testosterone)
    assert ec > 0
    assert tc == 0.0


def test_resolve_params_oral_e2():
    p = resolve_params(DoseEvent(Compound.E2, Route.oral, 0, 2.0))
    assert approx(p.k1_fast, 0.32)
    assert approx(p.k3, 0.41)
    assert approx(p.F, 0.03)
    assert approx(p.frac_fast, 1.0)


def test_resolve_params_injection_ev():
    p = resolve_params(DoseEvent(Compound.EV, Route.injection, 0, 5.0))
    assert approx(p.k3, 0.041)
    assert approx(p.frac_fast, 0.4)
    assert approx(p.k2, 0.07)
    assert approx(p.F, 0.062258288229969413)


def test_resolve_params_patch_zero_order():
    e = DoseEvent(Compound.E2, Route.patch_apply, 0, 0.1, release_rate_ug_per_day=100)
    p = resolve_params(e)
    assert p.rate_mg_h > 0
    assert approx(p.rate_mg_h, 100.0 / 24000.0)
    assert approx(p.k3, 0.41)


def test_resolve_params_gel():
    p = resolve_params(DoseEvent(Compound.E2, Route.gel, 0, 1.0))
    assert approx(p.k1_fast, 0.022)
    assert approx(p.F, 0.06)


def test_resolve_params_sublingual_default():
    p = resolve_params(DoseEvent(Compound.E2, Route.sublingual, 0, 2.0))
    assert approx(p.frac_fast, 0.11)
    assert approx(p.k1_fast, 1.8)
    assert approx(p.F_fast, 1.0)
    assert approx(p.F_slow, 0.03)


def test_resolve_params_sublingual_custom_theta():
    e = DoseEvent(Compound.E2, Route.sublingual, 0, 2.0, sublingual_theta=0.3)
    p = resolve_params(e)
    assert approx(p.frac_fast, 0.3)


def test_compound_info():
    info = COMPOUND_INFO[Compound.E2]
    assert info.full_name == "Estradiol"
    assert approx(info.molecular_weight, 272.38)
    assert info.to_active_factor == 1.0
    assert not info.is_prodrug

    info = COMPOUND_INFO[Compound.EV]
    assert info.full_name == "Estradiol Valerate"
    assert info.is_prodrug
    assert approx(info.to_active_factor, 272.38 / 356.50)

    info = COMPOUND_INFO[Compound.T]
    assert info.full_name == "Testosterone"
    assert approx(info.molecular_weight, 288.42)

    info = COMPOUND_INFO[Compound.TU]
    assert info.full_name == "Testosterone Undecanoate"
    assert approx(info.molecular_weight, 456.70)
    assert info.is_prodrug


def test_compound_hormone_mapping():
    assert compound_hormone(Compound.E2) == Hormone.estradiol
    assert compound_hormone(Compound.EV) == Hormone.estradiol
    assert compound_hormone(Compound.T) == Hormone.testosterone
    assert compound_hormone(Compound.TC) == Hormone.testosterone
    assert compound_hormone(Compound.TU) == Hormone.testosterone


def test_hormone_params():
    e2 = HORMONE_PARAMS[Hormone.estradiol]
    assert e2.vd_per_kg == 2.0
    assert e2.k_clear == 0.41
    assert e2.k_clear_injection == 0.041

    t = HORMONE_PARAMS[Hormone.testosterone]
    assert t.vd_per_kg == 2.0
    assert t.k_clear == 0.6
    assert t.k_clear_injection == 0.03


def test_simulate_timeline_empty():
    r = simulate_timeline([], Hormone.estradiol)
    assert len(r.time_h) == 0


def test_simulate_timeline_single():
    e = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    r = simulate_timeline([e], Hormone.estradiol, number_of_steps=100)
    assert len(r.time_h) == 100
    assert all(c >= -1e-12 for c in r.concentrations)
    assert r.auc > 0


def test_simulate_timeline_multiple():
    es = [
        DoseEvent(Compound.E2, Route.oral, 0, 2.0),
        DoseEvent(Compound.EV, Route.injection, 0, 5.0),
    ]
    r = simulate_timeline(es, Hormone.estradiol, number_of_steps=200)
    assert len(r.time_h) == 200
    assert r.auc > 0


def test_simulate_timeline_interpolation():
    e = DoseEvent(Compound.E2, Route.oral, 0, 2.0)
    r = simulate_timeline([e], Hormone.estradiol, number_of_steps=1000)
    c = r.concentration_at(2.0)
    direct = compute_concentration_at(2.0, [e], Hormone.estradiol)
    assert c is not None
    assert approx(c, direct, rel_tol=1e-2)


if __name__ == "__main__":
    test_fns = [n for n in dir() if n.startswith("test_") and callable(globals()[n])]
    passed = failed = 0
    for name in sorted(test_fns):
        fn = globals()[name]
        try:
            fn()
            passed += 1
        except AssertionError as e:
            failed += 1
            print(f"FAIL: {name} — {e}")
        except Exception as e:
            failed += 1
            import traceback

            print(f"ERROR: {name} — {e}")
            traceback.print_exc()
    print(f"\n{'=' * 50}")
    print(f"Results: {passed}/{passed + failed} passed"
          + (f", {failed} FAILED" if failed else ", ALL PASSED"))
    sys.exit(1 if failed else 0)
