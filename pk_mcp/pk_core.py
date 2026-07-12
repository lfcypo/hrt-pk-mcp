from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional

from .pk_params import (
    Compound, Hormone, Route,
    DoseEvent, PKParams,
    HORMONE_PARAMS, TWO_PART_DEPOT, FORMATION_FRACTION,
    HYDROLYSIS_K2, ORAL_KABS, ORAL_BIOAVAILABILITY,
    ORAL_DUAL_ABSORPTION, KABS_SL,
    compound_hormone, concentration_scale, vd_ml,
)


def resolve_params(event: DoseEvent) -> PKParams:
    hormone = compound_hormone(event.compound)
    core = HORMONE_PARAMS[hormone]
    k3 = core.k_clear_injection if event.route == Route.injection else core.k_clear

    # 注射给药
    if event.route == Route.injection:
        depot = TWO_PART_DEPOT.get(event.compound)
        if depot is None:
            return PKParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        k1corr = core.depot_k1_corr
        k1_fast = depot.k1_fast * k1corr
        k1_slow = depot.k1_slow * k1corr
        form = FORMATION_FRACTION.get(event.compound, 1.0)
        return PKParams(
            k1_fast=k1_fast, k1_slow=k1_slow,
            k2=HYDROLYSIS_K2.get(event.compound, 0),
            k3=k3,
            F=form,
            frac_fast=depot.frac_fast,
            F_fast=form, F_slow=form,
            lag_fast_h=0, lag_slow_h=0,
            rate_mg_h=0,
        )

    # 贴片
    if event.route == Route.patch_apply:
        if event.release_rate_ug_per_day is not None and event.release_rate_ug_per_day > 0:
            rate_mg_h = event.release_rate_ug_per_day / 24_000.0
            return PKParams(
                k1_fast=0, k1_slow=0, k2=0, k3=k3,
                F=core.patch_release_scale,
                frac_fast=1.0,
                F_fast=core.patch_release_scale, F_slow=core.patch_release_scale,
                lag_fast_h=0, lag_slow_h=0,
                rate_mg_h=rate_mg_h,
            )
        else:
            return PKParams(
                k1_fast=core.patch_fallback_k1, k1_slow=0, k2=0, k3=k3,
                F=1.0,
                frac_fast=1.0,
                F_fast=1.0, F_slow=1.0,
                lag_fast_h=0, lag_slow_h=0,
                rate_mg_h=0,
            )

    # 凝胶
    if event.route == Route.gel:
        return PKParams(
            k1_fast=core.gel_k1, k1_slow=0, k2=0, k3=k3,
            F=core.gel_f_max,
            frac_fast=1.0,
            F_fast=core.gel_f_max, F_slow=core.gel_f_max,
            lag_fast_h=0, lag_slow_h=0,
            rate_mg_h=0,
        )

    # 口服
    if event.route == Route.oral:
        dual = ORAL_DUAL_ABSORPTION.get(event.compound)
        if dual is not None:
            return PKParams(
                k1_fast=dual.k_abs_fast, k1_slow=dual.k_abs_slow,
                k2=0, k3=dual.k_clear,
                F=dual.bioavailability_fast,
                frac_fast=dual.frac_fast,
                F_fast=dual.bioavailability_fast, F_slow=dual.bioavailability_slow,
                lag_fast_h=dual.lag_hours_fast, lag_slow_h=dual.lag_hours_slow,
                rate_mg_h=0,
            )
        k1 = ORAL_KABS.get(event.compound, 0)
        k2 = HYDROLYSIS_K2.get(event.compound, 0) if event.compound == Compound.EV else 0
        F = ORAL_BIOAVAILABILITY.get(event.compound, 0)
        return PKParams(
            k1_fast=k1, k1_slow=0, k2=k2, k3=k3,
            F=F,
            frac_fast=1.0,
            F_fast=F, F_slow=F,
            lag_fast_h=0, lag_slow_h=0,
            rate_mg_h=0,
        )

    # 舌下含服
    if event.route == Route.sublingual:
        if hormone != Hormone.estradiol:
            return PKParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        theta = event.sublingual_theta if event.sublingual_theta is not None else 0.11
        theta = max(0.0, min(1.0, theta))
        k1_fast = KABS_SL
        k1_slow = ORAL_KABS.get(event.compound, 0)
        F_fast = 1.0
        F_slow = ORAL_BIOAVAILABILITY.get(event.compound, 0)
        k2_val = HYDROLYSIS_K2.get(event.compound, 0) if event.compound == Compound.EV else 0
        return PKParams(
            k1_fast=k1_fast, k1_slow=k1_slow,
            k2=k2_val, k3=k3,
            F=1.0,
            frac_fast=theta,
            F_fast=F_fast, F_slow=F_slow,
            lag_fast_h=0, lag_slow_h=0,
            rate_mg_h=0,
        )

    # 贴片移除
    if event.route == Route.patch_remove:
        return PKParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

    return PKParams(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


_TOLERANCE = 1e-9


def _bateman_amount(dose_mg: float, F: float, ka: float, ke: float, t: float) -> float:
    """
    单室 Bateman
    """
    if t < 0 or dose_mg <= 0 or ka <= 0:
        return 0.0
    if abs(ka - ke) < _TOLERANCE:
        return dose_mg * F * ka * t * math.exp(-ke * t)
    return dose_mg * F * ka / (ka - ke) * (math.exp(-ke * t) - math.exp(-ka * t))


def _analytic_3c(tau: float, dose_mg: float, F: float,
                 k1: float, k2: float, k3: float) -> float:
    """
    三室链 A --k1--> B --k2--> C --k3--> 消除

    返回 tau 时 C 室的药量
    """
    if tau < 0 or k1 <= 0 or k2 <= 0 or k3 <= 0 or dose_mg <= 0 or F <= 0:
        return 0.0

    scaled_dose = dose_mg * F

    k1_eq_k2 = abs(k1 - k2) < _TOLERANCE
    k1_eq_k3 = abs(k1 - k3) < _TOLERANCE
    k2_eq_k3 = abs(k2 - k3) < _TOLERANCE

    if k1_eq_k2 and k1_eq_k3:
        return scaled_dose * k1 * k1 * tau * tau * 0.5 * math.exp(-k1 * tau)

    if k2_eq_k3:
        delta = k1 - k2
        if abs(delta) < _TOLERANCE:
            return 0.0
        n = scaled_dose * k1 * k2
        return n * (math.exp(-k1 * tau) + math.exp(-k2 * tau) * (delta * tau - 1)) / (delta * delta)

    if k1_eq_k2:
        delta = k1 - k3
        if abs(delta) < _TOLERANCE:
            return 0.0
        n = scaled_dose * k1 * k1
        return n * (math.exp(-k3 * tau) - math.exp(-k1 * tau) * (1 + delta * tau)) / (delta * delta)

    if k1_eq_k3:
        delta = k2 - k1
        if abs(delta) < _TOLERANCE:
            return 0.0
        n = scaled_dose * k1 * k2
        return n * (math.exp(-k2 * tau) + math.exp(-k1 * tau) * (delta * tau - 1)) / (delta * delta)

    k1k2 = k1 - k2
    k1k3 = k1 - k3
    k2k3 = k2 - k3

    term1 = math.exp(-k1 * tau) / (k1k2 * k1k3)
    term2 = math.exp(-k2 * tau) / (-k1k2 * k2k3)
    term3 = math.exp(-k3 * tau) / (k1k3 * k2k3)

    return scaled_dose * k1 * k2 * (term1 + term2 + term3)


def inj_amount(tau: float, dose_mg: float, p: PKParams) -> float:
    """
    双室储库注射模型
    """
    dose_fast = dose_mg * p.frac_fast
    dose_slow = dose_mg * (1.0 - p.frac_fast)
    amt_fast = _analytic_3c(tau, dose_fast, p.F, p.k1_fast, p.k2, p.k3)
    amt_slow = _analytic_3c(tau, dose_slow, p.F, p.k1_slow, p.k2, p.k3)
    return amt_fast + amt_slow


def one_comp_amount(tau: float, dose_mg: float, p: PKParams) -> float:
    """
    单室 Bateman
    """
    return _bateman_amount(dose_mg, p.F, p.k1_fast, p.k3, tau)


def dual_abs_amount(tau: float, dose_mg: float, p: PKParams) -> float:
    """
    双路径 Bateman
    两个分支均为单室
    无水解
    """
    if dose_mg <= 0:
        return 0.0
    f = max(0.0, min(1.0, p.frac_fast))
    dose_f = dose_mg * f
    dose_s = dose_mg * (1.0 - f)
    amt_f = _bateman_amount(dose_f, p.F_fast, p.k1_fast, p.k3, tau - p.lag_fast_h)
    amt_s = _bateman_amount(dose_s, p.F_slow, p.k1_slow, p.k3, tau - p.lag_slow_h)
    return amt_f + amt_s


def dual_abs_mixed_amount(tau: float, dose_mg: float, p: PKParams) -> float:
    """
    双路径
    快速 -> 三室（水解）
    慢速 -> Bateman
    """
    if dose_mg <= 0:
        return 0.0
    f = max(0.0, min(1.0, p.frac_fast))
    dose_f = dose_mg * f
    dose_s = dose_mg * (1.0 - f)
    amt_f = _analytic_3c(tau - p.lag_fast_h, dose_f, p.F_fast, p.k1_fast, p.k2, p.k3)
    amt_s = _bateman_amount(dose_s, p.F_slow, p.k1_slow, p.k3, tau - p.lag_slow_h)
    return amt_f + amt_s


def patch_amount(tau: float, dose_mg: float, wear_h: float, p: PKParams) -> float:
    """
    贴片模型
    """
    if p.rate_mg_h > 0:
        effective_rate = p.rate_mg_h * p.F
        if tau <= wear_h:
            return effective_rate / p.k3 * (1.0 - math.exp(-p.k3 * tau))
        else:
            amt_at_removal = effective_rate / p.k3 * (1.0 - math.exp(-p.k3 * wear_h))
            dt = tau - wear_h
            return amt_at_removal * math.exp(-p.k3 * dt)
    else:
        one_comp_p = PKParams(
            k1_fast=p.k1_fast, k1_slow=0,
            k2=p.k2, k3=p.k3, F=p.F,
            frac_fast=1.0,
            F_fast=p.F, F_slow=p.F,
            lag_fast_h=p.lag_fast_h, lag_slow_h=p.lag_slow_h,
            rate_mg_h=0,
        )
        amt_under = one_comp_amount(tau, dose_mg, one_comp_p)
        if tau > wear_h:
            amt_at_removal = one_comp_amount(wear_h, dose_mg, one_comp_p)
            return amt_at_removal * math.exp(-p.k3 * (tau - wear_h))
        return amt_under


def amount_at_time(tau: float, event: DoseEvent, params: PKParams, wear_h: float) -> float:
    """
    tau 时中央室的药量, mg
    """
    if tau < 0:
        return 0.0
    if event.route == Route.injection:
        return inj_amount(tau, event.dose_mg, params)
    elif event.route == Route.gel:
        return one_comp_amount(tau, event.dose_mg, params)
    elif event.route == Route.oral:
        if params.k1_slow > 0 and params.frac_fast > 0:
            return dual_abs_amount(tau, event.dose_mg, params)
        return one_comp_amount(tau, event.dose_mg, params)
    elif event.route == Route.sublingual:
        if params.k2 > 0:
            return dual_abs_mixed_amount(tau, event.dose_mg, params)
        else:
            return dual_abs_amount(tau, event.dose_mg, params)
    elif event.route == Route.patch_apply:
        return patch_amount(tau, event.dose_mg, wear_h, params)
    elif event.route == Route.patch_remove:
        return 0.0
    return 0.0


@dataclass
class SimulationResult:
    time_h: List[float]
    concentrations: List[float]
    auc: float
    hormone: Hormone
    concentration_unit: str

    def concentration_at(self, hour: float) -> Optional[float]:
        """
        任意时间浓度
        """
        if not self.time_h or len(self.time_h) != len(self.concentrations):
            return None
        t = self.time_h
        c = self.concentrations
        if hour <= t[0]:
            return c[0]
        if hour >= t[-1]:
            return c[-1]
        lo, hi = 0, len(t) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if t[mid] == hour:
                return c[mid]
            elif t[mid] < hour:
                lo = mid
            else:
                hi = mid
        t0, t1 = t[lo], t[hi]
        c0, c1 = c[lo], c[hi]
        if t1 <= t0:
            return c0
        ratio = (hour - t0) / (t1 - t0)
        return c0 + (c1 - c0) * ratio


def simulate_timeline(
        events: List[DoseEvent],
        hormone: Hormone,
        body_weight_kg: float = 70.0,
        history_padding_hours: float = 24.0,
        forecast_hours: float = 24.0 * 14.0,
        number_of_steps: int = 1000,
) -> SimulationResult:
    """
    跨多个给药事件模拟
    """
    if not events:
        return SimulationResult([], [], 0.0, hormone, "")

    sim_events = [e for e in events if compound_hormone(e.compound) == hormone]
    if not sim_events:
        return SimulationResult([], [], 0.0, hormone, "")

    start_time = min(e.time_h for e in sim_events) - history_padding_hours
    end_time = max(e.time_h for e in sim_events) + forecast_hours

    if start_time >= end_time or number_of_steps <= 1:
        return SimulationResult([], [], 0.0, hormone, "")

    precomputed = []
    for event in sim_events:
        if event.route == Route.patch_remove:
            continue
        params = resolve_params(event)
        wear_h = float('inf')
        if event.route == Route.patch_apply:
            for other in sim_events:
                if other.route == Route.patch_remove and other.time_h > event.time_h:
                    wear_h = other.time_h - event.time_h
                    break
        precomputed.append((event, params, wear_h))

    vd = vd_ml(body_weight_kg, hormone)
    if vd <= 0:
        return SimulationResult([], [], 0.0, hormone, "")

    unit = HORMONE_PARAMS[hormone].concentration_unit
    scale = concentration_scale(hormone, unit)

    step_size = (end_time - start_time) / (number_of_steps - 1)
    times: List[float] = []
    concs: List[float] = []
    auc = 0.0
    prev_conc = 0.0

    for i in range(number_of_steps):
        t = start_time + i * step_size
        total_mg = 0.0
        for event, params, wear_h in precomputed:
            tau = t - event.time_h
            total_mg += amount_at_time(tau, event, params, wear_h)
        conc = total_mg * scale / vd
        times.append(t)
        concs.append(conc)
        if i > 0:
            auc += 0.5 * (conc + prev_conc) * step_size
        prev_conc = conc

    return SimulationResult(
        time_h=times,
        concentrations=concs,
        auc=auc,
        hormone=hormone,
        concentration_unit=unit.value,
    )


def compute_concentration_at(
        hour: float,
        events: List[DoseEvent],
        hormone: Hormone,
        body_weight_kg: float = 70.0,
) -> float:
    """
    特定时间点的浓度
    """
    sim_events = [e for e in events if compound_hormone(e.compound) == hormone]
    if not sim_events:
        return 0.0

    vd = vd_ml(body_weight_kg, hormone)
    if vd <= 0:
        return 0.0

    unit = HORMONE_PARAMS[hormone].concentration_unit
    scale = concentration_scale(hormone, unit)
    total_mg = 0.0

    for event in sim_events:
        if event.route == Route.patch_remove:
            continue
        tau = hour - event.time_h
        if tau < 0:
            continue
        params = resolve_params(event)
        wear_h = float('inf')
        if event.route == Route.patch_apply:
            for other in sim_events:
                if other.route == Route.patch_remove and other.time_h > event.time_h:
                    wear_h = other.time_h - event.time_h
                    break
        total_mg += amount_at_time(tau, event, params, wear_h)

    return total_mg * scale / vd


__all__ = [
    "resolve_params", "PKParams",
    "inj_amount", "one_comp_amount", "dual_abs_amount",
    "dual_abs_mixed_amount", "patch_amount", "amount_at_time",
    "SimulationResult", "simulate_timeline", "compute_concentration_at",
]
