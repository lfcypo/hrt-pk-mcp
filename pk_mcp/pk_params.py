from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Optional


class Hormone(str, Enum):
    estradiol = "estradiol"  # 雌二醇
    testosterone = "testosterone"  # 睾酮


class Compound(str, Enum):
    E2 = "E2"  # Estradiol 17β-雌二醇
    EB = "EB"  # Estradiol Benzoate苯甲酸雌二醇
    EV = "EV"  # Estradiol Valerate 戊酸雌二醇
    EC = "EC"  # Estradiol Cypionate 环戊丙酸雌二醇
    EN = "EN"  # Ethinyl Estradiol 炔雌醇

    T = "T"  # Testosterone 睾酮
    TC = "TC"  # Testosterone Cypionate 环戊丙酸睾酮
    TE = "TE"  # Testosterone Enanthate 庚酸睾酮
    TU = "TU"  # Testosterone Undecanoate 十一酸睾酮


class Route(str, Enum):
    injection = "injection"  # 注射给药
    gel = "gel"  # 凝胶
    oral = "oral"  # 口服
    sublingual = "sublingual"  # 舌下含服
    patch_apply = "patch_apply"  # 贴片
    patch_remove = "patch_remove"  # 移除贴片


class ConcentrationUnit(str, Enum):
    pgPerML = "pg/ML"
    pmolPerL = "pmol/L"
    ngPerDL = "ng/DL"
    ngPerML = "ng/ML"
    nmolPerL = "nmol/L"


@dataclass(frozen=True)
class CompoundInfo:
    full_name: str
    hormone: Hormone
    molecular_weight: float  # 总分子量
    active_molecular_weight: float  # 等效生物活性分子量
    is_prodrug: bool  # 前体药物

    @property
    def to_active_factor(self) -> float:
        """
        折算系数
        """
        return self.active_molecular_weight / self.molecular_weight


COMPOUND_INFO: Dict[Compound, CompoundInfo] = {
    # 雌二醇 C18H24O2 分子量272.38
    Compound.E2: CompoundInfo("Estradiol", Hormone.estradiol, 272.38, 272.38, False),
    Compound.EB: CompoundInfo("Estradiol Benzoate", Hormone.estradiol, 376.50, 272.38, True),
    Compound.EV: CompoundInfo("Estradiol Valerate", Hormone.estradiol, 356.50, 272.38, True),
    Compound.EC: CompoundInfo("Estradiol Cypionate", Hormone.estradiol, 396.58, 272.38, True),
    Compound.EN: CompoundInfo("Estradiol Enanthate", Hormone.estradiol, 384.56, 272.38, True),

    # 睾酮 C19H28O2 分子量288.42
    Compound.T: CompoundInfo("Testosterone", Hormone.testosterone, 288.42, 288.42, False),
    Compound.TC: CompoundInfo("Testosterone Cypionate", Hormone.testosterone, 412.61, 288.42, True),
    Compound.TE: CompoundInfo("Testosterone Enanthate", Hormone.testosterone, 400.59, 288.42, True),
    Compound.TU: CompoundInfo("Testosterone Undecanoate", Hormone.testosterone, 456.70, 288.42, True),
}


@dataclass(frozen=True)
class HormoneParams:
    concentration_unit: ConcentrationUnit
    vd_per_kg: float  # 表观分布容积 (L/kg)
    k_clear: float  # h⁻¹

    k_clear_injection: float
    depot_k1_corr: float  # 缓释储库吸收速率校正系数 default=1.0

    patch_fallback_k1: float
    patch_release_scale: float  # 贴皮零级释放校准缩放系数

    gel_k1: float
    gel_f_max: float  # 凝胶最大生物利用度分数 [0, 1]


HORMONE_PARAMS: Dict[Hormone, HormoneParams] = {
    Hormone.estradiol: HormoneParams(
        concentration_unit=ConcentrationUnit.pgPerML,
        vd_per_kg=2.0,
        k_clear=0.41,
        k_clear_injection=0.041,
        depot_k1_corr=1.0,
        patch_fallback_k1=0.0075,
        patch_release_scale=1.0,
        gel_k1=0.022,
        gel_f_max=0.06,
    ),
    Hormone.testosterone: HormoneParams(
        concentration_unit=ConcentrationUnit.ngPerDL,
        vd_per_kg=2.0,
        k_clear=0.6,
        k_clear_injection=0.03,
        depot_k1_corr=1.0,
        patch_fallback_k1=0.0051,
        patch_release_scale=3.5078419552373226,
        gel_k1=0.05534590723252352,
        gel_f_max=0.22613930825011333,
    ),
}


@dataclass(frozen=True)
class TwoPartDepotParams:
    frac_fast: float
    k1_fast: float  # h⁻¹
    k1_slow: float  # h⁻¹


TWO_PART_DEPOT: Dict[Compound, TwoPartDepotParams] = {
    # 雌二醇
    Compound.EB: TwoPartDepotParams(frac_fast=0.9, k1_fast=0.144, k1_slow=0.114),
    Compound.EV: TwoPartDepotParams(frac_fast=0.4, k1_fast=0.0216, k1_slow=0.0138),
    Compound.EC: TwoPartDepotParams(frac_fast=0.229164549, k1_fast=0.005035046, k1_slow=0.004510574),
    Compound.EN: TwoPartDepotParams(frac_fast=0.05, k1_fast=0.001, k1_slow=0.005),
    # 睾酮
    Compound.TC: TwoPartDepotParams(frac_fast=0.35, k1_fast=0.016, k1_slow=0.0018),
    Compound.TE: TwoPartDepotParams(frac_fast=0.35, k1_fast=0.022, k1_slow=0.0035),
    Compound.TU: TwoPartDepotParams(frac_fast=0.3, k1_fast=0.005, k1_slow=0.001127743154530867),
}

# 注射专用经验缩放系数
# 活性当量给药下拟合说明书 Cmax/Tmax 参考值
FORMATION_FRACTION: Dict[Compound, float] = {
    Compound.EB: 0.10922376473734707,
    Compound.EV: 0.062258288229969413,
    Compound.EC: 0.117255838,
    Compound.EN: 0.12,
    Compound.TC: 0.06775603562678995,
    Compound.TE: 0.09963018136697789,
    Compound.TU: 0.12940928580278235,
}

# 酯化物水解速率常数
HYDROLYSIS_K2: Dict[Compound, float] = {
    Compound.EB: 0.09,
    Compound.EV: 0.07,
    Compound.EC: 0.045,
    Compound.EN: 0.015,
    Compound.TC: 0.06,
    Compound.TE: 0.12,
    Compound.TU: 0.015,
}

# 口服一级吸收速率常数 kAbs (h⁻¹).
ORAL_KABS: Dict[Compound, float] = {
    Compound.E2: 0.32,
    Compound.EV: 0.05,
    Compound.TU: 0.2162055136986597,
}

# 口服生物利用度分数
ORAL_BIOAVAILABILITY: Dict[Compound, float] = {
    Compound.E2: 0.03,
    Compound.EV: 0.03,
    Compound.TU: 0.02698781505574721,
}


@dataclass(frozen=True)
class OralDualAbsorptionParams:
    frac_fast: float  # 快速吸收相剂量占比
    k_abs_fast: float  # h⁻¹ 快速吸收相
    k_abs_slow: float  # h⁻¹ 慢速吸收相
    bioavailability_fast: float  # 快速吸收相生物利用度分数 [0, 1]
    bioavailability_slow: float  # 慢速吸收相生物利用度分数 [0, 1]
    k_clear: float  # h⁻¹ 口服
    lag_hours_fast: float  # 快速吸收相给药滞后时间 (h)
    lag_hours_slow: float  # 慢速吸收相给药滞后时间 (h)


ORAL_DUAL_ABSORPTION: Dict[Compound, OralDualAbsorptionParams] = {
    Compound.TU: OralDualAbsorptionParams(
        frac_fast=1.0,
        k_abs_fast=0.450550912583251,
        k_abs_slow=0.0142806935343998,
        bioavailability_fast=0.025919316335729803,
        bioavailability_slow=0.0,
        k_clear=0.44024417908217306,
        lag_hours_fast=2.75,
        lag_hours_slow=0.0,
    ),
}

# 舌下含服 h⁻¹.
KABS_SL: float = 1.8


@dataclass
class DoseEvent:
    """
    单次给药
    """
    compound: Compound
    route: Route  # 给药途径
    time_h: float  # 仿真起始后的时长 (h)
    dose_mg: float  # 等效活性激素给药剂量 (mg)
    release_rate_ug_per_day: Optional[float] = None  # 透皮贴零级释放速率 (μg/d)
    area_cm2: Optional[float] = None  # 外用凝胶涂抹面积 (cm²)
    sublingual_theta: Optional[float] = None  # 舌下给药快速黏膜通路剂量占比


@dataclass
class PKParams:
    """
    单次给药 药动参数集
    """
    k1_fast: float  # h⁻¹ 快速吸收速率 / 缓释储库快速释放速率
    k1_slow: float  # h⁻¹ 缓释储库慢速释放速率
    k2: float  # h⁻¹ 酯前药水解速率
    k3: float  # h⁻¹ 消除速率常数
    F: float  # 整体生物利用度 / 有效入血生成分数
    frac_fast: float  # 分配至快速通路的剂量占比
    F_fast: float  # 快速吸收通路生物利用度分数
    F_slow: float  # 慢速吸收通路生物利用度分数
    lag_fast_h: float  # 快速通路吸收滞后时间 (h)
    lag_slow_h: float  # 慢速通路吸收滞后时间 (h)
    rate_mg_h: float  # 零级输入释放速率 (mg/h)


def compound_hormone(compound: Compound) -> Hormone:
    return COMPOUND_INFO[compound].hormone


def supported_routes(compound: Compound) -> list[Route]:
    h = compound_hormone(compound)
    if compound == Compound.E2:
        return [Route.patch_apply, Route.patch_remove, Route.gel, Route.oral, Route.sublingual]
    if compound == Compound.EV:
        return [Route.injection, Route.oral, Route.sublingual]
    if compound in (Compound.EB, Compound.EC, Compound.EN):
        return [Route.injection]
    if compound == Compound.T:
        return [Route.patch_apply, Route.patch_remove, Route.gel]
    if compound == Compound.TU:
        return [Route.injection, Route.oral]
    if compound in (Compound.TC, Compound.TE):
        return [Route.injection]
    return []


def concentration_scale(hormone: Hormone, unit: ConcentrationUnit) -> float:
    """
    中央室总药量(mg) 换算的显示血药浓度的缩放系数
    """
    if unit == ConcentrationUnit.pgPerML:
        return 1e9
    if unit == ConcentrationUnit.ngPerDL:
        return 1e8
    if unit == ConcentrationUnit.ngPerML:
        return 1e6
    mw = COMPOUND_INFO[Compound.E2 if hormone == Hormone.estradiol else Compound.T].molecular_weight
    if unit == ConcentrationUnit.pmolPerL:
        return 1e12 / mw
    if unit == ConcentrationUnit.nmolPerL:
        return 1e9 / mw
    return 1e9  # default = pg/mL-like


def default_concentration_unit(hormone: Hormone) -> ConcentrationUnit:
    return HORMONE_PARAMS[hormone].concentration_unit


def vd_ml(body_weight_kg: float, hormone: Hormone) -> float:
    """血浆表观分布容积 (mL)"""
    return HORMONE_PARAMS[hormone].vd_per_kg * body_weight_kg * 1000


__all__ = [
    "Hormone", "Compound", "Route", "ConcentrationUnit",
    "CompoundInfo", "COMPOUND_INFO",
    "HormoneParams", "HORMONE_PARAMS",
    "TwoPartDepotParams", "TWO_PART_DEPOT",
    "FORMATION_FRACTION", "HYDROLYSIS_K2",
    "ORAL_KABS", "ORAL_BIOAVAILABILITY",
    "OralDualAbsorptionParams", "ORAL_DUAL_ABSORPTION",
    "KABS_SL",
    "DoseEvent", "PKParams",
    "compound_hormone", "supported_routes",
    "concentration_scale", "default_concentration_unit", "vd_ml",
]
