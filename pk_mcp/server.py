from __future__ import annotations

import logging
from typing import Optional

from mcp.server.fastmcp import FastMCP

from .pk_core import compute_concentration_at
from .pk_params import (
    Compound, Hormone, Route,
    compound_hormone, default_concentration_unit,
)
from .storage import (
    log_dose as storage_log_dose,
    get_events_as_dose_events,
    list_events as storage_list_events,
    clear_events as storage_clear_events,
    remove_event as storage_remove_event,
)

logger = logging.getLogger(__name__)

mcp = FastMCP(
    name="hrt-pk-mcp",
    instructions="""# HRT PK MCP — Agent Instructions

## 1. Purpose

This MCP estimates blood plasma concentrations for estradiol and testosterone
by applying a pharmacokinetic model to logged dosing events.

It does NOT perform medical diagnosis, provide dosage recommendations, or
interpret laboratory results.

---

## 2. When to Use

Call this MCP when the user asks about:
- Predicted current or future hormone level
- Effect of a dose schedule change on concentration
- Timing of blood draws relative to doses
- Logging or reviewing their medication schedule
- Comparing concentration profiles across different regimens

---

## 3. When NOT to Use

Do NOT call this MCP when:
- The drug is not supported (progesterone, spironolactone, bicalutamide, finasteride, etc.)
- The route is not supported (testosterone sublingual, implants, intranasal)
- The user asks for medical advice or dosage recommendations — decline explicitly
- The user provides only calendar timestamps without a defined time reference

---

## 4. Data Model

**dose event** — A single medication administration. Stored persistently.
Each event records: compound, route, time_h, dose_mg, and optional extras.

**simulation** — The engine sums contributions from all logged events for
the queried hormone to produce a concentration estimate at the requested time.

**concentration** — Model output in standard units:
- Estradiol: pg/mL
- Testosterone: ng/dL

Events logged for estradiol do NOT affect testosterone, and vice versa.

---

## 5. Time System

All time values (`time_h`) are **relative hours** since a single reference
point chosen by the Agent in conversation with the user.

**Examples:**

| User statement | Your computation |
|----------------|------------------|
| "I last injected 7 days ago" | Last injection = t=0. Now = t=168 |
| "I take pills at 8am and 8pm, it's now 2pm" | 8am dose = t=0. Now = t=6. Next dose = t=12 |
| "I started therapy on July 1" | July 1 00:00 = t=0 |

- Choose a convenient t=0 (first dose, last dose, or therapy start).
- Convert ALL calendar times to hours from t=0.
- NEVER pass Unix timestamps, ISO 8601 strings, or clock times as time_h.
- If the user gives ambiguous timing, ask clarifying questions before proceeding.

---

## 6. Dose Specification

`dose_mg` = amount of the **formulation mass** actually administered.
The model internally converts prodrugs to active hormone equivalents.

**Do NOT apply molecular weight conversion.** Pass the real-world dose as-is.

| Compound | What to pass as dose_mg |
|----------|------------------------|
| E2 | mg of E2 as applied/taken |
| EB, EV, EC, EN | mg of ester as injected/taken |
| T | mg of T as applied |
| TC, TE, TU | mg of ester as injected/taken |

**Patch example:** A "100 ug/day system" → pass dose_mg=0.1 and
release_rate_ug_per_day=100 separately. The zero-order release model is
more accurate when the release rate is known.

**Gel example:** dose_mg reflects the active hormone content of the gel
packet. Area (area_cm2) is not currently used by the model.

**Sublingual:** The fraction of dose absorbed via the fast mucosal path
(sublingual_theta) defaults to 0.11 (standard technique). Acceptable range
is 0-1. If the user describes an unusual technique, you may pass a custom
value; otherwise omit.

---

## 7. Workflow

### Phase 1 — Assess
Determine whether the request is in scope (§2). If the drug or route is
unsupported, inform the user and stop.

### Phase 2 — Gather
Collect from the user (ask, do not guess):

| Required | Notes |
|----------|-------|
| Hormone | estradiol or testosterone |
| Compound | Which specific compound (E2, EV, TC, etc.) |
| Route | Must be valid for the compound (§12) |
| Dose amount | mg per administration |
| Schedule | Times of doses and frequency |
| Time reference | What is t=0? |

| Optional | Default | Notes |
|----------|---------|-------|
| Body weight | 70 kg | Affects distribution volume |
| Release rate | — | For patch; enables zero-order model |
| Sublingual theta | 0.11 | For sublingual only |

### Phase 3 — Validate
Before calling any tool, confirm:
- Compound is one of: E2, EB, EV, EC, EN, T, TC, TE, TU
- Route is valid for this compound (§12)
- time_h >= 0 and is a relative hour offset
- dose_mg > 0
- release_rate_ug_per_day > 0 if provided
- sublingual_theta in [0,1] if provided

If validation fails, return the error to the user. Do not silently correct.

### Phase 4 — Log
Call `log_dose` once per real administration event.

**Rules:**
- Each call must correspond to a real dose the user described.
- Do NOT fabricate or extrapolate doses the user did not explicitly state.
- Log ALL described doses before querying concentration.
- For patches: log BOTH `patch_apply` and `patch_remove` events. The engine
  pairs each apply with the next remove to stop drug input.

### Phase 5 — Query
Call `query_concentration` at the desired time point with the target hormone.
Include body_weight_kg if the user provided it.

### Phase 6 — Report
Present the result to the user. See §10 for output standards.

---

## 8. Validation and Error Handling

| Condition | Response |
|-----------|----------|
| Unknown compound | List valid compounds |
| Valid compound, invalid route | List valid routes for that compound |
| Unknown hormone | Only estradiol and testosterone are supported |
| time_h < 0 | Reject, must be >= 0 |
| dose_mg <= 0 | Reject, must be positive |
| release_rate_ug_per_day <= 0 when provided | Reject, must be positive |
| sublingual_theta outside [0,1] | Reject |
| Missing required info | Ask the user. Do not guess. |

**Route restrictions:**
- `injection`: ONLY prodrug esters (EB, EV, EC, EN / TC, TE, TU)
- `oral`: E2, EV (estradiol); TU (testosterone)
- `gel`: E2, T only
- `patch_apply` / `patch_remove`: E2, T only
- `sublingual`: E2, EV only (estradiol only — NOT for testosterone)

---

## 9. History Management

- Events persist across sessions until explicitly removed.
- Call `list_dose_events` to review current state.
- To correct an error: call `remove_event` with the event ID, then re-log
  the corrected event.
- `clear_all_events` removes ALL data. Require explicit user confirmation
  before calling this tool. Do NOT call it as a routine reset.

---

## 10. Output Interpretation

Results are **model estimates**, not laboratory measurements.

### Required language
- "The estimated concentration is approximately X pg/mL"
- "At t=Y hours, the predicted level is around Z ng/dL"
- NOT "Your level is X" or "Your blood test would show Y"

### Always include
- The numeric estimate
- The unit (pg/mL or ng/dL)
- Number of dose events included
- That this is a model prediction

### Never
- Describe results as if they were lab values
- Give medical advice or dosage recommendations based on results
- Omit units
- Overstate precision beyond 1 decimal place

---

## 11. Model Limitations

The following are NOT accounted for:
- Endogenous hormone production
- Age, SHBG, liver function, BMI
- Drug-drug interactions
- Injection site or oil volume variability
- Skin permeability differences
- First-pass metabolism variability
- T to DHT conversion

Results are population PK model estimates. Individual results may vary
substantially. Recommend lab confirmation for clinical decisions.

---

## 12. Supported Compounds and Routes

### Estradiol
| Compound | Inj | Oral | SL | Gel | Patch |
|----------|-----|------|----|-----|-------|
| E2 | — | yes | yes | yes | yes |
| EB | yes | — | — | — | — |
| EV | yes | yes | yes | — | — |
| EC | yes | — | — | — | — |
| EN | yes | — | — | — | — |

SL = sublingual. Inj = injection.

### Testosterone
| Compound | Inj | Oral | Gel | Patch |
|----------|-----|------|-----|-------|
| T | — | — | yes | yes |
| TC | yes | — | — | — |
| TE | yes | — | — | — |
| TU | yes | yes | — | — |

---

## 13. Tool Signatures

| Tool | Required | Optional | Returns |
|------|----------|----------|---------|
| log_dose | compound, route, time_h, dose_mg | release_rate_ug_per_day, area_cm2, sublingual_theta | Confirmation + event ID |
| query_concentration | time_h, hormone | body_weight_kg | Estimated concentration, event count, total dose |
| list_dose_events | — | — | Formatted list of all events |
| remove_event | event_id | — | Confirmation |
| clear_all_events | — | — | Confirmation + count removed |

Must have explicit user confirmation before clear_all_events.
""",

)


@mcp.tool()
def log_dose(
        compound: str,
        route: str,
        time_h: float,
        dose_mg: float,
        release_rate_ug_per_day: Optional[float] = None,
        area_cm2: Optional[float] = None,
        sublingual_theta: Optional[float] = None,
) -> str:
    """Record a hormone dosing event.

    Call this whenever the user mentions taking a dose. Stored events are
    used by query_concentration to predict future plasma levels via superposition.

    Args:
        compound: Drug compound — E2, EB, EV, EC, EN (estradiols) or T, TC, TE, TU (testosterones)
        route: Route — injection, gel, patch_apply, patch_remove, oral, sublingual
        time_h: Dose time in hours since user-defined reference (e.g., start of treatment plan)
        dose_mg: Dose in active-hormone-equivalent milligrams (prodrugs converted internally)
        release_rate_ug_per_day: For patch_apply — zero-order release rate in µg/day (optional)
        area_cm2: For gel — application area in cm² (optional, default ~750)
        sublingual_theta: For sublingual — fraction absorbed via fast mucosal path [0-1] (optional, default 0.11)
    """
    try:
        c = Compound(compound)
        r = Route(route)
    except ValueError as e:
        return f"Error: invalid compound or route — {e}"

    if r == Route.sublingual and compound_hormone(c) != Hormone.estradiol:
        return "Error: sublingual route is only supported for estradiol compounds"

    event = storage_log_dose(
        compound=c,
        route=r,
        time_h=time_h,
        dose_mg=dose_mg,
        release_rate_ug_per_day=release_rate_ug_per_day,
        area_cm2=area_cm2,
        sublingual_theta=sublingual_theta,
    )
    return (
        f"Dose recorded: {dose_mg} mg {compound} ({route}) at t={time_h}h. "
        f"Event ID: {event['id']}"
    )


@mcp.tool()
def query_concentration(
        time_h: float,
        hormone: str,
        body_weight_kg: float = 70.0,
) -> str:
    """Query predicted blood plasma concentration at a specific time.

    Computes concentration based on all previously logged dosing events
    for the specified hormone. Uses superposition to combine contributions
    from multiple doses.

    Args:
        time_h: Time point to query (hours since reference, same time base as log_dose)
        hormone: Target hormone -- "estradiol" or "testosterone"
        body_weight_kg: Body weight in kg (default 70.0, affects Vd)
    """
    try:
        h = Hormone(hormone)
    except ValueError:
        return f'Error: invalid hormone "{hormone}". Use "estradiol" or "testosterone".'

    events = get_events_as_dose_events()
    hormone_events = [e for e in events if compound_hormone(e.compound) == h]

    if not hormone_events:
        return (
            f"No dose events recorded for {hormone}. "
            f"Use log_dose to record events first."
        )

    conc = compute_concentration_at(
        hour=time_h,
        events=events,
        hormone=h,
        body_weight_kg=body_weight_kg,
    )
    unit = default_concentration_unit(h).value

    num_doses = len(hormone_events)
    nearest_dose_h = min(
        (e.time_h for e in hormone_events),
        key=lambda t: abs(t - time_h),
    )
    total_dose_mg = sum(e.dose_mg for e in hormone_events)

    return (
        f"At t={time_h}h, predicted {hormone} concentration: "
        f"{conc:.1f} {unit}\n"
        f"(based on {num_doses} dose event(s), "
        f"total {total_dose_mg:.2f} mg active-hormone-equivalent, "
        f"body weight {body_weight_kg:.0f} kg)"
    )


@mcp.tool()
def list_dose_events() -> str:
    """List all recorded dose events, ordered by time."""
    events = storage_list_events()
    if not events:
        return "No dose events recorded."

    lines = ["Recorded dose events (ordered by time):"]
    sorted_events = sorted(events, key=lambda e: e["time_h"])
    for e in sorted_events:
        extras = ""
        if e.get("release_rate_ug_per_day"):
            extras += f", release={e['release_rate_ug_per_day']} µg/day"
        if e.get("area_cm2"):
            extras += f", area={e['area_cm2']} cm²"
        if e.get("sublingual_theta") is not None:
            extras += f", θ={e['sublingual_theta']}"
        lines.append(
            f"  [{e['id']}] t={e['time_h']}h "
            f"{e['dose_mg']} mg {e['compound']} ({e['route']}){extras}"
        )

    return "\n".join(lines)


@mcp.tool()
def clear_all_events() -> str:
    """Clear ALL recorded dose events. Use with caution."""
    count = len(storage_list_events())
    storage_clear_events()
    return f"Cleared {count} dose event(s)."


@mcp.tool()
def remove_event(event_id: str) -> str:
    """Remove a single dose event by its ID.

    Args:
        event_id: The event ID returned by log_dose (or shown in list_dose_events)
    """
    if storage_remove_event(event_id):
        return f"Event {event_id} removed."
    return f"Event {event_id} not found."


def main() -> None:
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
