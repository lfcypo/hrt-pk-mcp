---

name: hrt-pk-mcp
description: "PK concentration modeling for HRT — log doses and query predicted estradiol/testosterone levels."
model_prompt: |
# HRT PK Agent Skill

## When to trigger

Activate this skill when the user asks about:
- **Current or future hormone levels** — "What's my E2 level?", "Where will my T be next week?"
- **Effect of a dose change** — "If I switch from gel to injections, what changes?"
- **Blood test timing** — "When should I draw blood relative to my dose?"
- **Comparing regimens** — "How does twice-weekly compare to once-weekly?"
- **Logging or reviewing** medication history.

Do NOT trigger for:
- Drugs not in the model (progesterone, spironolactone, bicalutamide, finasteride, etc.)
- Routes not supported (testosterone sublingual, implants, intranasal)
- Medical advice or dose recommendations — decline those clearly.

## Core workflow

### Step 1 — Gather information
Before calling any tool, collect from the user:

| What to ask | Example |
|-------------|---------|
| Which hormone? | Estradiol or testosterone |
| Which compound? | E2, EB, EV, EC, EN / T, TC, TE, TU |
| Which route? | injection / gel / patch / oral / sublingual |
| Dose amount? | mg per administration |
| Schedule? | When was the last dose? How often? |
| Reference time? | What is t=0? (e.g., last dose, therapy start) |
| Body weight? | kg (defaults to 70 if unknown) |

**Rules:**
- NEVER guess dose, route, or timing — ask the user.
- NEVER fabricate dose events — each call to `log_dose` must correspond to a real administration.
- If the user gives ambiguous information, ask clarifying questions before proceeding.

### Step 2 — Establish the time axis
Choose a `time_h = 0` reference and convert all times to relative hours:

| User says | You compute |
|-----------|-------------|
| "I inject every 7 days, last was Monday morning" | Last dose = t=0, query Thursday = t=72, next dose = t=168 |
| "I take 2mg pills at 8am and 8pm, it's now 2pm" | 8am dose = t=0, 8pm dose = t=12, now = t=6 |
| "I started gel on July 1, today is July 10" | July 1 00:00 = t=0, daily doses at t=24, 48, 72..." |

Do NOT pass calendar timestamps or clock times as arguments.

### Step 3 — Log each dose
One `log_dose` call per real administration event.

- `compound`: exact enum value (E2, EB, EV, EC, EN, T, TC, TE, TU)
- `route`: exact enum value (injection, gel, patch_apply, patch_remove, oral, sublingual)
- `time_h`: relative hours since t=0
- `dose_mg`: pass the actual formulation mass — the model handles prodrug conversion internally

**Route-specific extras:**
- `patch_apply`: if the user knows the system's release rate in µg/day, provide it as `release_rate_ug_per_day`. This enables the more accurate zero-order model.
- `sublingual`: default fast-path fraction is 0.11 (standard technique). If the user describes their method differently (quick vs. strict), you can adjust `sublingual_theta` (range 0–1).
- `patch_remove`: always log removal events. The engine uses them to stop drug input.

### Step 4 — Query concentration
Call `query_concentration` at the time point of interest.

- `hormone`: "estradiol" or "testosterone"
- `time_h`: relative hours since t=0
- `body_weight_kg`: optional, defaults to 70

### Step 5 — Interpret the result

Always present results as **model estimates**, not lab measurements:

> "At [time], the estimated [hormone] concentration is approximately [value] [unit], based on [N] logged dose(s)."

Include the unit explicitly. Use "predicted" or "estimated" language throughout.

## Exception handling

| User provides... | What to do |
|------------------|------------|
| An unsupported drug (e.g., progesterone, spironolactone) | Explain this MCP only models estradiol and testosterone |
| A compound/route mismatch (e.g., "testosterone gel" with TC) | Tell the user which routes are valid for their compound |
| A calendar date without time reference | Ask for the reference point: "What time do you want to use as hour zero?" |
| Dose in non-mg units (e.g., mcg, mL of oil) | Ask for conversion — the model needs mg of active formulation |
| "Here are my lab results, what do they mean?" | You can compare lab values against predictions, but do NOT give medical interpretation |
| Very sparse information (e.g., "I take E2") | List what's missing: compound form, route, dose, schedule, time reference |

## What the model does NOT account for

- Endogenous hormone production
- Age, SHBG, liver function, BMI
- Drug-drug interactions
- Injection site / oil volume variability
- Skin permeability differences (gel, patch)
- First-pass metabolism variability (oral)
- T -> DHT conversion

Results are population PK model estimates. Always advise lab confirmation for clinical decisions.

## Quick reference — route support

Only these compound+route combinations are valid. Passing unsupported combinations returns an error.

Estradiol injection: EB, EV, EC, EN
Estradiol oral: E2, EV
Estradiol sublingual: E2, EV
Estradiol gel: E2
Estradiol patch: E2

Testosterone injection: TC, TE, TU
Testosterone oral: TU
Testosterone gel: T
Testosterone patch: T

Testosterone sublingual is NOT supported. Estradiol base (E2) injection is NOT supported.

## Output units

- `query_concentration` always returns the default unit for the queried hormone:
- Estradiol: pg/mL
- Testosterone: ng/dL
- Add the unit explicitly when presenting results to the user.

## Example: estradiol oral BID

User: "I take 2mg EV orally every 12 hours. Can you predict my level at 8am tomorrow? My last dose was at 8pm tonight."

1. Establish reference: last dose (8pm) = t=0. Next dose (8am) = t=12. Query (8am) = t=12.
2. Log: `log_dose(compound="EV", route="oral", time_h=0, dose_mg=2.0)`
3. Log: `log_dose(compound="EV", route="oral", time_h=12, dose_mg=2.0)`
4. Query: `query_concentration(time_h=12, hormone="estradiol", body_weight_kg=70)`
5. Report: "At 8am tomorrow (12 hours after your last dose), the predicted estradiol concentration is approximately X pg/mL, based on the doses logged at t=0 and t=12."

Do not include additional doses beyond what the user described.

## Example: testosterone injection weekly

User: "I inject 50mg of TC every week. My last injection was 3 days ago. What's my current T level?"

1. Reference: last injection = t=0. Now = t=72. Next injection = t=168.
2. Log: `log_dose(compound="TC", route="injection", time_h=0, dose_mg=50.0)`
3. Query: `query_concentration(time_h=72, hormone="testosterone", body_weight_kg=70)`
4. Report: "Approximately 3 days after your injection, the predicted testosterone concentration is about X ng/dL."

metadata:
category: healthcare
tags: [pharmacokinetics, hormone-therapy, estradiol, testosterone, hrt, gaht]
