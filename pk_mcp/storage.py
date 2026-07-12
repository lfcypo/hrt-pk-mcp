from __future__ import annotations

import json
import time
from pathlib import Path
from typing import List, Optional

from .pk_params import Compound, DoseEvent, Route

_PROJECT_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = _PROJECT_DIR / ".hrt-pk-data"
EVENTS_FILE = DATA_DIR / "dose_events.json"


def _ensure_dir() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)


def _load_events() -> List[dict]:
    _ensure_dir()
    if not EVENTS_FILE.exists():
        return []
    try:
        with open(EVENTS_FILE, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return []


def _save_events(events: List[dict]) -> None:
    _ensure_dir()
    with open(EVENTS_FILE, "w") as f:
        json.dump(events, f, indent=2, default=str)


def log_dose(
        compound: Compound,
        route: Route,
        time_h: float,
        dose_mg: float,
        release_rate_ug_per_day: Optional[float] = None,
        area_cm2: Optional[float] = None,
        sublingual_theta: Optional[float] = None,
) -> dict:
    event = {
        "id": f"{int(time.time() * 1000000):x}",
        "compound": compound.value,
        "route": route.value,
        "time_h": time_h,
        "dose_mg": dose_mg,
        "recorded_at": time.time(),
    }
    if release_rate_ug_per_day is not None:
        event["release_rate_ug_per_day"] = release_rate_ug_per_day
    if area_cm2 is not None:
        event["area_cm2"] = area_cm2
    if sublingual_theta is not None:
        event["sublingual_theta"] = sublingual_theta

    events = _load_events()
    events.append(event)
    _save_events(events)
    return event


def list_events() -> List[dict]:
    return _load_events()


def get_events_as_dose_events() -> List[DoseEvent]:
    raw = _load_events()
    result: List[DoseEvent] = []
    for e in raw:
        try:
            compound = Compound(e["compound"])
            route = Route(e["route"])
            event = DoseEvent(
                compound=compound,
                route=route,
                time_h=e["time_h"],
                dose_mg=e["dose_mg"],
                release_rate_ug_per_day=e.get("release_rate_ug_per_day"),
                area_cm2=e.get("area_cm2"),
                sublingual_theta=e.get("sublingual_theta"),
            )
            result.append(event)
        except (KeyError, ValueError):
            continue
    return result


def clear_events() -> None:
    _ensure_dir()
    _save_events([])


def remove_event(event_id: str) -> bool:
    events = _load_events()
    for i, e in enumerate(events):
        if e.get("id") == event_id:
            events.pop(i)
            _save_events(events)
            return True
    return False


__all__ = [
    "log_dose", "list_events", "get_events_as_dose_events",
    "clear_events", "remove_event",
]
