# -*- coding: utf-8 -*-
# final_full_year_day_forecast_mundane_events.py
# Builds 7-day + annual JSONs with aligned start date and writes them to folders.
# Uses PyEphem for lunation dates, your ephemeris JSON for signs.

from __future__ import annotations

import json, re, os, requests
from datetime import datetime, timedelta, timezone, date
from typing import Dict, List, Any, Optional
from pathlib import Path

# ---------- Third-party ----------
try:
    import ephem as _ephem
except Exception as e:
    raise RuntimeError("PyEphem is required. Install in your venv:  python3 -m pip install ephem requests") from e

# ---------- Live feed URLs (your repos) ----------
RULES_URL   = "https://raw.githubusercontent.com/Stp155906/mundane-cycles/refs/heads/main/rules.json"
ECHOES_URL  = "https://raw.githubusercontent.com/Stp155906/mundane-cycles/refs/heads/main/echo_catalog.json"
ASPECTS_URL = "https://raw.githubusercontent.com/Stp155906/weekly-aspects/refs/heads/main/weekly_aspects.json"
EPHEMERIS_URL_TMPL = "https://raw.githubusercontent.com/Stp155906/ephemeris/refs/heads/main/{year}_ephemeris_with_signs.json"

# ---------- Eclipse heuristics ----------
ECLIPSE_NODE_THRESHOLD_NEW_DEG  = 17.0  # ~solar eclipse window
ECLIPSE_NODE_THRESHOLD_FULL_DEG = 11.0  # ~lunar eclipse window
NODE_KEYS = ["true node", "mean node", "north node", "node", "lunar node", "rahu"]

# ---------- Paths ----------
def _script_dir() -> Path:
    return Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()

def _default_base_dir() -> Path:
    # Prefer /content if present (keeps parity with Colab), else local folder
    return Path("/content") if Path("/content").exists() else _script_dir()

BASE_DIR = Path(os.environ.get("FORECAST_BASE_DIR", str(_default_base_dir())))
DIR_7D     = BASE_DIR / "forecast_7d"
DIR_ANNUAL = BASE_DIR / "forecast_annual"

def ensure_dir(p: str | Path):
    Path(p).mkdir(parents=True, exist_ok=True)

# ---------- JSON helpers ----------
def _strip_json_comments(s: str) -> str:
    s = re.sub(r'(?m)^[ \t]*//.*$', '', s)
    s = re.sub(r'/\*[\s\S]*?\*/', '', s)
    s = re.sub(r',\s*([}\]])', r'\1', s)
    return s

def load_json5(url: str) -> Any:
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    txt = r.text
    try:
        return json.loads(txt)
    except json.JSONDecodeError:
        return json.loads(_strip_json_comments(txt))

# ---------- Canonical + UI utils ----------
_ASPECT_CANON = {
    "□": "square", "square":"square",
    "△": "trine", "trine":"trine",
    "☍": "opposition", "opposition":"opposition",
    "☌": "conjunction", "conjunction":"conjunction",
    "✶": "sextile", "sextile":"sextile",
    "semi-sextile":"semi-sextile", "semisextile":"semi-sextile"
}
_ASPECT_SYMBOL = {"square":"□","trine":"△","opposition":"☍","conjunction":"☌","sextile":"✶","semi-sextile":"⚺"}

def canonical_aspect(s: str) -> str:
    return _ASPECT_CANON.get((s or "").strip().lower(), (s or "").strip().lower())

def symbol_for(aspect_name: str) -> str:
    return _ASPECT_SYMBOL.get(canonical_aspect(aspect_name), aspect_name)

def tcase(x: str) -> str:
    return (x or "").strip().title()

def pretty_date_str(yyyy_mm_dd: str) -> str:
    try:
        dt = datetime.strptime(yyyy_mm_dd, "%Y-%m-%d").date()
        try:    return dt.strftime("%b %-d, %Y")   # Linux/mac
        except: return dt.strftime("%b %#d, %Y")   # Windows
    except:
        return yyyy_mm_dd

def render_template(s: str, d: str) -> str:
    dl = pretty_date_str(d)
    return (s or "").replace("{{date_local}}", dl).replace("{date_local}", dl)

def build_subtitle_from_aspect(a: Dict[str,Any]) -> str:
    b1 = tcase(a.get("body1","")); b2 = tcase(a.get("body2",""))
    pos = a.get("positions") or {}
    s1 = (pos.get(a.get("body1",""),{}) or {}).get("sign","")
    s2 = (pos.get(a.get("body2",""),{}) or {}).get("sign","")
    sym = symbol_for(a.get("aspect_name",""))
    orb = float(a.get("orb_deg", 0.0)); ph = (a.get("phase") or "")
    left  = f"{b1}" + (f" in {s1}" if s1 else "")
    right = f"{b2}" + (f" in {s2}" if s2 else "")
    tail  = f" — orb {orb:.1f}°" + (f" ({ph})" if ph else "")
    return f"{left} {sym} {right}{tail}"

# ---------- Rule matching (for aspect cards) ----------
def rule_matches_aspect(rule: Dict[str,Any], a: Dict[str,Any]) -> bool:
    if rule.get("type") != "aspect_theme":
        return False
    m = (rule.get("match") or {})

    # aspect name must match
    want_aspect = canonical_aspect(m.get("aspect_name_exact",""))
    have_aspect = canonical_aspect(a.get("aspect_name",""))
    if want_aspect and want_aspect != have_aspect:
        return False

    # exact body pair (order-insensitive)
    rb = [tcase(x) for x in (m.get("bodies_exact") or [])]
    if rb:
        pair = sorted([tcase(a.get("body1","")), tcase(a.get("body2",""))])
        if pair != sorted(rb):
            return False

    # numeric gates
    if float(a.get("orb_deg", 999.0)) > float(m.get("orb_max_deg", 999.0)):
        return False
    if float(a.get("importance_score", 0.0)) < float(m.get("importance_min", 0.0)):
        return False

    # optional sign requirements
    if m.get("require_signs"):
        need1 = tcase(m.get("body1_sign","")); need2 = tcase(m.get("body2_sign",""))
        pos = a.get("positions") or {}
        p1 = tcase((pos.get(a.get("body1",""),{}) or {}).get("sign",""))
        p2 = tcase((pos.get(a.get("body2",""),{}) or {}).get("sign",""))
        if not ((p1==need1 and p2==need2) or (p1==need2 and p2==need1)):
            return False

    return True

# ---------- Data fetchers ----------
def fetch_rules_and_echoes():
    rules_doc = load_json5(RULES_URL) or {}
    echoes_doc = load_json5(ECHOES_URL) or {}
    rules  = rules_doc.get("rules", []) if isinstance(rules_doc, dict) else []
    echoes = echoes_doc.get("echoes", []) if isinstance(echoes_doc, dict) else []
    print(f"TE ▶︎ loaded rules={len(rules)} echoes={len(echoes)}")
    return rules, echoes

def fetch_weekly_aspects_preserve() -> List[Dict[str,Any]]:
    data = load_json5(ASPECTS_URL)
    weekly: List[Dict[str,Any]] = []
    if isinstance(data, dict) and "weekly_aspects" in data and isinstance(data["weekly_aspects"], list):
        for item in data["weekly_aspects"]:
            if isinstance(item, dict) and "date" in item and isinstance(item.get("aspects"), list):
                weekly.append({"date": item["date"], "aspects": item["aspects"]})
    elif isinstance(data, dict):
        for d, lst in sorted(data.items()):
            if isinstance(lst, list):
                weekly.append({"date": d, "aspects": lst})
    elif isinstance(data, list):
        by_day: Dict[str, List[Dict[str,Any]]] = {}
        for a in data:
            if isinstance(a, dict) and "date" in a:
                by_day.setdefault(a["date"], []).append(a)
        for d in sorted(by_day.keys()):
            weekly.append({"date": d, "aspects": by_day[d]})
    print(f"TE ▶︎ weekly_aspects days={len(weekly)}")
    return weekly

# ---------- Ephemeris cache + loaders ----------
_DATE_RE  = re.compile(r"^\d{4}-\d{2}-\d{2}$")
_MONTH_RE = re.compile(r"^\d{4}-\d{2}$")
EPHEMERIS_CACHE: Dict[int, Dict[str, Dict[str,str]]] = {}

def _month_name_to_num(mname: str) -> Optional[int]:
    s = (mname or "").strip()
    for cand in (s, s.title(), s.capitalize()):
        for fmt in ("%B", "%b"):  # October / Oct
            try:
                return datetime.strptime(cand, fmt).month
            except:
                pass
    return None

def fetch_ephemeris_for_year(year: int) -> Dict[str, Dict[str,str]]:
    """Load and normalize your ephemeris for a year; return from cache on repeats."""
    if year in EPHEMERIS_CACHE:
        cached = EPHEMERIS_CACHE[year]
        print(f"TE ▶︎ ephemeris (cache) days={len(cached)} (year={year})")
        return cached

    url = EPHEMERIS_URL_TMPL.format(year=year)
    data = load_json5(url) or {}
    norm: Dict[str, Dict[str,str]] = {}

    if isinstance(data, dict):
        for key, val in data.items():
            # Pattern A: {"YYYY-MM-DD": {...}}
            if isinstance(val, dict) and _DATE_RE.match(key):
                norm[key] = {k.lower(): v for k, v in val.items()}
                continue
            # Pattern B: {"YYYY-MM": {"1": {...}, ...}}
            if isinstance(val, dict) and _MONTH_RE.match(key):
                month_prefix = key  # YYYY-MM
                for day_key, bodies in val.items():
                    if not isinstance(bodies, dict): continue
                    try:
                        dd = int(str(day_key).strip())
                        iso = f"{month_prefix}-{dd:02d}"
                        norm[iso] = {k.lower(): v for k, v in bodies.items()}
                    except:
                        pass
                continue
            # Pattern C: {"January": {"1": {...}}}
            mnum = _month_name_to_num(key)
            if isinstance(val, dict) and mnum is not None:
                for day_key, bodies in val.items():
                    if not isinstance(bodies, dict): continue
                    try:
                        dd = int(str(day_key).strip())
                        iso = f"{year}-{mnum:02d}-{dd:02d}"
                        norm[iso] = {k.lower(): v for k, v in bodies.items()}
                    except:
                        pass

    EPHEMERIS_CACHE[year] = norm
    print(f"TE ▶︎ ephemeris days={len(norm)} (year={year})")
    return norm

def ensure_ephemeris_loaded(eph_by_day: Dict[str,Dict[str,str]], start: date, days: int):
    """Merge year data from cache for the whole window (+buffer)."""
    buffer_days = 90
    end_for_eph = start + timedelta(days=days + buffer_days)
    for y in range(start.year, end_for_eph.year + 1):
        eph_by_day.update(fetch_ephemeris_for_year(y))

def extract_sign_from_ephemeris(pos: str) -> str:
    return (pos or "").split(" ")[0].title().strip()

# ---------- Geometry helpers ----------
SIGNS = ["Aries","Taurus","Gemini","Cancer","Leo","Virgo","Libra","Scorpio","Sagittarius","Capricorn","Aquarius","Pisces"]
SIGN_TO_INDEX = {s:i for i,s in enumerate(SIGNS)}

def _parse_lon_deg(pos: str) -> Optional[float]:
    if not pos: return None
    clean = pos.replace("’"," ").replace("'", " ").replace("°"," ")
    m = re.search(r"([A-Za-z]+)\s+(\d+)(?:\s+(\d+))?", clean)
    if not m: return None
    sign = m.group(1).title()
    deg = int(m.group(2))
    minute = int(m.group(3)) if m.group(3) else 0
    idx = SIGN_TO_INDEX.get(sign)
    if idx is None: return None
    return idx*30 + deg + minute/60.0

def _body_lon(eph_by_day: Dict[str,Dict[str,str]], day: str, body: str) -> Optional[float]:
    row = eph_by_day.get(day) or {}
    pos = row.get(body.lower())
    return _parse_lon_deg(pos)

def _angle_diff(a: float, b: float) -> float:
    d = (b - a) % 360.0
    if d >= 180.0: d -= 360.0
    return d

def _node_lon(eph_by_day: Dict[str,Dict[str,str]], date_iso: str) -> Optional[float]:
    row = eph_by_day.get(date_iso) or {}
    for k in NODE_KEYS:
        if k in row:
            return _parse_lon_deg(row[k])
    return None

# Eclipse tagging with Node if present; otherwise leave as-is (fallback added later in list tagging)
ECLIPTIC_LAT_THRESH_NEW_DEG  = 1.6
ECLIPTIC_LAT_THRESH_FULL_DEG = 1.0

def _maybe_tag_eclipse(eph_by_day: Dict[str,Dict[str,str]], ev: Dict[str,Any]) -> Dict[str,Any]:
    node = _node_lon(eph_by_day, ev["date"])
    sun  = _body_lon(eph_by_day, ev["date"], "sun")
    if node is not None and sun is not None:
        sep = abs(_angle_diff(sun, node))
        out = dict(ev)
        if ev["phase"] == "new"  and sep <= ECLIPSE_NODE_THRESHOLD_NEW_DEG:  out["kind"] = "solar_eclipse"
        if ev["phase"] == "full" and sep <= ECLIPSE_NODE_THRESHOLD_FULL_DEG: out["kind"] = "lunar_eclipse"
        return out

    # Fallback: use PyEphem ecliptic latitude at local noon of the event date
    try:
        dt_mid = datetime.strptime(ev["date"], "%Y-%m-%d") + timedelta(hours=12)
        m = _ephem.Moon(dt_mid)
        lat_deg = abs(float(_ephem.Ecliptic(m).lat) * 180.0 / 3.141592653589793)
        out = dict(ev)
        if ev["phase"] == "new"  and lat_deg <= ECLIPTIC_LAT_THRESH_NEW_DEG:  out["kind"] = "solar_eclipse"
        if ev["phase"] == "full" and lat_deg <= ECLIPTIC_LAT_THRESH_FULL_DEG: out["kind"] = "lunar_eclipse"
        return out
    except Exception:
        return ev

def _phase_name_from_angle(delta_deg: float) -> str:
    a = delta_deg % 360.0
    if   a < 22.5 or a >= 337.5: return "new"
    elif a < 67.5:  return "waxing crescent"
    elif a < 112.5: return "first quarter"
    elif a < 157.5: return "waxing gibbous"
    elif a < 202.5: return "full"
    elif a < 247.5: return "waning gibbous"
    elif a < 292.5: return "last quarter"
    else:           return "waning crescent"

# ---------- Lunations via PyEphem ----------
def _iso(d: str) -> date:
    return datetime.strptime(d, "%Y-%m-%d").date()

def _ephem_next_on_or_after(start_iso: str, target: str) -> dict:
    """Next new/full moon on or AFTER the start date (UTC)."""
    start_dt = datetime.strptime(start_iso, "%Y-%m-%d")
    func = _ephem.next_new_moon if target == "new" else _ephem.next_full_moon
    dt = func(start_dt).datetime()  # UTC
    if dt.date() < start_dt.date():
        dt = func(start_dt + timedelta(days=1)).datetime()
    return {"date": dt.date().isoformat(), "phase": target, "kind": target, "estimated": False}

def _closest_phase_on_or_after(start_iso: str, eph_by_day: Dict[str,Dict[str,str]], target: str, horizon: int = 60) -> dict:
    # Signature kept for compatibility; eph_by_day/horizon unused here
    return _ephem_next_on_or_after(start_iso, target)

def list_lunations_for_window(eph_by_day: Dict[str,Dict[str,str]], start_iso: str, end_iso: str) -> List[Dict[str,Any]]:
    """Walk forward selecting whichever (new/full) comes next; tag eclipses; dedup; sort."""
    start = _iso(start_iso)
    end   = _iso(end_iso)
    events, seen = [], set()
    cur = start
    while cur <= end + timedelta(days=1):
        new_d  = _ephem.next_new_moon(cur).datetime().date()
        full_d = _ephem.next_full_moon(cur).datetime().date()
        pick_date, pick_phase = (new_d, "new") if new_d <= full_d else (full_d, "full")
        if pick_date > end: break
        key = (pick_date.isoformat(), pick_phase)
        if key not in seen:
            seen.add(key)
            ev = {"date": pick_date.isoformat(), "phase": pick_phase, "kind": pick_phase, "estimated": False}
            ev = _maybe_tag_eclipse(eph_by_day, ev)
            events.append(ev)
        cur = pick_date + timedelta(days=2)
    events.sort(key=lambda e: (e["date"], 0 if e["phase"]=="new" else 1))
    return events

# ---------- Moon rule copy helpers ----------
def find_moon_rule(rules: List[Dict[str,Any]], phase_for_rules: str, sign: str) -> Optional[Dict[str,Any]]:
    phase_for_rules = (phase_for_rules or "").lower()
    sign  = tcase(sign)
    best = None
    for r in rules:
        if r.get("type") != "moon_phase": continue
        m = r.get("match", {})
        if (m.get("phase","").lower() != phase_for_rules): continue
        need_sign = tcase(m.get("sign",""))
        if need_sign and sign == need_sign:
            if (best is None) or (float(r.get("priority",0)) > float(best.get("priority",0))):
                best = r
    if best: return best
    for r in sorted([x for x in rules if x.get("type")=="moon_phase"], key=lambda z: float(z.get("priority",0)), reverse=True):
        m = r.get("match", {})
        if (m.get("phase","").lower() == phase_for_rules) and not m.get("sign"):
            return r
    return None

def build_moon_copy_from_rules(rules, echoes, phase_for_rules, sign, date_iso, default_label) -> Dict[str,Any]:
    r = find_moon_rule(rules, phase_for_rules, sign)
    if r:
        title_tmpl = (r.get("copy", {}) or {}).get("title_template") or (
            f"{tcase(sign)} {default_label} → {{date_local}}" if sign else f"{default_label} → {{date_local}}"
        )
        title = render_template(title_tmpl, date_iso)
        future_theme = (r.get("copy", {}) or {}).get("future_theme", "")
        echo_ids = r.get("echo_ids") or []
        echo_map = {e.get("id"): e for e in echoes}
        echo_objs = [echo_map[eid] for eid in echo_ids if eid in echo_map]
        return {"title": title, "future_theme": future_theme or "Major tide change. Honor the phase’s tone.", "echoes": echo_objs, "rule_id": r.get("id","")}
    base = f"{tcase(sign)} {default_label}" if sign else default_label
    return {"title": f"{base} → {pretty_date_str(date_iso)}", "future_theme": "This upcoming lunation is a reset/culmination point. Notice what’s ending and what wants to begin.", "echoes": [], "rule_id": f"meta.{phase_for_rules}.default"}

# ---------- Build helpers ----------
def daterange(start: date, days: int) -> List[str]:
    return [(start + timedelta(days=i)).isoformat() for i in range(days)]

def event_sign(eph_by_day, date_iso: str, phase: str, kind: Optional[str]) -> str:
    # For New or Solar Eclipse → Sun sign; for Full/Lunar → Moon sign
    use_sun = (phase=="new") or (kind=="solar_eclipse")
    pos = (eph_by_day.get(date_iso, {}) or {}).get("sun" if use_sun else "moon")
    return extract_sign_from_ephemeris(pos) if pos else ""

def phase_label_for_day(eph_by_day, d: str) -> str:
    sun = _body_lon(eph_by_day, d, "sun")
    moon = _body_lon(eph_by_day, d, "moon")
    if sun is None or moon is None: return ""
    return _phase_name_from_angle((moon - sun) % 360.0)

# ---------- Main builder ----------
def build_window(start_date_str: str|None, days: int) -> Dict[str,Any]:
    from datetime import date  # already imported at top
    today = date.today()  # local date, not UTC
    start = datetime.strptime(start_date_str, "%Y-%m-%d").date() if start_date_str else today
    window = daterange(start, days)

    rules, echoes = fetch_rules_and_echoes()
    weekly = fetch_weekly_aspects_preserve()

    # Load ephemeris for all years touched by (window + buffer)
    eph_by_day: Dict[str, Dict[str,str]] = {}
    ensure_ephemeris_loaded(eph_by_day, start, days)

    # Index aspects by day
    by_day: Dict[str, List[Dict[str,Any]]] = {d: [] for d in window}
    for item in weekly:
        d = item.get("date")
        if d in by_day and isinstance(item.get("aspects"), list):
            by_day[d] = [dict(a) for a in item["aspects"]]

    # ----- META: next New & Full from window start (PyEphem) -----
    next_new  = _closest_phase_on_or_after(window[0], eph_by_day, "new")
    next_full = _closest_phase_on_or_after(window[0], eph_by_day, "full")
    next_new  = _maybe_tag_eclipse(eph_by_day, next_new)
    next_full = _maybe_tag_eclipse(eph_by_day, next_full)

    def _meta_block(nxt: Dict[str,Any]) -> Dict[str,Any]:
        d0 = datetime.strptime(window[0], "%Y-%m-%d").date()
        dN = datetime.strptime(nxt["date"], "%Y-%m-%d").date()
        label = "New Moon" if (nxt.get("phase")=="new") else "Full Moon"
        if nxt.get("kind") == "solar_eclipse": label = "Solar Eclipse"
        if nxt.get("kind") == "lunar_eclipse": label = "Lunar Eclipse"
        phase_for_rules = "new" if (nxt.get("phase")=="new" or nxt.get("kind")=="solar_eclipse") else "full"
        ev_sign = event_sign(eph_by_day, nxt["date"], nxt["phase"], nxt.get("kind"))
        copy_pack = build_moon_copy_from_rules(rules, echoes, phase_for_rules, ev_sign, nxt["date"], label)
        return {
            "phase": nxt.get("phase"),
            "kind": nxt.get("kind", nxt.get("phase")),
            "date": nxt["date"],
            "days_out": (dN - d0).days,
            "sign": ev_sign,
            "title": copy_pack["title"],
            "future_theme": copy_pack["future_theme"],
            "echoes": copy_pack["echoes"],
            "rule_id": copy_pack["rule_id"],
            "estimated": bool(nxt.get("estimated", False))
        }

    meta: Dict[str,Any] = {}
    meta["next_new_moon"]  = _meta_block(next_new)
    meta["next_full_moon"] = _meta_block(next_full)
    earlier = meta["next_new_moon"] if meta["next_new_moon"]["days_out"] <= meta["next_full_moon"]["days_out"] else meta["next_full_moon"]
    meta["next_lunation"] = {k: earlier[k] for k in earlier}

    result: Dict[str,Any] = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00","Z"),
        "window": {"start": window[0], "days": days},
        "weekly_aspects": [],
        "meta": meta
    }

    # Build per-day cards
    for d in window:
        aspects = by_day.get(d, [])
        aspects_with_cards: List[Dict[str,Any]] = []
        for a in aspects:
            matches = [r for r in rules if r.get("type")=="aspect_theme" and rule_matches_aspect(r, a)]
            matches.sort(key=lambda r: float(r.get("priority", 0)), reverse=True)
            if matches:
                r0 = matches[0]
                title = r0.get("copy",{}).get("title") or render_template(r0.get("copy",{}).get("title_template",""), d)
                card = {
                    "state": "matched",
                    "rule_id": r0.get("id",""),
                    "title": title,
                    "subtitle": build_subtitle_from_aspect(a),
                    "future_theme": r0.get("copy",{}).get("future_theme",""),
                    "echoes": [e for e in echoes if e.get("id") in (r0.get("echo_ids") or [])]
                }
            else:
                sym = symbol_for(a.get("aspect_name",""))
                suffix = " — Harmony Flow" if canonical_aspect(a.get("aspect_name",""))=="trine" else ""
                title = f"{tcase(a.get('body1',''))} {sym} {tcase(a.get('body2',''))}{suffix}"
                card = {
                    "state": "unmatched",
                    "rule_id": "dev.aspect.unmatched",
                    "title": title,
                    "subtitle": build_subtitle_from_aspect(a),
                    "future_theme": "📜 No echoes linked for this aspect. Focus on the current geometry and vibe.",
                    "echoes": []
                }
            b = dict(a); b["card"] = card
            aspects_with_cards.append(b)

        # Moon anchor
        sun_lon = _body_lon(eph_by_day, d, "sun")
        moon_lon = _body_lon(eph_by_day, d, "moon")
        moon_phase = _phase_name_from_angle((moon_lon - sun_lon) % 360.0) if (sun_lon is not None and moon_lon is not None) else ""
        moon_sign = extract_sign_from_ephemeris((eph_by_day.get(d, {}) or {}).get("moon") or "")
        moon_blob = {"phase": moon_phase, "sign": moon_sign}

        moon_anchor_card = {
            "state": "moon_anchor",
            "rule_id": "meta.moon.anchor",
            "title": f"{(moon_phase or 'Moon').title()} (Anchor)",
            "subtitle": None,
            "future_theme": "Collective mood anchor. 📜 No echoes linked today. Focus on the present lunation tone.",
            "echoes": []
        }

        # Per-day countdowns
        def _countdown(day_iso, target: str):
            nxt = _closest_phase_on_or_after(day_iso, eph_by_day, target)
            nxt = _maybe_tag_eclipse(eph_by_day, nxt)
            d0 = datetime.strptime(day_iso, "%Y-%m-%d").date()
            dN = datetime.strptime(nxt["date"], "%Y-%m-%d").date()
            days_out = (dN - d0).days
            label = "New Moon" if (nxt.get("phase")=="new") else "Full Moon"
            if nxt.get("kind") == "solar_eclipse": label = "Solar Eclipse"
            if nxt.get("kind") == "lunar_eclipse": label = "Lunar Eclipse"
            phase_for_rules = "new" if (nxt.get("phase")=="new" or nxt.get("kind")=="solar_eclipse") else "full"
            ev_sign = event_sign(eph_by_day, nxt["date"], nxt["phase"], nxt.get("kind"))
            copy_pack = build_moon_copy_from_rules(rules, echoes, phase_for_rules, ev_sign, nxt["date"], label)
            when_label = "Today" if days_out == 0 else f"in {days_out} days"
            sign_label = f"{ev_sign} " if ev_sign else ""
            return {
                "state": "countdown",
                "rule_id": copy_pack["rule_id"] or f"meta.countdown.{phase_for_rules}",
                "title": f"Next {label} {when_label} — {sign_label}{label}",
                "subtitle": f"Exact on {pretty_date_str(nxt['date'])}",
                "future_theme": copy_pack["future_theme"],
                "echoes": copy_pack["echoes"],
                "estimated": bool(nxt.get("estimated", False)),
                "event_date": nxt["date"]
            }

        cd_new  = _countdown(d, "new")
        cd_full = _countdown(d, "full")
        d0 = datetime.strptime(d, "%Y-%m-%d").date()
        d_new  = datetime.strptime(cd_new["event_date"],  "%Y-%m-%d").date()
        d_full = datetime.strptime(cd_full["event_date"], "%Y-%m-%d").date()
        earlier_cd = cd_new if (d_new - d0).days <= (d_full - d0).days else cd_full

        result["weekly_aspects"].append({
            "date": d,
            "moon": moon_blob,
            "aspects": aspects_with_cards,
            "moon_anchor_card": moon_anchor_card,
            "countdown_new_card": cd_new,
            "countdown_full_card": cd_full,
            "countdown_card": earlier_cd
        })

    # ----- Add complete lunation/eclipses list (PyEphem + ephemeris signs) -----
    all_events = list_lunations_for_window(eph_by_day, window[0], window[-1])
    decorated = []
    for ev in all_events:
        label = "New Moon" if ev.get("phase")=="new" else "Full Moon"
        if ev.get("kind") == "solar_eclipse": label = "Solar Eclipse"
        if ev.get("kind") == "lunar_eclipse": label = "Lunar Eclipse"
        phase_for_rules = "new" if (ev.get("phase")=="new" or ev.get("kind")=="solar_eclipse") else "full"
        ev_sign = event_sign(eph_by_day, ev["date"], ev["phase"], ev.get("kind"))
        copy_pack = build_moon_copy_from_rules(rules, echoes, phase_for_rules, ev_sign, ev["date"], label)
        decorated.append({
            "date": ev["date"],
            "phase": ev["phase"],
            "kind": ev.get("kind", ev["phase"]),
            "sign": ev_sign,
            "title": copy_pack["title"],
            "future_theme": copy_pack["future_theme"],
            "echoes": copy_pack["echoes"],
            "rule_id": copy_pack["rule_id"],
            "estimated": bool(ev.get("estimated", False))
        })
    result["lunation_events"] = decorated
    return result

# ---------- Exports ----------
def write_main(result: Dict[str,Any], base_path: str | Path) -> str:
    ensure_dir(base_path)
    p = str(Path(base_path) / "forecast_mundaneevents.json")
    with open(p, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print("✅ wrote", p)
    return p

def export_jsons(result: Dict[str,Any], base_path: str | Path) -> Dict[str,str]:
    ensure_dir(base_path)
    base_path = str(base_path)
    paths = {}

    # 1) All aspects
    all_payload = {
        "generated_utc": result["generated_utc"],
        "window": result["window"],
        "days": [{"date": d["date"], "aspects": d["aspects"]} for d in result["weekly_aspects"]]
    }
    p_all = f"{base_path}/forecast_aspects_all.json"
    with open(p_all, "w", encoding="utf-8") as f:
        json.dump(all_payload, f, ensure_ascii=False, indent=2)
    paths["all_aspects"] = p_all

    # 2) Matched
    matched_days = []
    for d in result["weekly_aspects"]:
        ev = [a for a in d["aspects"] if (a.get("card",{}).get("state")=="matched")]
        matched_days.append({"date": d["date"], "aspects": ev})
    matched_payload = {"generated_utc": result["generated_utc"], "window": result["window"], "days": matched_days}
    p_matched = f"{base_path}/forecast_aspects_matched.json"
    with open(p_matched, "w", encoding="utf-8") as f:
        json.dump(matched_payload, f, ensure_ascii=False, indent=2)
    paths["matched_aspects"] = p_matched

    # 3) Unmatched
    unmatched_days = []
    for d in result["weekly_aspects"]:
        ev = [a for a in d["aspects"] if (a.get("card",{}).get("state")=="unmatched")]
        unmatched_days.append({"date": d["date"], "aspects": ev})
    unmatched_payload = {"generated_utc": result["generated_utc"], "window": result["window"], "days": unmatched_days}
    p_unmatched = f"{base_path}/forecast_aspects_unmatched.json"
    with open(p_unmatched, "w", encoding="utf-8") as f:
        json.dump(unmatched_payload, f, ensure_ascii=False, indent=2)
    paths["unmatched_aspects"] = p_unmatched

    # 4) Moon events (meta + per-day anchor & countdowns)
    moon_days = []
    for d in result["weekly_aspects"]:
        moon_days.append({
            "date": d["date"],
            "moon": d["moon"],
            "moon_anchor_card": d.get("moon_anchor_card"),
            "countdown_new_card": d.get("countdown_new_card"),
            "countdown_full_card": d.get("countdown_full_card"),
            "countdown_card": d.get("countdown_card")
        })
    moon_payload = {"generated_utc": result["generated_utc"], "window": result["window"], "meta": result.get("meta", {}), "days": moon_days}
    p_moon = f"{base_path}/forecast_moon_events.json"
    with open(p_moon, "w", encoding="utf-8") as f:
        json.dump(moon_payload, f, ensure_ascii=False, indent=2)
    paths["moon_events"] = p_moon

    # 5) Eclipses-only list (subset of lunation_events)
    eclipses = [e for e in (result.get("lunation_events") or []) if e.get("kind") in {"solar_eclipse","lunar_eclipse"}]
    eclipse_payload = {"generated_utc": result["generated_utc"], "window": result["window"], "eclipses": eclipses}
    p_eclipses = f"{base_path}/forecast_eclipses.json"
    with open(p_eclipses, "w", encoding="utf-8") as f:
        json.dump(eclipse_payload, f, ensure_ascii=False, indent=2)
    paths["eclipses"] = p_eclipses

    # Console summary
    total_all = sum(len(d["aspects"]) for d in all_payload["days"])
    total_matched = sum(len(d["aspects"]) for d in matched_payload["days"])
    total_unmatched = sum(len(d["aspects"]) for d in unmatched_payload["days"])
    print(f"🗂 wrote {p_all}      • events={total_all}")
    print(f"🗂 wrote {p_matched}  • matched={total_matched}")
    print(f"🗂 wrote {p_unmatched}• unmatched={total_unmatched}")
    print(f"🗂 wrote {p_moon}     • days={len(moon_payload['days'])}  (meta keys={list(moon_payload.get('meta',{}).keys())})")
    print(f"🗂 wrote {p_eclipses} • eclipses={len(eclipses)}")
    return paths

# ---------- Alignment self-check (prints only) ----------
def verify_alignment(result_7d: Dict[str,Any], result_annual: Dict[str,Any], today_iso: str):
    m7 = result_7d.get("meta", {})
    n7 = m7.get("next_new_moon", {})
    f7 = m7.get("next_full_moon", {})
    print("\n— 7-DAY META —")
    print("Next New :", {k:n7.get(k) for k in ["date","sign","title","estimated"]})
    print("Next Full:", {k:f7.get(k) for k in ["date","sign","title","estimated"]})

    def iso(d: str) -> date:
        return datetime.strptime(d, "%Y-%m-%d").date()
    evs = [e for e in (result_annual.get("lunation_events") or []) if "date" in e and iso(e["date"]) >= iso(today_iso)]
    evs.sort(key=lambda e: e["date"])

    print(f"\n— ANNUAL (next 8 lunations on/after {today_iso}) —")
    for e in evs[:8]:
        print(f"  {e['date']}: {e.get('sign','')} {e['phase']}  • {e.get('title','')}  (kind={e.get('kind')}, est={e.get('estimated')})")

    probs = []
    def _check(tag, d):
        if not d:
            probs.append(f"{tag}: missing date"); return
        if not any(x["date"] == d for x in evs):
            probs.append(f"{tag}: {d} not found among upcoming annual lunations")

    _check("Next New (7-day)", n7.get("date"))
    _check("Next Full (7-day)", f7.get("date"))

    missing_sign = [e["date"] for e in evs[:8] if not (e.get("sign") or "").strip()]
    if missing_sign:
        probs.append("Annual lunations missing sign for: " + ", ".join(missing_sign))

    print("\n— VERDICT —")
    if probs:
        for p in probs: print("⚠️", p)
    else:
        print("✅ 7-day and annual dates line up; no missing signs detected in upcoming items.")

# ---------- Main ----------
def main():
    ensure_dir(DIR_7D)
    ensure_dir(DIR_ANNUAL)

    start_env = os.environ.get("FORECAST_START_DATE", "").strip()
    # If set, both windows anchor to that date; if empty, uses today's UTC date (automatic)
    start_for_both = start_env if start_env else None

    print("TE ▶︎ building 7-day window…")
    result_7d = build_window(start_for_both, 7)
    write_main(result_7d, DIR_7D)
    export_jsons(result_7d, DIR_7D)

    print("TE ▶︎ building annual window (366d)…")
    result_annual = build_window(start_for_both, 366)
    write_main(result_annual, DIR_ANNUAL)
    export_jsons(result_annual, DIR_ANNUAL)

    # Alignment self-check (uses the SAME anchor; if not set, uses today's start)
    today_iso = result_annual["window"]["start"]
    verify_alignment(result_7d, result_annual, today_iso)

    # small peek
    def peek(result, label):
        print(f"\n— Preview ({label}) —")
        print("window:", result["window"])
        print("meta.next_new_moon:", result["meta"].get("next_new_moon"))
        print("meta.next_full_moon:", result["meta"].get("next_full_moon"))
        if result.get("lunation_events"):
            first_ev = result["lunation_events"][0]
            print("first lunation event:", {k:first_ev[k] for k in ["date","phase","kind","sign","title"]})

    peek(result_7d, "7d")
    peek(result_annual, "annual")

if __name__ == "__main__":
    # Leave FORECAST_START_DATE unset to use today's date automatically
    # Examples (set via env, not here):
    #   FORECAST_START_DATE=2025-09-30
    #   FORECAST_BASE_DIR=/some/path
    main()
