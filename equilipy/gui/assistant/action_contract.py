"""Provider-neutral GUI action contract for assistant integrations."""

from __future__ import annotations

import json
from typing import Any

ACTION_BLOCK_START = "<EQUILIPY_ACTIONS>"
ACTION_BLOCK_END = "</EQUILIPY_ACTIONS>"

ALLOWED_ACTION_NAMES = (
    "calculation.create_session",
    "calculation.load_database",
    "calculation.add_module",
    "calculation.set_condition",
    "calculation.set_units",
    "calculation.set_type",
    "calculation.set_batch_condition",
    "calculation.set_batch_grid",
    "solidification.set_batch_grid",
    "calculation.set_nucleation_undercooling",
    "calculation.set_result_columns",
    "calculation.select_phases",
    "calculation.calculate",
    "calculation.calculate_all",
)

ACTION_RESPONSE_SCHEMA: dict[str, Any] = {
    "type": "object",
    "properties": {
        "reply": {
            "type": "string",
            "description": "Brief user-visible response. Do not claim actions ran.",
        },
        "actions": {
            "type": "array",
            "description": "Whitelisted Equilipy GUI actions to request.",
            "items": {
                "type": "object",
                "properties": {
                    "name": {
                        "type": "string",
                        "enum": list(ALLOWED_ACTION_NAMES),
                    },
                    "args": {
                        "type": "object",
                        "additionalProperties": True,
                    },
                },
                "required": ["name"],
                "additionalProperties": False,
            },
        },
    },
    "required": ["reply", "actions"],
    "additionalProperties": False,
}


def action_contract_prompt() -> str:
    """Return provider-neutral instructions for requesting GUI actions."""
    return (
        "You are the right-side assistant inside the Equilipy GUI. "
        "Help with thermodynamic database editing and calculation workflows. "
        "You do not directly control the GUI. You may request only the "
        "whitelisted GUI actions listed below through Equilipy's action bridge. "
        "Do not claim that a GUI action has happened; the Equilipy app will "
        "execute the action and report the result. Do not emit an action block "
        "unless the user explicitly asks for one of these direct GUI actions. "
        "If the user asks to create, set up, or populate a calculation, create "
        "or update the calculation module and leave it visible for review; do "
        "not launch the calculation. If the user asks ambiguously for a "
        "calculation, set up the module first and ask whether to run it. Only "
        "request calculation.calculate or calculation.calculate_all when the "
        "user explicitly asks to run, start, launch, or calculate now in the "
        "current request. Those run actions must include args.confirmed=true. "
        "When a direct GUI action is requested, include a final action block "
        "with valid JSON and no comments. Put the JSON between action tags: the "
        "opening tag is '<' + 'EQUILIPY_ACTIONS' + '>' and the closing tag is "
        "'</' + 'EQUILIPY_ACTIONS' + '>'. The JSON shape is:\n"
        '{"actions":[{"name":"calculation.create_session"},'
        '{"name":"calculation.load_database"},'
        '{"name":"calculation.add_module","args":{"kind":"equilibrium"}},'
        '{"name":"calculation.add_module","args":{"kind":"solidification"}},'
        '{"name":"calculation.set_condition","args":{"T":1000,"P":1,'
        '"composition":{"Al":0.5,"Fe":0.5}}},'
        '{"name":"calculation.set_batch_grid","args":{"T":700,"P":1,'
        '"balance":"Al","axes":{"Cu":[0,100,10],"Mg":[0,100,10],'
        '"Si":[0,100,10]},"total":100}},'
        '{"name":"calculation.select_phases","args":{"phases":["LIQUID"]}}]}\n'
        f"Allowed action names are {', '.join(ALLOWED_ACTION_NAMES)}. "
        "For calculation.add_module, allowed kind values are equilibrium and "
        "solidification. To target an existing calculation session without "
        "asking the user to select it first, include args.session, "
        'args.session_name, or args.session_id, for example {"session":"Session#2"}. '
        "Session targeting is supported by calculation.load_database, "
        "calculation.add_module, and module-edit actions such as "
        "calculation.set_condition. For calculation.set_condition, use args.T, args.P, "
        "and args.composition as a species-to-amount object; solidification may "
        "also use args.delta_t, args.liquid_phase, and args.from_liquidus. "
        "For calculation.set_units, allowed temperature units are K/C/F/R, "
        "pressure units are atm/psi/bar/Pa/kPa, and amount units are the GUI "
        "amount-unit labels. For calculation.set_type, allowed modes are "
        "single and batch; solidification models are scheil, nucleoscheil, "
        "and equilibrium_cooling. For calculation.set_batch_condition, use "
        'args.condition as a column object like {"T":[...],"P":[...],"Al":[...]}. '
        "For calculation.set_batch_grid, use args.balance, args.axes, optional "
        "args.fixed, args.total, args.T, and args.P; axis ranges use "
        "[start, stop, step] and the balance species receives the residual. "
        "For calculation.set_nucleation_undercooling, use args.undercooling as "
        "a phase-to-temperature object. calculation.load_database loads a .dat "
        "or .tdb file from the active calculation directory's database folder. "
        "If no active calculation directory is set, Equilipy will prompt the "
        "user to set one when the load action runs. Equilipy will create the "
        "database folder when the load action runs. If multiple database files "
        "are in the folder, include args.file with the requested file name. To "
        "run the current module after an explicit user request, use "
        'calculation.calculate with args {"confirmed":true}. To run every '
        "module after an explicit user request, use calculation.calculate_all "
        'with args {"confirmed":true}. For all other data/file/calculation '
        "changes, explain what the user should do or ask for confirmation in "
        "the GUI."
    )


def action_response_schema_json() -> str:
    """Return the structured response schema as compact JSON."""
    return json.dumps(ACTION_RESPONSE_SCHEMA, separators=(",", ":"))


def render_action_block(actions: list[object]) -> str:
    """Render action objects into the tagged text bridge used by the GUI."""
    payload = {"actions": actions}
    action_text = json.dumps(payload, separators=(",", ":"))
    return f"{ACTION_BLOCK_START}\n{action_text}\n{ACTION_BLOCK_END}"


def structured_response_to_text(value: object) -> str:
    """Render a structured assistant response into visible text plus actions."""
    if not isinstance(value, dict):
        return ""
    reply = str(value.get("reply") or "").strip()
    actions = value.get("actions")
    if not actions:
        return reply
    if not isinstance(actions, list):
        return reply
    action_block = render_action_block(actions)
    if reply:
        return f"{reply}\n{action_block}"
    return action_block


def strip_json_code_fence(text: str) -> str:
    """Return JSON payload text, tolerating Markdown fences from providers."""
    stripped = str(text or "").strip()
    if not stripped.startswith("```"):
        return stripped
    lines = stripped.splitlines()
    if len(lines) < 2:
        return stripped
    first = lines[0].strip().lower()
    if first not in {"```", "```json", "```javascript", "```js"}:
        return stripped
    if lines[-1].strip() != "```":
        return stripped
    return "\n".join(lines[1:-1]).strip()
