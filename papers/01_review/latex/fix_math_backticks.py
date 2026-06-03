#!/usr/bin/env python3
"""Convert backtick expressions in manuscript.md to LaTeX inline math
`$...$` when the content is mathematical, and leave code identifiers
as backticks. Idempotent: re-applies safely.

Code patterns (keep backticks):
- File extensions (.py, .rs, .md, .pdf, .R, .tex, .bib, .tif, .csv, .json,
  .toml, .yaml, .yml, .sh, .log)
- Rust generics with angle brackets (<T>, <f64>, etc.)
- Snake_case multi-word identifier (e.g. forward_euler_step) — all
  lowercase ASCII letters + digits + underscores
- CamelCase type identifiers (Conserved2D, Mesh2DG, Dual, Real, etc.)
- The named Rust types we care about explicitly

Math patterns (convert to $...$):
- Anything with unicode math chars (η, Δ, ², ·, ∞, ≈, ‖, ≤, ≥, ∂, ←, →, −)
- LaTeX commands inside (\\X or \X)
- Single-letter math variables (h, g, F, n) or short (dt, dx, dy, hu, hv)
- Tuples of variables ((u, v), (h, hu, hv), (i, j))
- Math subscript like S_{fx}, z_b, h_n
- Numbers with units (1.07 s, 200 m, 6.5×)
- Bare numbers in a backtick (likely a quantitative result)
"""

import re
import sys

MATH_CHARS = "ηΔασ²³⁰⁴⁵⁶⁷⁸⁹⁻·∞≈‖≤≥∂←→−"
MATH_PAT = re.compile(rf"[{MATH_CHARS}]")
LATEX_CMD = re.compile(r"\\\\?[A-Za-z]+")

UNI_TO_LATEX = {
    "η": r"\eta ", "Δ": r"\Delta ", "α": r"\alpha ", "σ": r"\sigma ",
    "∞": r"\infty ", "≈": r"\approx ", "‖": r"\|",
    "≤": r"\leq ", "≥": r"\geq ", "∂": r"\partial ",
    "←": r"\leftarrow ", "→": r"\rightarrow ",
    "−": r"-", "·": r"\cdot ",
}

# Known Rust types / generic patterns that must stay as code.
RUST_KEYWORDS = {
    "f64", "f32", "i32", "i64", "u8", "u32", "u64", "usize", "bool", "str",
    "Real", "Dual", "Conserved2D", "Conserved2DG", "Mesh2D", "Mesh2DG",
    "FluxX", "FluxY", "FluxXG", "FluxYG", "Primitive2D", "Primitive2DG",
    "Boundaries2D", "Boundary", "Side", "Array2", "Vec",
    "manning", "manning_friction_step", "forward_euler_step", "ssprk2_step",
    "cfl_time_step", "max_wave_speeds", "hllc_flux_x", "hllc_flux_y",
    "ghost_cell", "well_balanced_x_face", "well_balanced_y_face",
    "compute_slopes_x", "compute_slopes_y", "build_z_face_x", "build_z_face_y",
    "abs", "max", "min", "sqrt", "powf", "powi", "boundary",
}

FILE_EXT_PAT = re.compile(r"\.[a-zA-Z]{1,5}(?:[\s,;:.)\]}'\"]|$)")
GENERIC_PAT = re.compile(r"<[A-Za-z0-9_]+>")
SNAKE_CASE_PAT = re.compile(r"^[a-z][a-z0-9]*(?:_[a-z0-9]+)+(?:\.\w+)?$")
CAMEL_CASE_PAT = re.compile(r"^[A-Z][A-Za-z0-9]+(?:<[A-Za-z0-9_]+>)?$")
SHELL_CMD_PAT = re.compile(r"^[a-z]+ [a-z]+")  # e.g. "cargo flamegraph"


def is_code(s: str) -> bool:
    s_strip = s.strip()
    # Empty or whitespace only
    if not s_strip:
        return True
    # Rust attribute / lint suppression: #![forbid(...)], #[derive(...)]
    if s_strip.startswith(("#[", "#!")):
        return True
    # File path or filename with extension
    if FILE_EXT_PAT.search(s_strip):
        return True
    # Generic type parameter
    if GENERIC_PAT.search(s_strip):
        return True
    # Known Rust identifier
    if s_strip in RUST_KEYWORDS:
        return True
    # Snake_case identifier
    if SNAKE_CASE_PAT.match(s_strip):
        return True
    # CamelCase type
    if CAMEL_CASE_PAT.match(s_strip):
        return True
    # Shell command like "cargo flamegraph"
    if SHELL_CMD_PAT.match(s_strip) and "=" not in s_strip and "·" not in s_strip:
        return True
    # Anything with the curly-brace syntax of Rust types like Dual {val, dval}
    if re.match(r"^[A-Z][A-Za-z0-9]+ \{", s_strip):
        return True
    return False


def superscript_run(m: re.Match) -> str:
    return "^{" + "".join(
        {"⁰":"0","¹":"1","²":"2","³":"3","⁴":"4",
         "⁵":"5","⁶":"6","⁷":"7","⁸":"8","⁹":"9","⁻":"-"}[c]
        for c in m.group(1)
    ) + "}"


def subscript_run(m: re.Match) -> str:
    return "_{" + "".join(
        {"₀":"0","₁":"1","₂":"2","₃":"3","₄":"4",
         "₅":"5","₆":"6","₇":"7","₈":"8","₉":"9"}[c]
        for c in m.group(1)
    ) + "}"


def to_math(s: str) -> str:
    # Unescape doubled backslashes (markdown escape) so \\partial -> \partial.
    s = re.sub(r"\\\\([A-Za-z]+)", r"\\\1", s)
    s = s.replace("\\\\!", "\\!").replace("\\\\,", "\\,")
    s = re.sub(r"([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", superscript_run, s)
    s = re.sub(r"([₀₁₂₃₄₅₆₇₈₉]+)", subscript_run, s)
    # Norm with infinity subscript: ‖X‖∞ → \|X\|_\infty
    s = re.sub(r"‖([^‖]+)‖\s*∞", r"\\|\1\\|_\\infty", s)
    out = []
    for ch in s:
        out.append(UNI_TO_LATEX.get(ch, ch))
    return "".join(out)


def convert_backticks(text: str) -> str:
    def repl(m: re.Match) -> str:
        content = m.group(1)
        if is_code(content):
            return f"`{content}`"
        return f"${to_math(content)}$"

    return re.sub(r"`([^`\n]+)`", repl, text)


def main():
    path = sys.argv[1]
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    new = convert_backticks(text)
    with open(path, "w", encoding="utf-8") as f:
        f.write(new)


if __name__ == "__main__":
    main()
