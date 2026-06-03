#!/usr/bin/env python3
"""Convert backtick expressions in manuscript.md that contain unicode
math characters or LaTeX commands into LaTeX inline math `$...$`.
Leave code identifiers (file paths, Rust types) alone.

Math triggers (any one is enough): η, Δ, α, σ, ², ³, ⁰⁴⁵⁶⁷⁸⁹, ⁻, ·,
∞, ≈, ‖, ≤, ≥, ∂, ←, →, − (en-dash math minus). Also: `\\` (latex command).
"""

import re
import sys

MATH_CHARS = "ηΔασ²³⁰⁴⁵⁶⁷⁸⁹⁻·∞≈‖≤≥∂←→−"
MATH_PAT = re.compile(rf"[{MATH_CHARS}]")
LATEX_CMD = re.compile(r"\\\\[A-Za-z]+")

# Unicode replacements that are valid inside math mode but need their
# LaTeX equivalents because the math font lacks the glyph.
UNI_TO_LATEX = {
    "η": r"\eta ",
    "Δ": r"\Delta ",
    "α": r"\alpha ",
    "σ": r"\sigma ",
    "∞": r"\infty ",
    "≈": r"\approx ",
    "‖": r"\|",
    "≤": r"\leq ",
    "≥": r"\geq ",
    "∂": r"\partial ",
    "←": r"\leftarrow ",
    "→": r"\rightarrow ",
    "−": r"-",       # en-dash minus -> ASCII minus
    "·": r"\cdot ",
    "²": r"^{2}",
    "³": r"^{3}",
    "⁰": r"^{0}",
    "⁴": r"^{4}",
    "⁵": r"^{5}",
    "⁶": r"^{6}",
    "⁷": r"^{7}",
    "⁸": r"^{8}",
    "⁹": r"^{9}",
    "⁻": r"^{-}",
}


def unicode_to_latex_math(s: str) -> str:
    # Handle superscript runs: e.g. ⁻¹⁶ → ^{-16} (combine the - and digits).
    s = re.sub(
        r"([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)",
        lambda m: "^{" + "".join(
            {"⁰":"0","¹":"1","²":"2","³":"3","⁴":"4",
             "⁵":"5","⁶":"6","⁷":"7","⁸":"8","⁹":"9","⁻":"-"}[c]
            for c in m.group(1)
        ) + "}",
        s,
    )
    # Handle subscript runs: ₀₁₂...₉ → _{0}_{1}... grouped together.
    s = re.sub(
        r"([₀₁₂₃₄₅₆₇₈₉]+)",
        lambda m: "_{" + "".join(
            {"₀":"0","₁":"1","₂":"2","₃":"3","₄":"4",
             "₅":"5","₆":"6","₇":"7","₈":"8","₉":"9"}[c]
            for c in m.group(1)
        ) + "}",
        s,
    )
    # Norm with infinity subscript: ‖X‖∞ → \|X\|_\infty (must run BEFORE
    # the generic ∞ replacement so we can detect the post-‖ position).
    s = re.sub(r"‖([^‖]+)‖\s*∞", r"\\|\1\\|_\\infty", s)
    # Then the remaining unicode mappings (η, Δ, ≈, etc).
    out = []
    for ch in s:
        out.append(UNI_TO_LATEX.get(ch, ch))
    return "".join(out)


def convert_backticks(text: str) -> str:
    def repl(m: re.Match) -> str:
        content = m.group(1)
        if MATH_PAT.search(content) or LATEX_CMD.search(content):
            # Math expression: convert unicode + escape backslashes properly
            # The markdown source has double backslash for LaTeX commands;
            # pandoc passes them through. For math mode we want single
            # backslash, so unescape.
            content = content.replace(r"\\", "\\")
            content = unicode_to_latex_math(content)
            return f"${content}$"
        return f"`{content}`"

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
