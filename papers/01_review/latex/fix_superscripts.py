#!/usr/bin/env python3
"""Convert unicode superscript runs in a LaTeX file into math-mode
superscripts so that Latin Modern (or any standard math font) can
render them. Operates in place on the file passed as argument.

Pattern: a single base character (letter or digit) followed by 1+
unicode superscript chars (⁰⁻¹²³⁴⁵⁶⁷⁸⁹) becomes
`base$^{ASCII_run}$`. Example: `L¹` → `L$^{1}$`, `10⁻¹⁶` → `10$^{-16}$`.

Run:
  python3 fix_superscripts.py body_clean.tex
"""

import re
import sys

SUP = {
    "⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4",
    "⁵": "5", "⁶": "6", "⁷": "7", "⁸": "8", "⁹": "9",
    "⁻": "-",
}

PAT = re.compile(r"([A-Za-z0-9])([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)")


def convert(m: re.Match) -> str:
    run = "".join(SUP[c] for c in m.group(2))
    return f"{m.group(1)}$^{{{run}}}$"


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: fix_superscripts.py <file>", file=sys.stderr)
        sys.exit(2)
    path = sys.argv[1]
    with open(path, "r", encoding="utf-8") as f:
        text = f.read()
    new = PAT.sub(convert, text)
    with open(path, "w", encoding="utf-8") as f:
        f.write(new)


if __name__ == "__main__":
    main()
