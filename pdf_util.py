"""Charge-sheet PDF export. Renders the existing plain-text charge sheets
(charge_sheet_text / blend_sheet_text / qa_solids_sheet_text) into a monospaced
PDF — a render/export utility, no chemistry. Uses fpdf2 (lightweight, pure
Python; add `fpdf2` to requirements.txt)."""
from fpdf import FPDF


def text_to_pdf(text, title="Charge sheet"):
    """Return PDF bytes of a fixed-width text block (Courier), letter size.
    The charge-sheet text is ASCII/Latin-1; any stray non-Latin-1 char is
    replaced rather than raising. Uses cell() per line (no wrapping) — the
    charge sheets are pre-formatted fixed-width and fit the page."""
    pdf = FPDF(format="letter", unit="mm")
    pdf.set_title(title)
    pdf.set_auto_page_break(auto=True, margin=12)
    pdf.set_margins(12, 12, 12)
    pdf.add_page()
    pdf.set_font("Courier", size=8)
    for raw in text.split("\n"):
        line = raw.encode("latin-1", "replace").decode("latin-1")
        pdf.cell(0, 4.0, line, new_x="LMARGIN", new_y="NEXT")
    return bytes(pdf.output())   # fpdf2 returns a bytearray
