# Lab Formulation Calculator

Molar ratios in -> scalable gram charges, condensate prediction, theoretical
end-group values (AV/AmV/OHV), and Carothers/stoichiometry feedback out.

## Features
- Components by CAS or name; MW auto-fills from PubChem (cached). Manual MW
  entry for UVCBs -- use effective MW, or equivalent weight x functionality.
- Assay/purity correction: "active g" vs "as-is g to weigh".
- Scale anchors: total batch mass, fixed component mass ("I have 247 g of
  dimer acid left"), or fixed component moles.
- Condensation accounting for amidation / esterification: predicted water,
  theoretical final mass and % solids, at any conversion p.
- Theoretical acid value / amine value / OH value of the resin -- titration
  targets for "did this batch behave".
- Carothers feedback: r, Xn, Mn, end-group character, and a gel warning at
  near-perfect stoichiometry + high conversion.
- Printable charge sheet with header block (project / experiment ID /
  chemist), add-checkboxes, and blanks to log actual condensate + titrations.

## Deploy (same as the TSCA tool)
Repo containing: formulation_app.py, formulation_core.py, requirements.txt
-> share.streamlit.io -> New app -> main file: formulation_app.py
(or Railway: start command
 `streamlit run formulation_app.py --server.port $PORT --server.address 0.0.0.0`)

## Notes
- Functionality is per the reaction being run (DETA ~2 primary amines for
  amidation; count secondary NH if your chemistry uses it -- enter 5 N-H eq
  if that's how you formulate).
- Theoretical values assume ideal step growth; real cooks deviate. Planning
  tool, not a CoA.
