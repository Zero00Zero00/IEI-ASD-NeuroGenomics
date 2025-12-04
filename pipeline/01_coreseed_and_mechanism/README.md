# 01_coreseed_and_mechanism

Scripts and rules to build the IEI–ASD CoreSeed panel and mechanism annotations.

This module builds the IEI–ASD CoreSeed panel and its basic annotations.

Steps:
- Harmonise gene symbols across IEI catalogues (IUIS and related lists) and ASD resources (e.g. SFARI) against HGNC.
- Annotate each candidate gene with clinical mechanism labels and constraint metrics (e.g. LOEUF, pLI, missense pathogenicity scores).
- Apply the rules described in the manuscript to select the final IEI–ASD CoreSeed gene set.
- Export unified tables that are later used by downstream modules (pathway analysis, axis definition, Set A construction).

This module corresponds to the “CoreSeed definition and mechanistic curation” part of the Supplementary Methods.
