# Talbot emulation of the witness frames — code, result files and manuscript source

This record accompanies the manuscript "Emulating the frames of a mutual-information witness with time-bin light,
and how precise the optics must be". It describes only what is in the record.

> **Important — licensing:** text and data under CC BY 4.0; code under the Apache License, Version 2.0. Please read [IMPORTANT.md](IMPORTANT.md).

## Contents

| Folder / file | What it is |
|---|---|
| `manuscript/` | LaTeX source of the manuscript (`talbot_manuscript.tex`, `talbot_refs.bib`) and Fig. 1 (`a_fig1_hero.pdf`, written by `code/a_fig1_hero.py`). |
| `code/` | Python scripts, the result file of run 1, and the console outputs of the check scripts (table below). |
| `LICENSE` | Apache License, Version 2.0 (applies to the code). |
| `CHECKSUMS.sha256` | SHA-256 of every file in this record except itself (text files with LF line endings, except the two run-1 files, which are kept byte for byte as written by the run and have CRLF; all checksums are over the deposited bytes), checkable with `sha256sum -c`. |

| File in `code/` | What it is |
|---|---|
| `a_kern.py` | Building blocks: SU(2)_k modular data, the map from labels to bins, frames, entropies. |
| `a_lauf.py` | The computation functions for the tests A1 and A2 and the fixed settings `PREREG`. |
| `a_tests.py` | The pre-registered run that produced run 1. It first verifies the locked pre-registration (see below). |
| `a_reproduzieren.py` | The same computation without the lock check; compares its result with run 1. |
| `a_tests_results_LAUF1_2026-10-03.json`, `a_tests_out_LAUF1_2026-10-03.txt` | Run 1: result file and console output. |
| `A_PREREG_LOCK.sha256` | SHA-256 of the locked pre-registration (LF form). |
| `a_auswertung.py` | Evaluates run 1 against the pre-registered criteria. |
| `a_quantil_exakt.py` | Sampling-free quantile of the phase-noise margin in frame `ST`. |
| `a_nachtrag.py` | Addendum, not part of the pre-registration: the sampling-free quantile on the whole sigma grid (80 values) and the sigma = 0 control with 40 000 repetitions. |
| `a_vergleich_zweite_route.py` | Comparison of run 1 with the values of an independently written second implementation. The values are transcribed in the script; their source is the next file. Its output reports one point outside the bounds set in advance (k = 2, sigma = 0.02: 0.0208 against 0.0213, 2.4 % against a bound of 2 %); the verdicts are the same. |
| `a2_second_route.py`, `a2_second_route_result.json`, `a2_second_route_out.txt` | The independently written second implementation of hypothesis A2 (semi-analytic for (b) and (c), separate Monte Carlo for (a)), its result file and the console output of a re-run. Needs numpy and scipy, reads no input file; run `python a2_second_route.py [output_directory]` (about 40 s). Its comments and printed text are English; the keys of the result file are kept as in the reference run (partly German) so that both result files compare entry by entry. The "sekunden" entry (running time) differs from run to run; all other entries, including the entry `k`, were reproduced exactly in a re-run. It is a second route to the same verdicts, not a result that preceded the lock of the pre-registration. |
| `a_nk.py` | Negative and positive controls of the computation tools (not a pre-registered test). Its output ends with "Kontrollen: NICHT bestanden" (controls not passed): the criterion set beforehand expected the wrong (even) embedding to break both identities at k = 2 and 4, but at k = 2 the entropy identity stays intact. This is the result reported in the manuscript (section "Mapping"), and the reason why the test A1 includes the matrix identities. |
| `a_vek_check.py` | Vectorised against scalar building blocks. |
| `a_paare.py` | Blind labels and blind product pairs for k = 2 and 4. |
| `a_paare_wortlaenge.py` | Shortest frame word that shows a non-blind pair as non-integer. |
| `a_tabelle_woerter.py` | Values for the table of frames in the manuscript. |
| `a_fig1_hero.py` | Draws Fig. 1 of the manuscript: panel (a) is a direct calculation for one example, panels (b) and (c) use run 1; checks the example identity and the thresholds against run 1. Needs matplotlib. |
| `a_zahlencheck_manuskript.py` | Checks the printed result numbers listed in the script against run 1 and the outputs, and that each appears in the manuscript. |
| `*_out.txt` | Console output of the script of the same name, produced in this folder. |

## Requirements and running

Tested with Python 3.12.10 and numpy 2.4.3 (Windows 11), the environment recorded in run 1; the figure script also needs
matplotlib (tested with 3.10.8). No other packages.

1. In `code/`, run `python a_reproduzieren.py`. It writes `a_tests_results_reproduziert.json` and compares it with run 1.
   In the environment above the result is identical to run 1 (all of A2; the A1 deviations bit for bit), see
   `a_reproduzieren_out.txt`. Running time about 2 minutes on one machine (measured: 106 s).
2. The check scripts read only run 1 or compute on their own; run each as `python <script>.py > <script>_out.txt`.
   `a_zahlencheck_manuskript.py` reads `a_nk_out.txt` and `a_nachtrag_out.txt`, so run `a_nk.py` and `a_nachtrag.py` first.
   The random seeds are fixed; the numpy random generator does not guarantee the same stream across versions or platforms.

## Pre-registration

The tests A1 and A2 were pre-registered, and the pre-registration was locked before the computation; its SHA-256 is in
`A_PREREG_LOCK.sha256` and in the manuscript. The pre-registration file itself is not part of this record. `a_tests.py` verifies the locked file before computing and therefore
stops in this version with a missing-file error (`A_PREREG_ENTWURF_v4.md`), before computing anything;
`a_reproduzieren.py` performs the same computation without that check.

## Changes for publication

The scripts are the files used for run 1 and for the checks. For publication, internal working names were removed from
comments and printed labels, file names and paths were adjusted, and a comparison printout was removed from
`a_paare.py`; no numerical code was changed. `a_reproduzieren.py` was written for this record. `a2_second_route.py` is a copy of the second implementation in which only the surroundings were changed (English comments and printed text, no file access except the result, output directory as argument); the calculation and the seeds are the same, and the entry `k` of its result file agrees exactly with the reference run. Comments in the other scripts are German with an English summary in the header; a complete English version follows in a later version.

## Licence and citation

Code: Apache License, Version 2.0 (file `LICENSE`; see also `NOTICE` and `IMPORTANT.md`). Text and data: CC BY 4.0.
Cite as: Sayim, B. Y. (2026), version 1.0, DOI 10.5281/zenodo.23148214.
