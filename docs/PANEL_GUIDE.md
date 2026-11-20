# Panel guide

The keys below belong to the original data archive. Use the [manuscript crosswalk](manuscript_panel_crosswalk.csv) and [final panel guide](MANUSCRIPT_PANELS.md) when looking up the current figure numbers.

See [Methods-to-code map](METHODS_CODE_MAP.md), [panel input counts](PANEL_INPUTS.csv) and the [Figure 1j-l guide](FIG01J_L.md). `python reproduce.py --list-panels` lists all 43 configured workflows.

| Panel/workflow | Companion panels | Starting point | Archived P checks | Configuration |
|---|---|---|---:|---|
| Fig01b-d |  | pixel exports from all eight settings; complete matched raw GSLab input sets | 0 | [JSON](../configs/Fig01b-d.json) |
| Fig01f-h |  | exact per-image object exports with retained row/ID selection | 6 | [JSON](../configs/Fig01f-h.json) |
| Fig01j-l |  | Selected object tables and uniquely recovered anchor IDs; archived projections and SVG curves verified | 0 | [JSON](../configs/Fig01j-l.json) |
| Fig02d_cells | Fig02d | exact per-image object exports with retained row/ID selection | 9 | [JSON](../configs/Fig02d_cells.json) |
| Fig02d_nuclei | Fig02d | deposited object or measurement exports | 9 | [JSON](../configs/Fig02d_nuclei.json) |
| Fig03b |  | deposited pooled-bin summary | 0 | [JSON](../configs/Fig03b.json) |
| Fig03c | FigS07d-g | deposited pixel exports | 36 | [JSON](../configs/Fig03c.json) |
| Fig03d-h |  | deposited pixel exports | 52 | [JSON](../configs/Fig03d-h.json) |
| Fig03i-m |  | deposited pixel exports | 78 | [JSON](../configs/Fig03i-m.json) |
| Fig03m_lumen |  | deposited pixel exports | 4 | [JSON](../configs/Fig03m_lumen.json) |
| Fig03n-p | FigS08m-o | deposited object or measurement exports | 2 | [JSON](../configs/Fig03n-p.json) |
| Fig04c |  | deposited pixel exports | 3 | [JSON](../configs/Fig04c.json) |
| Fig04d | FigS09 | deposited pixel exports | 10 | [JSON](../configs/Fig04d.json) |
| Fig04g-h |  | deposited pixel exports | 18 | [JSON](../configs/Fig04g-h.json) |
| Fig05b |  | deposited pixels and final measurement masks | 6 | [JSON](../configs/Fig05b.json) |
| Fig05d |  | deposited pixels and final measurement masks | 6 | [JSON](../configs/Fig05d.json) |
| Fig05f |  | deposited pixels and final measurement masks | 9 | [JSON](../configs/Fig05f.json) |
| Fig05h |  | deposited pixels and final measurement masks | 9 | [JSON](../configs/Fig05h.json) |
| Fig05j |  | deposited pixels and final measurement masks | 9 | [JSON](../configs/Fig05j.json) |
| Fig06b |  | deposited pixel exports | 18 | [JSON](../configs/Fig06b.json) |
| Fig06d |  | deposited pixel exports | 3 | [JSON](../configs/Fig06d.json) |
| Fig06g |  | deposited object or measurement exports | 13 | [JSON](../configs/Fig06g.json) |
| FigS01b-c |  | deposited pixel exports | 6 | [JSON](../configs/FigS01b-c.json) |
| FigS01d |  | deposited object or measurement exports | 0 | [JSON](../configs/FigS01d.json) |
| FigS02 |  | exact per-image object exports with retained row/ID selection | 15 | [JSON](../configs/FigS02.json) |
| FigS02_ROI |  | deposited object or measurement exports | 0 | [JSON](../configs/FigS02_ROI.json) |
| FigS03 |  | deposited object or measurement exports | 0 | [JSON](../configs/FigS03.json) |
| FigS05c-d |  | deposited pixel exports | 18 | [JSON](../configs/FigS05c-d.json) |
| FigS07b |  | deposited source-sample module scores | 0 | [JSON](../configs/FigS07b.json) |
| FigS07c |  | deposited source-sample module scores | 0 | [JSON](../configs/FigS07c.json) |
| FigS07j-k |  | deposited pixel exports | 3 | [JSON](../configs/FigS07j-k.json) |
| FigS07l-n |  | deposited pixel exports | 2 | [JSON](../configs/FigS07l-n.json) |
| FigS08a-e |  | deposited pixel exports | 52 | [JSON](../configs/FigS08a-e.json) |
| FigS08g-i |  | deposited time-resolved image summary | 0 | [JSON](../configs/FigS08g-i.json) |
| FigS08j-l |  | deposited time-resolved image summary | 0 | [JSON](../configs/FigS08j-l.json) |
| FigS10b | FigS10a | deposited pixel exports | 4 | [JSON](../configs/FigS10b.json) |
| FigS10c-h |  | deposited pixel exports | 24 | [JSON](../configs/FigS10c-h.json) |
| FigS10i-l |  | deposited pixel exports | 16 | [JSON](../configs/FigS10i-l.json) |
| FigS11b |  | deposited pixels and final measurement masks | 9 | [JSON](../configs/FigS11b.json) |
| FigS11c |  | deposited pixels and final measurement masks | 9 | [JSON](../configs/FigS11c.json) |
| FigS11d |  | deposited pixels and final measurement masks | 9 | [JSON](../configs/FigS11d.json) |
| FigS12b |  | deposited pixel exports | 24 | [JSON](../configs/FigS12b.json) |
| FigS12d |  | deposited pixel exports | 4 | [JSON](../configs/FigS12d.json) |
