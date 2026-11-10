"""Qualify the native ambient summary before its R process starts."""
from pathlib import Path
import sys


def qualify(text):
    assert "old_q_primary" not in text
    text = text.replace("_70", "_analysis_donors")
    anchor = 'stopifnot(nrow(inc) == 69L)'
    assert anchor in text
    text = text.replace(anchor, anchor + '\nqc <- qc[donor %in% inc$donor]\ndl <- dl[donor %in% inc$donor]')
    text = text.replace('old_q_bh = q_bh, old_q_emp',
                        'old_q_bh = q_bh, old_q_primary = q_primary, old_primary_basis = primary_basis, old_q_emp')
    text = text.replace('new_q_bh = q_bh,\n',
                        'new_q_bh = q_bh, new_q_primary = q_primary, new_primary_basis = primary_basis,\n')
    start = text.index('concl <- function(i)')
    end = text.index('cc <- t(vapply', start)
    body = text[start:end].replace('old_q_bh', 'old_q_primary').replace('new_q_bh', 'new_q_primary')
    body = body.replace('BH q', 'primary FDR q').replace('N3 q', 'raw primary FDR q')
    text = text[:start] + body + text[end:]
    anchor = '  n_q_lt_0.05 = list('
    assert anchor in text
    text = text.replace(anchor,
        '  n_primary_FDR_lt_0.05 = list(old = sum(tab$old_q_primary < 0.05, na.rm = TRUE), new = sum(tab$new_q_primary < 0.05, na.rm = TRUE)),\n'
        '  conclusion_basis = "Prespecified primary FDR, empirical BH when ordinary-BH calibration fails; ordinary BH and contamination-covariate q values are diagnostic only",\n'
        '  n_ordinary_BH_lt_0.05_diagnostic = list(')
    return text


if __name__ == "__main__":
    path = Path(sys.argv[1])
    text = qualify(path.read_text())
    path.write_text(text)
