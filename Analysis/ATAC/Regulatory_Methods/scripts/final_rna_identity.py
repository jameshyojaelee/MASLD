"""Source-declared RNA identity and diagnosis, independent of the ATAC crosswalk."""
import json


def source_rna_libraries(candidate):
    workbook = json.loads((candidate/'rna_source_audit/source_workbook.json').read_text())
    result = {}
    for rows in workbook.values():
        for row in rows[1:]:
            if not row[1]:
                continue
            specimen, library, donor, condition, gsm = row[:5]
            assert library not in result and condition in ['NORMAL', 'MASL', 'MASH']
            assert all(isinstance(value, str) and value for value in [specimen, library, donor, gsm])
            result[library] = dict(specimen=specimen, library_id=library, source_donor_id=donor,
                                   condition=condition, gsm=gsm)
    assert len(result) == 18
    assert len({row['source_donor_id'] for row in result.values()}) == 18
    assert len({row['gsm'] for row in result.values()}) == 18
    return result
