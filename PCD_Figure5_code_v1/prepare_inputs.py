"""Generate restricted, coded analysis inputs from the frozen Confirmed worksheet."""
import argparse
import json
import math
from pathlib import Path
from openpyxl import load_workbook

PHENOTYPES = {'Brain disease': 'Brain disease_new',
 'Bronchiectasis': 'bronchiectasis',
 'Bronchitis': 'bronchitis',
 'Chronic respiratory infections': 'chronic_respiratory_infections',
 'Chronic wet cough': 'chronic_wet_cough',
 'Congenital Heart Disease': 'CHD_new',
 'Fetal ventriculomegaly': 'Fetal lateral ventricle widening_new',
 'Hearing impairment': 'Hearing_impairment_new',
 'Infertility': 'infertility_new',
 'Kidney disease': 'Kindey disease_new',
 'Maternal miscarriage': 'maternal miscarriage_new',
 'Nasosinusitis': 'nasosinusitis',
 'Neonatal pulmonary symptoms': 'Neonatal pulmonary symptoms_new',
 'Otitis media mastoiditis': 'Otitis_media/mastoiditis',
 'Pneumonia': 'Pneumonia',
 'Pulmonary atelectasis pulmonary': 'pulmonary atelectasis/pulmonary '
                                    'consolidation_new',
 'Situs inversus': 'Situs_inversus',
 'Skeletal abnormalities': 'Skeletal abnormalities'}
parser = argparse.ArgumentParser()
parser.add_argument('workbook', type=Path)
parser.add_argument('private_output', type=Path)
args = parser.parse_args()
book = load_workbook(args.workbook, read_only=True, data_only=True)
rows = iter(book['Confirmed'].values)
header = next(rows)
index = {name: header.index(name) for name in [*PHENOTYPES.values(), 'Names',
         'Gene_new', 'Enrollment age_new', 'Gender_new', 'EFV1%']}

def number(value):
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None

individuals = {}
for row_id, row in enumerate(rows, start=2):
    if not row[index['Names']]:
        continue
    gene = str(row[index['Gene_new']]).strip()
    record = dict(enrollment_age=number(row[index['Enrollment age_new']]),
                  sex=row[index['Gender_new']], fev=number(row[index['EFV1%']]),
                  phenotypes={p: {'yes': 1, 'no': 0}.get(str(row[index[c]]).strip().lower())
                              for p, c in PHENOTYPES.items()})
    # Confirmed duplicate in this specific frozen workbook. Fail if the row layout changes.
    key = 52 if row_id == 53 else row_id
    if key in individuals:
        previous = individuals[key]
        assert row_id == 53 and set(previous['genes'] + [gene]) == {'DNAH5','DNAH11'}
        assert all(previous[k] == v for k, v in record.items())
        previous['genes'].append(gene)
    else:
        individuals[key] = dict(record, patient_id=f'C{key-1:03d}',
                                genes=[] if gene in ('NA','None','') else [gene])
book.close()
assert len(individuals) == 130 and sum(bool(v['genes']) for v in individuals.values()) == 125
a, b = [], []
for record in individuals.values():
    a.append({k:v for k,v in record.items() if k != 'fev'})
    if record['fev'] is not None and record['enrollment_age'] is not None:
        b.append(dict(patient_id=record['patient_id'],gene=' + '.join(record['genes']),
                      age_at_spirometry_years=record['enrollment_age'],sex=record['sex'],
                      FEV1_percent_predicted=record['fev'],multigene_patient=len(record['genes'])>1))
args.private_output.mkdir(parents=True, exist_ok=True)
for name,records in [('a',a),('b',b)]:
    (args.private_output/(name+'.json')).write_text(json.dumps(records))
