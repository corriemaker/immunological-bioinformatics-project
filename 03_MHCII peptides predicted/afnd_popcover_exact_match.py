#!/usr/bin/env python3
"""Fetch AFND HLA allele frequencies and estimate PopCover inputs.

Requires Python 3 and wget. All queried molecules are fixed in QUERY_ALLELES.
Uses the verified AFND size= parameter and checks all records were returned.
If a complete response is unavailable, FAIL instead of using partial results.
"""
import json
import math
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from urllib.parse import urlencode

QUERY_ALLELES = [
    'DRB1_0101', 'DRB1_0301', 'DRB1_0401',
    'HLA-DPA10103-DPB10201', 'HLA-DPA10103-DPB10301',
    'HLA-DPA10103-DPB10401', 'HLA-DQA10102-DQB10501',
    'HLA-DQA10103-DQB10301', 'HLA-DQA10401-DQB10301',
]
BASE = 'https://www.allelefrequencies.net/api/hla/allele-frequencies/search'
DATA_DIR = Path('afnd_downloads')
OUTPUT_FILE = Path('popcover_allele_frequencies.txt')
REPORT_FILE = Path('afnd_weighted_summary.tsv')
CHAIN_ID = '4G3Y_2_Chain'  # Change if required.
POPULATION_STANDARD = None  # or 'goldOnly'
REQUEST_DELAY = 0.25


def components_for(query):
    return [query] if query.startswith('DRB1_') else query.removeprefix('HLA-').split('-')


def api_name(component):
    m = re.fullmatch(r'(DRB1_|DPA1|DPB1|DQA1|DQB1)(\d\d)(\d\d)', component)
    if not m:
        raise ValueError(f'Invalid component: {component}')
    return f"{m[1].rstrip('_')}*{m[2]}:{m[3]}"


def request(allele, extra=None, allow_failure=False):
    params = {'allele': allele}
    if POPULATION_STANDARD:
        params['populationStandard'] = POPULATION_STANDARD
    if extra:
        params.update(extra)
    url = f'{BASE}?{urlencode(params)}'
    with tempfile.TemporaryDirectory() as d:
        outfile = Path(d) / 'response.json'
        proc = subprocess.run(
            ['wget', '-q', '--server-response', '--tries=1', '--timeout=25',
             '-O', str(outfile), url], capture_output=True, text=True,
        )
        time.sleep(REQUEST_DELAY)
        if proc.returncode != 0:
            if allow_failure:
                return None
            raise RuntimeError(f'wget failed for {allele} (exit={proc.returncode})\n'
                               f'URL: {url}\n{proc.stderr[-900:]}')
        try:
            data = json.loads(outfile.read_text())
        except (json.JSONDecodeError, OSError):
            if allow_failure:
                return None
            raise RuntimeError(f'Invalid JSON from {url}')
        if not isinstance(data, dict) or not isinstance(data.get('content'), list):
            if allow_failure:
                return None
            raise RuntimeError(f'Unexpected JSON format from {url}')
        return data


def fetch_all(allele):
    # AFND has been verified by the user to honour size=1000 and return
    # all 587 DRB1*01:01 entries in one response.
    data = request(allele, {'size': 1000})
    total = int(data.get('totalElements', -1))
    if total < 0:
        raise RuntimeError(f'Missing totalElements for {allele}')

    # If the allele has >1000 records, request one page big enough for all.
    if total > 1000:
        print(f'  {total} records reported; requesting a larger page ...', flush=True)
        data = request(allele, {'size': total})

    rows = data['content']
    if (int(data.get('totalElements', -1)) != total
            or int(data.get('totalPages', -1)) != 1
            or int(data.get('currentPage', -1)) != 1
            or len(rows) != total):
        raise RuntimeError(
            f'Incomplete AFND response for {allele}: downloaded {len(rows)} '
            f'of {total} records (totalPages={data.get("totalPages")}). '
            'Not calculating a frequency from partial results.'
        )
    # Some AFND search responses may include additional allele names. Keep only
    # exact two-field matches; report what was excluded for auditability.
    matching = [r for r in rows if r.get('allele') == allele]
    other = sorted({str(r.get('allele')) for r in rows if r.get('allele') != allele})
    if other:
        print(f'  WARNING: AFND returned {len(rows) - len(matching)} records with '
              f'non-matching allele names; excluded.', flush=True)
        print('  Other allele names: ' + ', '.join(other[:12]) +
              (' ...' if len(other) > 12 else ''), flush=True)
    if not matching:
        raise RuntimeError(f'No exact allele matches for {allele}; '
                           'inspect the downloaded JSON and query conventions')

    safe = allele.replace('*', '_').replace(':', '_')
    DATA_DIR.mkdir(exist_ok=True)
    path = DATA_DIR / f'{safe}.json'
    path.write_text(json.dumps(data, indent=2) + '\n')
    return matching, path, len(rows)


def weighted_frequency(records, allele):
    numerator = 0.0
    denominator = 0
    valid = 0
    for record in records:
        try:
            f = float(record['geneFrequency'])
            n = int(record['sampleSize'])
        except (KeyError, ValueError, TypeError):
            continue
        if not math.isfinite(f) or not 0 <= f <= 1 or n <= 0:
            continue
        numerator += f * n
        denominator += n
        valid += 1
    if denominator == 0:
        raise RuntimeError(f'No usable data for {allele}')
    return numerator / denominator, numerator, denominator, valid


def main():
    components = list(dict.fromkeys(p for q in QUERY_ALLELES for p in components_for(q)))
    print(f'{len(QUERY_ALLELES)} MHC-II molecules; {len(components)} unique component alleles', flush=True)
    frequencies = {}
    report = ['component\tallele\tfrequency\tweighted_sum\tsample_size_sum\tvalid_records\tall_records\tjson_file']
    for comp in components:
        allele = api_name(comp)
        print(f'Fetching {allele} ...', flush=True)
        records, path, total_records = fetch_all(allele)
        f, numer, denom, valid = weighted_frequency(records, allele)
        frequencies[comp] = f
        report.append(f'{comp}\t{allele}\t{f:.8f}\t{numer:.6f}\t{denom}\t{valid}\t{total_records}\t{path}')
        print(f'  {total_records} downloaded, {len(records)} exact matches, {valid} usable, frequency={f:.6f}', flush=True)
    final_lines = []
    for query in QUERY_ALLELES:
        freq = math.prod(frequencies[c] for c in components_for(query))
        locus = ('DRB1' if query.startswith('DRB1_') else
                 'HLA-DP' if 'DPA1' in query else 'HLA-DQ')
        final_lines.append(f'{query} {locus} {CHAIN_ID} {freq:.6f}')
    REPORT_FILE.write_text('\n'.join(report) + '\n')
    OUTPUT_FILE.write_text('\n'.join(final_lines) + '\n')
    print(f'\nSaved {len(components)} allele JSON files in {DATA_DIR}/')
    print(f'Saved: {REPORT_FILE}\nSaved: {OUTPUT_FILE}')
    print('\n'.join(final_lines))


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError) as exc:
        sys.exit(f'ERROR: {exc}')
