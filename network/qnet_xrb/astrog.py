#!/usr/bin/env python3
"""Replace theoretical rates with SMARAGD fits in the existing XRB format.

Based on the direction-aware search in astrog_mod.py (ntwk_gen 8ddd698).
SMARAGD's numeric Z,N fields, never its product-name field, identify reactions.
Reverse fits omit partition functions: the XRB solver applies them when the
last integer (indinv) is 1. Original Q values and partition data are retained.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import io
import json
import math
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

# Preserve the original user-editable defaults.
Rate_mod = True
Ra26 = True
Ra26_all = False
n_selec = 0

THEORY_LABELS = frozenset(('rath', 'ths8', 'thra'))
PARTICLES = {'g': (0, 0), 'n': (0, 1), 'p': (1, 0), 'a': (2, 2)}
CHANNELS = frozenset(('pg', 'gp', 'ng', 'gn', 'ag', 'ga',
                      'pn', 'np', 'pa', 'ap', 'na', 'an'))
MODE_CHANNELS = {2: {'pg', 'gp'}, 3: {'ap', 'pa', 'ag', 'ga'},
                 4: {'ng', 'gn', 'np', 'pn'}}
ELEMENTS = ('x h he li be b c n o f ne na mg al si p s cl ar k ca sc ti v cr mn fe co ni cu zn '
            'ga ge as se br kr rb sr y zr nb mo tc ru rh pd ag cd in sn sb te i xe cs ba la ce '
            'pr nd pm sm eu gd tb dy ho er tm yb lu hf ta w re os ir pt au hg tl pb bi po at '
            'rn fr ra ac th pa u np pu am cm bk cf es fm md no lr rf db sg bh hs mt ds rg '
            'cn nh fl mc lv ts og').split()
ZN = tuple[int, int]
Reaction = tuple[ZN, ZN, ZN, ZN]  # target, projectile, ejectile, product


class RateError(ValueError):
    """An input or ambiguous match that must stop conversion before output."""


def require(condition, message):
    if not condition:
        raise RateError(message)


def parse_name(name: str) -> ZN:
    name = name.strip().lower()
    special = {'gam': (0, 0), 'n': (0, 1), 'p': (1, 0), 'd': (1, 1), 't': (1, 2)}
    if name in special:
        return special[name]
    match = re.fullmatch(r'([a-z]+)(\d+)', name)
    require(match is not None and match[1] in ELEMENTS, f'Unknown nuclide: {name!r}')
    z = ELEMENTS.index(match[1])
    n = int(match[2]) - z
    require(z > 0 and n >= 0, f'Invalid nuclide: {name!r}')
    return z, n


def format_name(zn: ZN) -> str:
    if zn == (0, 0):
        return 'gam'
    if zn == (0, 1):
        return 'n'
    return f'{ELEMENTS[zn[0]]}{sum(zn)}'


def numbers(text: str, count: int, where: str) -> tuple[float, ...]:
    try:
        values = tuple(float(v.replace('D', 'E').replace('d', 'e')) for v in text.split())
    except ValueError as exc:
        raise RateError(f'{where}: invalid number') from exc
    require(len(values) == count and all(math.isfinite(v) for v in values),
            f'{where}: expected {count} finite numbers, got {len(values)}')
    return values


@dataclass(frozen=True)
class Nuclide:
    name: str
    a: int
    zn: ZN
    mass_excess: float


def read_partition(path: Path) -> list[Nuclide]:
    lines = path.read_text().splitlines()
    count = int(lines[0])
    require(len(lines) == 1 + 4*count, f'{path}: expected 4 lines per nuclide')
    nuclides = []
    for offset in range(1, len(lines), 4):
        row = lines[offset].split()
        require(len(row) == 6, f'{path}:{offset+1}: invalid nuclide row')
        name, a, z, n = row[:4]
        zn = int(z), int(n)
        require(parse_name(name) == zn and int(a) == sum(zn),
                f'{path}:{offset+1}: nuclide name, A, Z, N disagree')
        pf = numbers(' '.join(lines[offset+1:offset+4]), 24, f'{path}:{offset+2}')
        require(all(v > 0 for v in pf), f'{path}:{offset+2}: partition function must be positive')
        nuclides.append(Nuclide(name, int(a), zn, float(row[5])))
    require(len({n.zn for n in nuclides}) == count, f'{path}: duplicate nuclide')
    return nuclides


@dataclass
class RateRecord:
    section: str
    line: int
    header: str
    ids: tuple[int, ...]
    raw: list[str]
    terms: list[tuple[float, ...]] = field(default_factory=list)
    q: float | None = None
    iaa: bool = False

    @property
    def label(self):
        return self.header[40:45].strip()


def read_rates(path: Path) -> list[str | RateRecord]:
    # Retain raw bytes (including line endings) for every untouched record.
    lines = path.read_bytes().decode('ascii').splitlines(keepends=True)
    items, section, i, ended = [], None, 0, False
    while i < len(lines):
        line = lines[i].rstrip('\r\n')
        words = line.split()
        if not words or line.startswith('#'):
            items.append(lines[i]); i += 1; continue
        if words[:2] == ['end', 'of']:
            items.extend(lines[i:]); ended = True; break
        if words[:2] == ['data', 'of']:
            section = words[2]
            items.append(lines[i]); i += 1; continue
        require(section is not None, f'{path}:{i+1}: rate outside reaction section')
        iaa = line[40:45].strip() == 'iaa'
        width = 6 if iaa else 4
        require(i+width <= len(lines), f'{path}:{i+1}: incomplete rate block')
        ids = tuple(map(int, lines[i+1].split()))
        require(len(ids) == 8, f'{path}:{i+2}: expected 8 integers')
        values = numbers(' '.join(lines[i+2:i+width]), 28 if iaa else 8, f'{path}:{i+3}')
        record = RateRecord(section, i+1, line, ids, lines[i:i+width],
                            [] if iaa else [values[:7]], None if iaa else values[7], iaa)
        # Consecutive identical headers contain additional fits to one rate.
        previous = items[-1] if items else None
        if (not iaa and isinstance(previous, RateRecord) and not previous.iaa
                and previous.section == section and previous.header == line and previous.ids == ids):
            require(previous.q == record.q, f'{path}:{i+1}: Q differs between fit terms')
            previous.raw.extend(record.raw)
            previous.terms.extend(record.terms)
        else:
            items.append(record)
        i += width
    require(ended, f'{path}: missing end marker')
    return items


@dataclass(frozen=True)
class Fit:
    identity: Reaction
    coefficients: tuple[float, ...]
    direction: str
    source: str
    line: int
    q: float  # Q for this oriented reaction
    original_names: str


def read_smaragd(directory: Path):
    candidates, mismatches = {}, []
    for projectile, particle in [('p', 'p'), ('n', 'n'), ('a', 'a')]:
        path = directory/f'fits_{projectile}.dat'
        lines = path.read_text().splitlines()
        require(len(lines) % 3 == 0, f'{path}: incomplete 3-line fit entry')
        for i in range(0, len(lines), 3):
            row = lines[i].split()
            require(len(row) == 15, f'{path}:{i+1}: expected 15 header fields')
            key = tuple(tuple(map(int, row[j:j+2])) for j in (4, 6, 8, 10))
            require(all(0 <= z < len(ELEMENTS) and n >= 0 for z, n in key),
                    f'{path}:{i+1}: invalid numeric Z,N')
            require(key[1] == PARTICLES[particle] and key[2] in PARTICLES.values(),
                    f'{path}:{i+1}: wrong projectile or ejectile')
            require(tuple(key[0][j]+key[1][j] for j in (0, 1)) ==
                    tuple(key[2][j]+key[3][j] for j in (0, 1)),
                    f'{path}:{i+1}: numeric Z,N are not conserved')
            for j, (label, zn) in enumerate(zip(row[:4], key)):
                if parse_name(label) != zn:
                    mismatches.append(dict(source=path.name, line=i+1, field=j,
                                           supplied_name=label, numeric_name=format_name(zn), z=zn[0], n=zn[1]))
            q = numbers(row[-1], 1, f'{path}:{i+1}')[0]
            for direction, oriented, sign, offset in [
                ('direct', key, 1, 1), ('reverse', (key[3], key[2], key[1], key[0]), -1, 2)
            ]:
                coeff = numbers(lines[i+offset], 7, f'{path}:{i+offset+1}')
                candidate = Fit(oriented, coeff, direction, path.name, i+1, sign*q, ' '.join(row[:4]))
                previous = candidates.get(oriented)
                require(previous is None,
                        f'{path}:{i+1}: duplicate reaction {oriented}; '
                        f'previous={previous.source+":"+str(previous.line) if previous else ""}')
                candidates[oriented] = candidate
    return candidates, mismatches


def reaction_identity(record: RateRecord, nuclides: list[Nuclide]) -> Reaction:
    slots = []
    for j, number in enumerate(record.ids[:6]):
        label = record.header[j*5:j*5+5].strip()
        if number == 0:
            require(not label, f'input:{record.line}: name without nuclide index')
            slots.append(None)
        else:
            require(0 < number <= len(nuclides), f'input:{record.line}: nuclide index out of range')
            zn = nuclides[number-1].zn
            require(parse_name(label) == zn, f'input:{record.line}: name/index disagree: {label}, {number}')
            slots.append(zn)
    rtype = record.section
    if rtype.startswith('g'):
        key = slots[0], PARTICLES['g'], slots[1], slots[2]
        unused = slots[3:]
    elif rtype.endswith('g'):
        key = slots[1], slots[0], PARTICLES['g'], slots[2]
        unused = slots[3:]
    else:
        key = slots[1], slots[0], slots[2], slots[3]
        unused = slots[4:]
    require(all(zn is not None for zn in key) and not any(unused),
            f'input:{record.line}: wrong reaction multiplicity')
    require(key[1] == PARTICLES[rtype[0]] and key[2] == PARTICLES[rtype[1]],
            f'input:{record.line}: particles disagree with section {rtype}')
    require(tuple(key[0][j]+key[1][j] for j in (0, 1)) ==
            tuple(key[2][j]+key[3][j] for j in (0, 1)),
            f'input:{record.line}: Z,N are not conserved')
    return key


def selected(record: RateRecord, key: Reaction, all_rates: bool, mode: int) -> bool:
    if all_rates:  # Retain the original Ra26_all override semantics.
        return True
    if record.label not in THEORY_LABELS:
        return False
    if mode == 0:
        return True
    if mode == 1:
        return key[0][0] < 47 and key[0][1] < 53
    return record.section in MODE_CHANNELS[mode]


def render_replacement(record: RateRecord, fit: Fit) -> str:
    reverse = int(fit.direction == 'reverse')
    # The textual marker and the solver's integer flag must agree.
    header = record.header.ljust(45)
    header = header[:39] + ('v' if reverse else ' ') + ' ra26' + header[45:]
    ids = record.ids[:7] + (reverse,)
    eol = '\r\n' if record.raw[0].endswith('\r\n') else '\n'
    values = fit.coefficients + (record.q,)
    return (header+eol + ''.join(f'{v:5d}' for v in ids)+eol
            + ''.join(f'{v:13.5e}' for v in values[:4])+eol
            + ''.join(f'{v:13.5e}' for v in values[4:])+eol)


def convert(items, nuclides, fits, *, modify=True, all_rates=False, mode=0):
    output, log, counts, seen = [], [], collections.defaultdict(collections.Counter), {}
    mass_by_zn = {n.zn: n.mass_excess for n in nuclides}
    mass_by_zn[(0, 0)] = 0.0
    for item in items:
        if not isinstance(item, RateRecord):
            output.append(item); continue
        if not modify or item.iaa or item.section not in CHANNELS:
            output.extend(item.raw); continue
        # Validate only the selectable records. Unrelated legacy entries stay intact.
        if not all_rates and item.label not in THEORY_LABELS:
            output.extend(item.raw); continue
        key = reaction_identity(item, nuclides)
        if not selected(item, key, all_rates, mode):
            output.extend(item.raw); continue
        require(key not in seen,
                f'input:{item.line}: repeated selected reaction, previous line={seen.get(key)}')
        seen[key] = item.line
        counts[item.section]['selected'] += 1
        fit = fits.get(key)
        event = dict(section=item.section, input_line=item.line,
                     reaction=' '.join(format_name(zn) for zn in key), old_label=item.label,
                     status='replaced' if fit else 'kept_no_match',
                     source=fit.source if fit else '', source_line=fit.line if fit else '',
                     direction=fit.direction if fit else '', old_indinv=item.ids[7],
                     new_indinv=int(fit.direction == 'reverse') if fit else item.ids[7],
                     original_q=item.q, smaragd_q=fit.q if fit else '',
                     mass_q='',
                     q_difference=item.q-fit.q if fit else '')
        event['mass_q'] = sum(mass_by_zn[zn] for zn in key[:2])-sum(mass_by_zn[zn] for zn in key[2:])
        if fit is None:
            counts[item.section]['kept_no_match'] += 1
            output.extend(item.raw)
        else:
            counts[item.section][fit.direction] += 1
            counts[item.section]['flag_changed'] += item.ids[7] != int(fit.direction == 'reverse')
            output.append(render_replacement(item, fit))
        log.append(event)
    return ''.join(output), log, {k: dict(v) for k, v in counts.items()}


def tsv(rows):
    out = io.StringIO(newline='')
    if rows:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    return out.getvalue()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=Path('in'))
    parser.add_argument('--smaragd-dir', type=Path, default=Path('../../lib_data/ra26'))
    parser.add_argument('--output-dir', type=Path, default=Path('out'))
    parser.add_argument('--selection', type=int, choices=range(5), default=n_selec,
                        help='0:all, 1:Z<47,N<53, 2:pg/gp, 3:ap/pa/ag/ga, 4:ng/gn/np/pn')
    parser.add_argument('--all-rates', action='store_true', default=Ra26_all,
                        help='Override theory labels and selection, as original Ra26_all did')
    parser.add_argument('--no-modify', action='store_true', default=not (Rate_mod and Ra26))
    parser.add_argument('--force', action='store_true', help='Replace existing generated output files')
    args = parser.parse_args(argv)
    try:
        part_path, rate_path = args.input_dir/'part-v3.z897', args.input_dir/'rate-v3.z897c'
        nuclides = read_partition(part_path)
        items = read_rates(rate_path)
        fits, mismatches = ({}, []) if args.no_modify else read_smaragd(args.smaragd_dir)
        result, events, counts = convert(items, nuclides, fits, modify=not args.no_modify,
                                         all_rates=args.all_rates, mode=args.selection)
        inputs = [part_path, rate_path]
        if not args.no_modify:
            inputs += [args.smaragd_dir/f'fits_{p}.dat' for p in ('p', 'n', 'a')]
        summary = dict(selection=args.selection, all_rates=args.all_rates, modify=not args.no_modify,
                       converter_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                       partition_nuclei=len(nuclides), source_name_mismatches=len(mismatches),
                       selected=len(events), replaced=sum(e['status']=='replaced' for e in events),
                       kept_no_match=sum(e['status']=='kept_no_match' for e in events),
                       by_section=counts, q_policy='Preserve input Q; solver recomputes Q from partition masses',
                       reverse_policy='indinv=1 only for the SMARAGD reverse fit; solver applies PF ratio',
                       inputs=[dict(path=str(p.resolve()), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs])
        outputs = {'rate.dat': result, 'net_nz.dat': ''.join(f'{n.zn[1]:5d}{n.zn[0]:5d}{n.a:5d}\n' for n in nuclides),
                   'replacement_log.tsv': tsv(events), 'source_label_mismatches.tsv': tsv(mismatches),
                   'summary.json': json.dumps(summary, indent=2)+'\n'}
        for name in outputs:
            path = args.output_dir/name
            require(path.resolve() not in {p.resolve() for p in inputs}, f'Output would overwrite input: {path}')
            require(args.force or not path.exists(), f'Output exists: {path}; use another directory or --force')
        args.output_dir.mkdir(parents=True, exist_ok=True)
        for name, content in outputs.items():
            (args.output_dir/name).write_bytes(content.encode('ascii'))
        print(json.dumps(summary, indent=2))
        return 0
    except (RateError, OSError, ValueError, IndexError) as exc:
        print(f'Conversion stopped: {exc}', file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
