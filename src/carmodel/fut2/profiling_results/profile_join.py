"""Join futhark profile .summary (cost / time) with .log (provenance).

Usage: python3 profile_join.py <prof_subdir> [--by-source] [--top N] [--out FILE]
  prof_subdir : directory containing the .log/.summary/.timeline triplet
                (typically `<name>.prof/<entry>/`)
  --by-source : aggregate by source provenance instead of kernel name
  --top N     : print only top N rows by total time (default 30)
  --out FILE  : write output to FILE instead of stdout
"""
import re
import sys
from collections import defaultdict
from contextlib import nullcontext
from pathlib import Path


def parse_log(path):
    """Return {kernel_name: provenance}."""
    mapping = {}
    pending = None
    with open(path) as f:
        for line in f:
            if line.startswith('Event: '):
                pending = line[len('Event: '):].rstrip()
            elif pending and line.lstrip().startswith('at: '):
                mapping.setdefault(pending, line.strip()[len('at: '):])
                pending = None
            else:
                pending = None
    return mapping


_TIME_RE = re.compile(r'([\d.]+)μs')
_ROW_RE = re.compile(
    r'^(\S+)\s+(\d+)\s+([\d.]+)μs\s+([\d.]+)μs\s+([\d.]+)μs\s+([\d.]+)μs\s+([\d.]+)\s*$'
)

def parse_summary(path):
    """Yield dicts {name, count, sum_us, avg_us, min_us, max_us, fraction}."""
    with open(path) as f:
        for line in f:
            m = _ROW_RE.match(line)
            if m:
                yield {
                    'name':     m.group(1),
                    'count':    int(m.group(2)),
                    'sum_us':   float(m.group(3)),
                    'avg_us':   float(m.group(4)),
                    'min_us':   float(m.group(5)),
                    'max_us':   float(m.group(6)),
                    'fraction': float(m.group(7)),
                }


def find_files(prof_dir):
    prof_dir = Path(prof_dir)
    logs = list(prof_dir.glob('*.log'))
    summaries = list(prof_dir.glob('*.summary'))
    if not logs or not summaries:
        sys.exit(f"no .log/.summary in {prof_dir}")
    return logs[0], summaries[0]


def main(argv):
    if len(argv) < 2:
        sys.exit(__doc__)
    prof_dir = argv[1]
    by_source = '--by-source' in argv
    top = 30
    if '--top' in argv:
        top = int(argv[argv.index('--top') + 1])
    out = None
    if '--out' in argv:
        out = argv[argv.index('--out') + 1]

    log_path, summary_path = find_files(prof_dir)
    prov = parse_log(log_path)
    rows = list(parse_summary(summary_path))

    with (open(out, 'w') if out else nullcontext(sys.stdout)) as fh:
        if by_source:
            agg = defaultdict(lambda: {'count': 0, 'sum_us': 0.0,
                                       'frac': 0.0, 'kernels': set()})
            for r in rows:
                key = prov.get(r['name'], '<unknown>')
                # Keep untagged kernels distinguishable by suffixing the kernel
                # name; otherwise every "unknown" kernel collapses into one bucket.
                if key in ('<unknown>', 'unknown'):
                    key = f'<unknown:{r["name"]}>'
                a = agg[key]
                a['count']   += r['count']
                a['sum_us']  += r['sum_us']
                a['frac']    += r['fraction']
                a['kernels'].add(r['name'])
            items = sorted(agg.items(), key=lambda kv: -kv[1]['sum_us'])[:top]
            print(f"{'source':<45} {'#kernels':>8} {'count':>8} {'sum (ms)':>12} {'frac':>7}",
                  file=fh)
            print('-' * 85, file=fh)
            for src, a in items:
                print(f"{src:<45} {len(a['kernels']):>8} {a['count']:>8} "
                      f"{a['sum_us']/1000:>12.2f} {a['frac']:>7.4f}",
                      file=fh)
        else:
            rows.sort(key=lambda r: -r['sum_us'])
            print(f"{'kernel':<45} {'count':>6} {'sum (ms)':>11} {'avg (μs)':>10} "
                  f"{'frac':>7}  source",
                  file=fh)
            print('-' * 110, file=fh)
            for r in rows[:top]:
                print(f"{r['name']:<45} {r['count']:>6} {r['sum_us']/1000:>11.2f} "
                      f"{r['avg_us']:>10.2f} {r['fraction']:>7.4f}  "
                      f"{prov.get(r['name'], '<unknown>')}",
                      file=fh)


if __name__ == '__main__':
    main(sys.argv)
