#!/usr/bin/env python3
"""Render the offline RST reading guide as a PDF using only Matplotlib.

This intentionally handles the guide's small RST subset: numbered sections,
paragraphs, bullets, inline literals, and indented literal blocks. It is not a
general RST converter. No native MEDICC modules, network, or LaTeX are needed.
"""
import argparse
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re
import textwrap

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
WIDTH, HEIGHT = 864, 648
INK, MUTED, ACCENT = '#183247', '#526878', '#087f8c'


def sections(source):
    headings = list(re.finditer(r'^(\d{2}) / (.+)\n-+\n', source, re.M))
    if not headings:
        raise ValueError('No numbered guide sections found')
    for i, match in enumerate(headings):
        if int(match[1]) != i + 1:
            raise ValueError('Guide sections must be consecutively numbered')
        end = headings[i + 1].start() if i + 1 < len(headings) else len(source)
        yield match[1], match[2], source[match.end():end]


def blocks(body):
    lines = body.strip('\n').splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        if line.startswith('    '):
            group = []
            while i < len(lines) and (lines[i].startswith('    ') or not lines[i].strip()):
                group.append(lines[i])
                i += 1
            yield 'code', textwrap.dedent('\n'.join(group).rstrip('\n'))
            continue
        kind = 'bullet' if line.startswith('* ') else 'paragraph'
        group = [line[2:] if kind == 'bullet' else line]
        i += 1
        while i < len(lines) and lines[i].strip() and not lines[i].startswith(('    ', '* ')):
            group.append(lines[i].strip())
            i += 1
        content = ' '.join(group).replace('``', '')
        if content.endswith('::'):
            content = content[:-1]
        yield kind, content


def text_width(renderer, text, size, mono=False, bold=False):
    prop = FontProperties(family='DejaVu Sans Mono' if mono else 'DejaVu Sans',
                          size=size, weight='bold' if bold else 'normal')
    return renderer.get_text_width_height_descent(text, prop, False)[0]


def wrap(renderer, text, size, width, bold=False):
    result, line = [], ''
    for word in text.split():
        candidate = f'{line} {word}'.strip()
        if text_width(renderer, candidate, size, bold=bold) <= width:
            line = candidate
        else:
            if line:
                result.append(line)
            if text_width(renderer, word, size, bold=bold) > width:
                raise ValueError(f'Unbreakable text too wide: {word}')
            line = word
    if line:
        result.append(line)
    return result


def layout(renderer, content, size):
    result = []
    for kind, text in content:
        if kind == 'code':
            font = size - 1.1
            lines = text.splitlines()
            available = WIDTH - 132
            while max(text_width(renderer, line, font, mono=True) for line in lines) > available:
                font -= 0.25
                if font < 9:
                    raise ValueError('Literal block too wide at readable font size')
            leading, padding = font * 1.34, 22
        else:
            font = size - 2 if text.startswith('Code:') else size
            leading, padding = font * 1.43, 10
            if kind == 'bullet':
                padding = 4
            lines = wrap(renderer, text, font, WIDTH - 96 - (18 if kind == 'bullet' else 0))
        result.append((kind, lines, font, leading, len(lines) * leading + padding))
    return result


def render(number, title, body, total, revision, edition):
    fig = plt.figure(figsize=(WIDTH / 72, HEIGHT / 72), dpi=72)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, WIDTH), ylim=(0, HEIGHT))
    ax.axis('off')
    renderer = fig.canvas.get_renderer()
    ax.add_patch(Rectangle((0, HEIGHT - 8), WIDTH, 8, color=ACCENT))
    ax.text(48, HEIGHT - 39, f'MEDICC  /  ALGORITHM READING GUIDE  /  {number}',
            fontsize=10, color=ACCENT, weight='bold', va='top')
    title_lines = wrap(renderer, title, 25, WIDTH - 96, bold=True)
    y = HEIGHT - 69
    for line in title_lines:
        ax.text(48, y, line, fontsize=25, color=INK, weight='bold', va='top')
        y -= 32
    y -= 19
    content = list(blocks(body))
    size = 14.5
    while True:
        plan = layout(renderer, content, size)
        if sum(row[-1] for row in plan) <= y - 57:
            break
        size -= 0.25
        if size < 11:
            raise ValueError(f'Page {number} overflows; shorten its source')
    for kind, lines, font, leading, height in plan:
        x = 48
        if kind == 'code':
            ax.add_patch(Rectangle((48, y - height + 7), WIDTH - 96, height - 7,
                                   facecolor='#eff5f8', edgecolor='none'))
            ax.add_patch(Rectangle((48, y - height + 7), 3, height - 7,
                                   facecolor=ACCENT, edgecolor='none'))
            x = 64
            text_y = y - 10
        else:
            text_y = y
        if kind == 'bullet':
            ax.text(49, text_y, '\u2022', fontsize=font, color=ACCENT, va='top')
            x += 18
        for line in lines:
            ax.text(x, text_y, line, fontsize=font, color=INK,
                    family='DejaVu Sans Mono' if kind == 'code' else 'DejaVu Sans',
                    va='top')
            text_y -= leading
        y -= height
    ax.plot([48, WIDTH - 48], [40, 40], color='#d7e2e8', linewidth=0.8)
    ax.text(48, 25, f'CAPPED WGD  /  reviewed source {revision}  /  {edition}',
            fontsize=8.5, color=MUTED, va='center')
    ax.text(WIDTH - 48, 25, f'{int(number):02d} / {total:02d}', fontsize=9,
            color=ACCENT, ha='right', va='center')
    # Catch actual font-layout issues rather than silently clipping PDF content.
    fig.canvas.draw()
    for item in ax.texts:
        bbox = item.get_window_extent(fig.canvas.get_renderer())
        if bbox.x0 < 0 or bbox.x1 > WIDTH or bbox.y0 < 0 or bbox.y1 > HEIGHT:
            raise ValueError(f'Page {number} has text outside the page: {item.get_text()}')
    return fig, size


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=ROOT / 'doc/algorithm_reading_guide.rst')
    parser.add_argument('--output', type=Path, default=ROOT / 'doc/full_fst_algorithms.pdf')
    parser.add_argument('--preview-dir', type=Path, help='Optionally render PNG page previews')
    args = parser.parse_args()
    source = args.source.read_text(encoding='utf-8')
    revision = re.search(r'implementation commit ``([0-9a-f]+)``', source).group(1)
    pages = list(sections(source))
    if args.preview_dir:
        args.preview_dir.mkdir(parents=True, exist_ok=True)
    matplotlib.rcParams.update({'pdf.fonttype': 42, 'font.family': 'DejaVu Sans',
                                'text.usetex': False})
    edition = re.search(r'as reviewed on\s+(\d{1,2} \w+ \d{4})', source).group(1)
    date = datetime.strptime(edition, '%d %B %Y').replace(tzinfo=timezone.utc)
    metadata = dict(Title='Understanding the capped-WGD algorithm', Author='MEDICC project / Codex',
                    Subject=f'Offline algorithm reading guide; reviewed implementation {revision}',
                    Keywords='source-sha256:' + hashlib.sha256(source.encode('utf-8')).hexdigest(),
                    CreationDate=date, ModDate=date)
    # Replace the output only after every page has rendered successfully.
    temporary = args.output.with_suffix('.tmp.pdf')
    try:
        with PdfPages(temporary, metadata=metadata) as pdf:
            for number, title, body in pages:
                fig, size = render(number, title, body, len(pages), revision, edition)
                try:
                    pdf.savefig(fig)
                    if args.preview_dir:
                        fig.savefig(args.preview_dir / f'page-{number}.png', dpi=110)
                finally:
                    plt.close(fig)
                print(f'Page {number}: body {size:g} pt; {title}')
        temporary.replace(args.output)
    finally:
        if temporary.exists():
            temporary.unlink()
    print(f'Wrote {len(pages)} pages to {args.output}')


if __name__ == '__main__':
    main()
