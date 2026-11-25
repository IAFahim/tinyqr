"""Pure-Python ISO/IEC 18004 QR Code encoder (Versions 1-6, ECC M/H) + PNG, SVG & ASCII renderers.

Zero external dependencies (uses only stdlib struct + zlib).
"""
from __future__ import annotations

import struct
import zlib
from typing import Dict, List, Optional, Sequence, Tuple

__all__ = [
    "encode",
    "encode_qr",
    "render_png",
    "render_svg",
    "render_ascii",
    "to_png",
    "to_svg",
    "to_ascii",
    "rs_verify",
]

# GF(256) with primitive polynomial x^8 + x^4 + x^3 + x^2 + 1 (0x11d)
_GF_EXP = [0] * 512
_GF_LOG = [0] * 256
_v = 1
for _i in range(255):
    _GF_EXP[_i] = _v
    _GF_LOG[_v] = _i
    _v <<= 1
    if _v & 0x100:
        _v ^= 0x11D
for _i in range(255, 512):
    _GF_EXP[_i] = _GF_EXP[_i - 255]


def _gf_mul(a: int, b: int) -> int:
    if a == 0 or b == 0:
        return 0
    return _GF_EXP[_GF_LOG[a] + _GF_LOG[b]]


def _rs_generator(ec_len: int) -> List[int]:
    poly = [1]
    for i in range(ec_len):
        root = _GF_EXP[i]
        nxt = [0] * (len(poly) + 1)
        for j, coef in enumerate(poly):
            nxt[j] ^= coef
            nxt[j + 1] ^= _gf_mul(coef, root)
        poly = nxt
    return poly


def _rs_encode(data: Sequence[int], ec_len: int) -> List[int]:
    gen = _rs_generator(ec_len)
    rem = list(data) + [0] * ec_len
    for i in range(len(data)):
        factor = rem[i]
        if factor != 0:
            for j in range(1, len(gen)):
                rem[i + j] ^= _gf_mul(gen[j], factor)
    return rem[len(data):]


def rs_verify(block: Sequence[int], ec_len: int) -> bool:
    """Verify a single Reed-Solomon block (data + ec) has zero syndrome."""
    for i in range(ec_len):
        root = _GF_EXP[i]
        val = 0
        for b in block:
            val = _gf_mul(val, root) ^ b
        if val != 0:
            return False
    return True


# Table: (version, ecc) -> (total_data_cw, ec_cw_per_block, group1_blocks, group1_data_cw, group2_blocks, group2_data_cw)
# Supports Versions 1..6 for ECC 'M' and 'H' (Byte mode covers up to 106 bytes)
_EC_SPECS: Dict[Tuple[int, str], Tuple[int, int, int, int, int, int]] = {
    (1, "M"): (16, 10, 1, 16, 0, 0),
    (1, "H"): (9, 17, 1, 9, 0, 0),
    (2, "M"): (28, 16, 1, 28, 0, 0),
    (2, "H"): (16, 28, 1, 16, 0, 0),
    (3, "M"): (44, 26, 1, 44, 0, 0),
    (3, "H"): (26, 22, 2, 13, 0, 0),
    (4, "M"): (64, 18, 2, 32, 0, 0),
    (4, "H"): (36, 16, 4, 9, 0, 0),
    (5, "M"): (86, 24, 2, 43, 0, 0),
    (5, "H"): (46, 22, 2, 11, 2, 12),
    (6, "M"): (108, 16, 4, 27, 0, 0),
    (6, "H"): (60, 28, 4, 15, 0, 0),
}

_ALIGNMENT_CENTERS = {
    1: [],
    2: [6, 18],
    3: [6, 22],
    4: [6, 26],
    5: [6, 30],
    6: [6, 34],
}

_ECC_FORMAT_BITS = {"L": 0b01, "M": 0b00, "Q": 0b11, "H": 0b10}


def _select_version(byte_len: int, ecc: Optional[str] = None) -> Tuple[int, str]:
    if ecc is not None:
        ecc = ecc.upper()
        if ecc not in ("M", "H"):
            raise ValueError(f"Unsupported ECC level: {ecc}")
        for ver in range(1, 7):
            cap = _EC_SPECS[(ver, ecc)][0] - 2  # 4-bit mode + 8-bit length + 4-bit terminator = 16 bits = 2 bytes
            if byte_len <= cap:
                return ver, ecc
        raise ValueError(f"Payload too long ({byte_len} bytes) for QR V1-6 ECC {ecc}")

    # Prefer ECC H for tiny codes (<= 14 bytes -> V1-H or V2-H), then M up to 106 bytes.
    if byte_len <= 7:
        return 1, "H"
    if byte_len <= 14:
        return 2, "H"
    for ver in range(2, 7):
        cap = _EC_SPECS[(ver, "M")][0] - 2
        if byte_len <= cap:
            return ver, "M"
    raise ValueError(f"Payload too long ({byte_len} bytes) for QR V1-6")


def _encode_codewords(payload: bytes, version: int, ecc: str) -> List[int]:
    total_data_cw, ec_cw, g1_blocks, g1_len, g2_blocks, g2_len = _EC_SPECS[(version, ecc)]
    bits: List[int] = []

    def push_bits(val: int, count: int) -> None:
        for shift in range(count - 1, -1, -1):
            bits.append((val >> shift) & 1)

    push_bits(0b0100, 4)  # Byte mode
    push_bits(len(payload), 8)
    for b in payload:
        push_bits(b, 8)

    max_bits = total_data_cw * 8
    term_len = min(4, max_bits - len(bits))
    push_bits(0, term_len)
    while len(bits) % 8 != 0:
        bits.append(0)

    data_cw: List[int] = []
    for i in range(0, len(bits), 8):
        byte = 0
        for b in bits[i:i + 8]:
            byte = (byte << 1) | b
        data_cw.append(byte)

    pad_bytes = (0xEC, 0x11)
    pad_idx = 0
    while len(data_cw) < total_data_cw:
        data_cw.append(pad_bytes[pad_idx % 2])
        pad_idx += 1

    blocks_data: List[List[int]] = []
    blocks_ec: List[List[int]] = []
    pos = 0
    for _ in range(g1_blocks):
        blk = data_cw[pos:pos + g1_len]
        pos += g1_len
        blocks_data.append(blk)
        blocks_ec.append(_rs_encode(blk, ec_cw))
    for _ in range(g2_blocks):
        blk = data_cw[pos:pos + g2_len]
        pos += g2_len
        blocks_data.append(blk)
        blocks_ec.append(_rs_encode(blk, ec_cw))

    out: List[int] = []
    max_data = max(len(b) for b in blocks_data)
    for col in range(max_data):
        for blk in blocks_data:
            if col < len(blk):
                out.append(blk[col])
    for col in range(ec_cw):
        for ec in blocks_ec:
            out.append(ec[col])
    return out


def _draw_finder(grid: List[List[int]], func: List[List[bool]], r0: int, c0: int) -> None:
    n = len(grid)
    for dr in range(-1, 8):
        for dc in range(-1, 8):
            r, c = r0 + dr, c0 + dc
            if not (0 <= r < n and 0 <= c < n):
                continue
            func[r][c] = True
            if dr in (-1, 7) or dc in (-1, 7):
                grid[r][c] = 0
            elif dr in (0, 6) or dc in (0, 6) or (2 <= dr <= 4 and 2 <= dc <= 4):
                grid[r][c] = 1
            else:
                grid[r][c] = 0


def _draw_alignment(grid: List[List[int]], func: List[List[bool]], cr: int, cc: int) -> None:
    for dr in range(-2, 3):
        for dc in range(-2, 3):
            r, c = cr + dr, cc + dc
            func[r][c] = True
            if abs(dr) == 2 or abs(dc) == 2 or (dr == 0 and dc == 0):
                grid[r][c] = 1
            else:
                grid[r][c] = 0


def _format_bits(ecc: str, mask: int) -> int:
    data = (_ECC_FORMAT_BITS[ecc] << 3) | mask
    rem = data << 10
    for shift in range(14, 9, -1):
        if (rem >> shift) & 1:
            rem ^= 0b10100110111 << (shift - 10)
    return ((data << 10) | rem) ^ 0b101010000010010


def _apply_format(grid: List[List[int]], ecc: str, mask: int) -> None:
    n = len(grid)
    fmt = _format_bits(ecc, mask)
    for i in range(15):
        bit = (fmt >> (14 - i)) & 1
        # Copy 1: around top-left finder
        if i <= 5:
            grid[8][i] = bit
        elif i == 6:
            grid[8][7] = bit
        elif i == 7:
            grid[8][8] = bit
        elif i == 8:
            grid[7][8] = bit
        else:
            grid[14 - i][8] = bit
        # Copy 2: bottom-left & top-right
        if i <= 6:
            grid[n - 1 - i][8] = bit
        else:
            grid[8][n - 15 + i] = bit


def _mask_cond(mask: int, r: int, c: int) -> bool:
    if mask == 0:
        return (r + c) % 2 == 0
    if mask == 1:
        return r % 2 == 0
    if mask == 2:
        return c % 3 == 0
    if mask == 3:
        return (r + c) % 3 == 0
    if mask == 4:
        return (r // 2 + c // 3) % 2 == 0
    if mask == 5:
        return (r * c) % 2 + (r * c) % 3 == 0
    if mask == 6:
        return ((r * c) % 2 + (r * c) % 3) % 2 == 0
    return ((r + c) % 2 + (r * c) % 3) % 2 == 0


def _penalty(grid: List[List[int]]) -> int:
    n = len(grid)
    score = 0
    # N1: runs of 5+
    for r in range(n):
        run = 1
        for c in range(1, n):
            if grid[r][c] == grid[r][c - 1]:
                run += 1
            else:
                if run >= 5:
                    score += 3 + (run - 5)
                run = 1
        if run >= 5:
            score += 3 + (run - 5)
    for c in range(n):
        run = 1
        for r in range(1, n):
            if grid[r][c] == grid[r - 1][c]:
                run += 1
            else:
                if run >= 5:
                    score += 3 + (run - 5)
                run = 1
        if run >= 5:
            score += 3 + (run - 5)
    # N2: 2x2 blocks
    for r in range(n - 1):
        for c in range(n - 1):
            v = grid[r][c]
            if v == grid[r][c + 1] == grid[r + 1][c] == grid[r + 1][c + 1]:
                score += 3
    # N3: finder-like pattern 1011101 with 4 zeros on either side
    p1 = [1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0]
    p2 = [0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1]
    for r in range(n):
        row = grid[r]
        for c in range(n - 10):
            sl = row[c:c + 11]
            if sl == p1 or sl == p2:
                score += 40
    for c in range(n):
        col = [grid[r][c] for r in range(n)]
        for r in range(n - 10):
            sl = col[r:r + 11]
            if sl == p1 or sl == p2:
                score += 40
    # N4: dark ratio
    dark = sum(sum(row) for row in grid)
    pct = (dark * 100) // (n * n)
    prev5 = abs((pct // 5) * 5 - 50) // 5
    next5 = abs(((pct // 5) + 1) * 5 - 50) // 5
    score += min(prev5, next5) * 10
    return score


def encode_qr(text: str, ecc: Optional[str] = None) -> List[List[int]]:
    """Return a 2D matrix of 0 (light) and 1 (dark) QR modules."""
    payload = text.encode("utf-8")
    version, chosen_ecc = _select_version(len(payload), ecc)
    codewords = _encode_codewords(payload, version, chosen_ecc)
    bits: List[int] = []
    for cw in codewords:
        for shift in range(7, -1, -1):
            bits.append((cw >> shift) & 1)

    n = 17 + 4 * version
    base = [[0] * n for _ in range(n)]
    func = [[False] * n for _ in range(n)]

    _draw_finder(base, func, 0, 0)
    _draw_finder(base, func, 0, n - 7)
    _draw_finder(base, func, n - 7, 0)

    centers = _ALIGNMENT_CENTERS[version]
    for cr in centers:
        for cc in centers:
            if (cr == 6 and cc == 6) or (cr == 6 and cc == n - 7) or (cr == n - 7 and cc == 6):
                continue
            _draw_alignment(base, func, cr, cc)

    for i in range(8, n - 8):
        func[6][i] = True
        base[6][i] = 1 if i % 2 == 0 else 0
        func[i][6] = True
        base[i][6] = 1 if i % 2 == 0 else 0

    func[n - 8][8] = True
    base[n - 8][8] = 1

    # Reserve format info cells
    for i in range(9):
        if i != 6:
            func[8][i] = True
            func[i][8] = True
    for i in range(1, 9):
        func[8][n - i] = True
    for i in range(1, 8):
        func[n - i][8] = True

    # Place data bits in zig-zag
    bit_idx = 0
    upward = True
    col = n - 1
    while col > 0:
        if col == 6:
            col -= 1
        rows = range(n - 1, -1, -1) if upward else range(n)
        for r in rows:
            for c in (col, col - 1):
                if not func[r][c]:
                    base[r][c] = bits[bit_idx] if bit_idx < len(bits) else 0
                    bit_idx += 1
        upward = not upward
        col -= 2

    best_grid: Optional[List[List[int]]] = None
    best_score: Optional[int] = None
    for mask in range(8):
        cand = [row[:] for row in base]
        for r in range(n):
            for c in range(n):
                if not func[r][c] and _mask_cond(mask, r, c):
                    cand[r][c] ^= 1
        _apply_format(cand, chosen_ecc, mask)
        sc = _penalty(cand)
        if best_score is None or sc < best_score:
            best_score = sc
            best_grid = cand

    assert best_grid is not None
    return best_grid


encode = encode_qr


def render_png(matrix: Sequence[Sequence[int]], scale: int = 8, border: int = 4) -> bytes:
    """Render QR matrix as an 8-bit grayscale PNG byte string."""
    n = len(matrix)
    size = (n + border * 2) * scale
    raw = bytearray()
    white_row = b"\x00" + (b"\xff" * size)
    for _ in range(border * scale):
        raw.extend(white_row)
    for row in matrix:
        line = bytearray(b"\x00")
        line.extend(b"\xff" * (border * scale))
        for cell in row:
            line.extend((b"\x00" if cell else b"\xff") * scale)
        line.extend(b"\xff" * (border * scale))
        row_bytes = bytes(line)
        for _ in range(scale):
            raw.extend(row_bytes)
    for _ in range(border * scale):
        raw.extend(white_row)

    def chunk(tag: bytes, data: bytes) -> bytes:
        crc = zlib.crc32(tag + data) & 0xFFFFFFFF
        return struct.pack(">I", len(data)) + tag + data + struct.pack(">I", crc)

    ihdr = struct.pack(">IIBBBBB", size, size, 8, 0, 0, 0, 0)
    idat = zlib.compress(bytes(raw), 6)
    return b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", idat) + chunk(b"IEND", b"")


def render_svg(matrix: Sequence[Sequence[int]], scale: int = 8, border: int = 4) -> str:
    """Render QR matrix as a crisp standalone SVG document."""
    n = len(matrix)
    total = n + border * 2
    px = total * scale
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {total} {total}" '
        f'width="{px}" height="{px}" shape-rendering="crispEdges">',
        f'<rect width="{total}" height="{total}" fill="#fff"/>',
        '<g fill="#000">',
    ]
    for r, row in enumerate(matrix):
        y = r + border
        c = 0
        while c < n:
            if row[c]:
                start = c
                while c < n and row[c]:
                    c += 1
                parts.append(f'<rect x="{start + border}" y="{y}" width="{c - start}" height="1"/>')
            else:
                c += 1
    parts.append("</g></svg>")
    return "".join(parts)


def render_ascii(matrix: Sequence[Sequence[int]], border: int = 2) -> str:
    """Render a compact half-block Unicode QR code for terminals."""
    n = len(matrix)
    padded = [[0] * (n + 2 * border) for _ in range(border)]
    for r in matrix:
        padded.append([0] * border + list(r) + [0] * border)
    padded.extend([[0] * (n + 2 * border) for _ in range(border)])
    lines: List[str] = []
    for y in range(0, len(padded), 2):
        top = padded[y]
        bot = padded[y + 1] if y + 1 < len(padded) else [0] * len(top)
        chars: List[str] = []
        for t, b in zip(top, bot):
            # Inverted colors look sharp on dark terminals (white border, black modules)
            if t == 0 and b == 0:
                chars.append("\u2588")
            elif t == 0 and b == 1:
                chars.append("\u2580")
            elif t == 1 and b == 0:
                chars.append("\u2584")
            else:
                chars.append(" ")
        lines.append("".join(chars))
    return "\n".join(lines)


def to_png(text: str, scale: int = 8, border: int = 4, ecc: Optional[str] = None) -> bytes:
    """Encode text and render directly to PNG bytes."""
    return render_png(encode_qr(text, ecc=ecc), scale=scale, border=border)


def to_svg(text: str, scale: int = 8, border: int = 4, ecc: Optional[str] = None) -> str:
    """Encode text and render directly to an SVG string."""
    return render_svg(encode_qr(text, ecc=ecc), scale=scale, border=border)


def to_ascii(text: str, border: int = 2, ecc: Optional[str] = None) -> str:
    """Encode text and render directly to compact terminal Unicode blocks."""
    return render_ascii(encode_qr(text, ecc=ecc), border=border)
