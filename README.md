# tinyqr

Zero-dependency, pure-Python ISO/IEC 18004 QR code encoder (Versions 1-6, Byte mode, ECC `M`/`H`) with PNG, SVG, and terminal half-block ASCII renderers.

## Usage

```python
import tinyqr

matrix = tinyqr.encode("P1")
png_bytes = tinyqr.to_png("P1", scale=8, border=4)
svg_text = tinyqr.to_svg("https://example.ts.net/", scale=8, border=4)
print(tinyqr.to_ascii("https://example.ts.net/"))
```
