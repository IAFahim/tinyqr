import unittest

import tinyqr


class TinyQRTests(unittest.TestCase):
    def test_qr_matrix_dimensions_and_finders(self) -> None:
        m = tinyqr.encode("P1")
        self.assertEqual(len(m), 21)
        self.assertEqual(len(m[0]), 21)
        # Top-left 7x7 finder outer ring
        self.assertEqual(m[1][:7], [1, 0, 0, 0, 0, 0, 1])
        self.assertEqual(m[3][2:5], [1, 1, 1])
        # Dark module at (4*1+9, 8) = (13, 8)
        self.assertEqual(m[13][8], 1)

    def test_qr_url_versions_png_svg_ascii(self) -> None:
        for txt in ("P2", "https://hub.example.ts.net/", "http://192.0.2.10:24808"):
            mat = tinyqr.encode(txt)
            self.assertEqual(len(mat), len(mat[0]))
            png = tinyqr.render_png(mat, scale=8, border=2)
            self.assertTrue(png.startswith(b"\x89PNG\r\n\x1a\n"))
            self.assertGreater(len(png), 150)

            svg = tinyqr.to_svg(txt, scale=8, border=2)
            self.assertIn("<svg", svg)
            self.assertIn("viewBox=", svg)
            self.assertIn("<rect", svg)

            ascii_qr = tinyqr.to_ascii(txt)
            self.assertTrue(any(ch in ascii_qr for ch in ("\u2580", "\u2584", "\u2588")))

    def test_rs_syndrome_zero(self) -> None:
        codewords = tinyqr._encode_codewords(b"P1", 1, "H")
        self.assertEqual(len(codewords), 26)
        self.assertTrue(tinyqr.rs_verify(codewords, 17))


if __name__ == "__main__":
    unittest.main()
