"""The built page is the template, and carries the views and controls the documentation describes.

Needs no data: it reads scripts/explorer_template.html and docs/index.html only.
"""
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = (ROOT / "scripts" / "explorer_template.html").read_text()


class ExplorerPage(unittest.TestCase):
    def test_built_page_is_the_template(self):
        self.assertEqual((ROOT / "docs" / "index.html").read_text(), TEMPLATE, "run scripts/09_build_explorer.py")

    def test_every_tab_has_a_route_and_a_renderer(self):
        tabs = re.findall(r'class="tab[^"]*" data-view="(\w+)"', TEMPLATE)
        self.assertEqual(tabs, ["map", "matrix", "sources", "notes", "documentation"])
        routes = re.search(r"S\.view = \[([^\]]+)\]\.includes\(h\)", TEMPLATE).group(1)
        renderers = re.search(r"\(\{ (map: renderMap[^}]+) \}\)\[S\.view\]", TEMPLATE).group(1)
        for tab in tabs:
            self.assertIn(f"'{tab}'", routes)
            self.assertIn(f"{tab}: ", renderers)

    def test_map_download_formats(self):
        self.assertIn("['png', 'svg', 'pdf'].map", TEMPLATE)
        for name in ("function toSVG", "async function toPDF", "function toPNG", "function exportFigure"):
            self.assertIn(name, TEMPLATE)

    def test_map_footer_reserves_every_variant(self):
        # the legend and caption of every colouring and estimate share one grid cell, so the map keeps its size
        self.assertIn("#mapfoot > div { grid-area: 1 / 1;", TEMPLATE)
        self.assertIn("[['gene', 'projected'], ['comp', est], ['type', est]]", TEMPLATE)


if __name__ == "__main__":
    unittest.main()
