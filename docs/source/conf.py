# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# The narrative pages are README.md (quick start) and GUIDE.md (in-depth guide) at the repository root,
# copied here at build time so that there is a single source. Do not edit docs/source/README.md or GUIDE.md.

import os
import re
import sys

root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, os.pardir))
sys.path.insert(0, root_dir)


def _copy_root_docs():
    """Copy README.md and GUIDE.md from the repository root, fixing image paths for docs/source."""
    here = os.path.dirname(os.path.abspath(__file__))
    for name in ('README.md', 'GUIDE.md'):
        with open(os.path.join(root_dir, name), encoding='utf-8') as f:
            text = f.read()
        text = text.replace('](docs/source/_static/', '](_static/')
        with open(os.path.join(here, name), 'w', encoding='utf-8') as f:
            f.write(text)

_copy_root_docs()

# -- Project information -----------------------------------------------------

def _read_version():
    with open(os.path.join(root_dir, 'asari', '__init__.py'), encoding='utf-8') as f:
        return re.search(r'__version__\s*=\s*["\']([^"\']+)["\']', f.read()).group(1)

project = 'asari'
copyright = '2023-2026, Shuzhao Li, Amnah Siddiqa, Maheshwor Thapa, Shujian Zheng'
author = 'Shuzhao Li, Amnah Siddiqa, Maheshwor Thapa, Shujian Zheng'
release = _read_version()
version = release

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'myst_parser',
]

# optional dependencies (extras viz, gc, plot) are not needed to build the API pages
autodoc_mock_imports = ['panel', 'holoviews', 'hvplot', 'bokeh', 'param', 'matplotlib', 'seaborn',
                        'plotly', 'sklearn', 'tkinter']

napoleon_custom_sections = [('Updates', 'params_style'), ('Input', 'params_style'), ('Outputs', 'rubric')]

myst_heading_anchors = 3            # GitHub-style anchors, e.g. GUIDE.md#joining-results

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
