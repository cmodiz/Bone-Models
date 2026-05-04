import os
import sys
# Reliably point to the root regardless of where Sphinx executes from
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Bone Models'
copyright = '2025, Corinna Modiz'
author = 'Corinna Modiz'
release = '1.0.0'

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# extensions = ['sphinx.ext.autodoc', 'sphinx.ext.napoleon', 'sphinx.ext.viewcode']
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', "sphinx.ext.viewcode", 'sphinx.ext.napoleon', 'sphinx.ext.mathjax', "sphinx.ext.intersphinx", "sphinx.ext.extlinks"]
extlinks = {
    'doi': ('https://doi.org/%s', 'DOI: %s')
}
templates_path = ['_templates']
exclude_patterns = []

autoclass_content = "both"

autodoc_mock_imports = ["numpy", "scipy", "pandas", "matplotlib", "fipy"]
