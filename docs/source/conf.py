# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import re
import sys
from datetime import date
from pathlib import Path

SETRIQ_HOME = Path(__file__).parent.parent.parent.resolve()
sys.path.insert(0, SETRIQ_HOME.as_posix())


# -- Project information -----------------------------------------------------

project = "setriq"
copyright = "2021-%s, Benjamin Tenmann" % date.today().year
author = "Benjamin Tenmann"


# The full version, including alpha/beta/rc tags
def get_version() -> str:
    return (
        re.findall(r"\[([\d.]+)]", (SETRIQ_HOME / "CHANGELOG.md").read_text())
        or ["0.1.0"]
    )[0]


release = get_version()


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.napoleon",  # enable numpy doc
    "sphinx.ext.autodoc",  # auto doc generation
    "sphinx.ext.autosummary",  # auto generate summaries
    "sphinx.ext.mathjax",  # enable LaTeX docs
    "sphinx.ext.autosectionlabel",  # allow auto section generation
]
autosectionlabel_prefix_document = True
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "pydata_sphinx_theme"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]


# -- Auto API -----------------------------------------------------------------
extensions.append("autoapi.extension")
autoapi_dirs = [(SETRIQ_HOME / "src").as_posix()]

autoapi_options = [
    "members",
    "undoc-members",
    "show-inheritance",
    "show-module-summary",
    "imported-members",
]


# -- MyST ---------------------------------------------------------------------
extensions.append("myst_parser")


# -- External links -----------------------------------------------------------
extensions.append("sphinx.ext.extlinks")
extlinks = {
    "doi": ("https://dx.doi.org/%s", "doi:"),
}

html_theme_options = {
    "github_url": "https://github.com/BenTenmann/setriq",
}
