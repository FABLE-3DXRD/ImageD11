# Sphinx configuration for the ImageD11 documentation.
#
# Build with:  pip install -r requirements.txt ; make html
# Nothing is compiled and no notebook is executed.

import glob
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
GITHUB = "https://github.com/FABLE-3DXRD/ImageD11/blob/master/"

# Use the c2ImageD11 wheel instead of compiling ImageD11._cImageD11
os.environ.setdefault("IMAGED11_USE_C2", "1")
# autodoc imports every module: never open plot windows
os.environ["MPLBACKEND"] = "Agg"
sys.path.insert(0, ROOT)
import ImageD11

project = "ImageD11"
copyright = "Jon Wright and ImageD11 contributors"
version = release = ImageD11.__version__

extensions = [
    "sphinx.ext.apidoc",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "myst_nb",
]

root_doc = "index"
exclude_patterns = ["_build", "**.ipynb_checkpoints"]
html_theme = "furo"
html_title = "ImageD11 " + release

# Notebooks are shown as saved in git (outputs, if any, are not refreshed)
nb_execution_mode = "off"

# API pages, generated at build time into api/ (gitignored)
apidoc_modules = [
    {
        "path": os.path.join(ROOT, "ImageD11"),
        "destination": "api",
        # absolute: sphinx takes relative patterns from this folder
        "exclude_patterns": [
            os.path.join(ROOT, "ImageD11", "sandbox"),
            os.path.join(ROOT, "ImageD11", "depreciated"),
        ],
        "separate_modules": True,
        "module_first": True,
    },
]
# Optional dependencies of a few modules (GUI, forward model)
autodoc_mock_imports = ["ipywidgets", "joblib", "pyFAI", "silx", "sklearn"]


def copy_notebook(src, dst):
    """Copy a notebook into the docs tree, adding a title if it has none
    and a note saying where the source lives."""
    with open(src) as f:
        nb = json.load(f)
    cells = nb["cells"]
    for cell in cells:  # myst-nb ignores images pasted into markdown cells
        for name, data in cell.pop("attachments", {}).items():
            mime, b64 = list(data.items())[0]
            uri = "data:%s;base64,%s" % (mime, "".join(b64))
            cell["source"] = "".join(cell["source"]).replace("attachment:" + name, uri)
    title_at = None
    for i, cell in enumerate(cells):
        if cell["cell_type"] == "markdown":
            lines = "".join(cell["source"]).splitlines()
            if any(line.startswith("# ") for line in lines):
                title_at = i
                break
    if title_at is None:
        name = os.path.splitext(os.path.basename(src))[0]
        cells.insert(0, {"cell_type": "markdown", "metadata": {}, "source": "# " + name})
        title_at = 0
    path = os.path.relpath(src, ROOT).replace(os.sep, "/")
    note = (
        "Source: [%s](%s%s). Shown as saved in git; "
        "the documentation build does not run it." % (path, GITHUB, path)
    )
    cells.insert(title_at + 1, {"cell_type": "markdown", "metadata": {}, "source": note})
    if not os.path.isdir(os.path.dirname(dst)):
        os.makedirs(os.path.dirname(dst))
    with open(dst, "w") as f:
        json.dump(nb, f)


# Copy the notebooks into the source tree (nbgui/ and notes/ are gitignored)
for folder in ("nbgui", "notes"):
    shutil.rmtree(os.path.join(HERE, folder), ignore_errors=True)
NBGUI = os.path.join(ROOT, "ImageD11", "nbGui")
for src in glob.glob(os.path.join(NBGUI, "**", "*.ipynb"), recursive=True):
    if ".ipynb_checkpoints" not in src:
        copy_notebook(src, os.path.join(HERE, "nbgui", os.path.relpath(src, NBGUI)))
for src in glob.glob(os.path.join(ROOT, "docs", "*.ipynb")):
    copy_notebook(src, os.path.join(HERE, "notes", os.path.basename(src)))

HISTORICAL_NOTE = """
.. note::

   This page is from the older ImageD11 documentation and has not been
   reviewed recently. Much of it still applies. The current processing
   workflows are the notebooks under :doc:`/workflows`.

"""


def add_historical_note(app, docname, source):
    if docname.startswith("historical/") and docname != "historical/index":
        # insert after the page title, which is underlined, maybe overlined
        lines = source[0].split("\n")
        ruled = [i for i, line in enumerate(lines) if line.strip() and set(line.strip()) <= set("=-~*#^")]
        end = ruled[0] + 1
        if ruled[0] == 0 or not lines[ruled[0] - 1].strip():  # overlined
            end = ruled[0] + 3
        source[0] = "\n".join(lines[:end] + [HISTORICAL_NOTE] + lines[end:])


# c2ImageD11 docstrings have Checks and Overloads sections holding C code
C2_SECTIONS = ["Checks", "Overloads"]
napoleon_custom_sections = C2_SECTIONS


def c2_literal_sections(app, what, name, obj, options, lines):
    """Make the C2_SECTIONS contents literal text (runs before napoleon)"""
    if not name.startswith("c2ImageD11"):
        return
    out = []
    literal = False
    for i, line in enumerate(lines):
        underline = i + 1 < len(lines) and lines[i + 1].strip() and set(lines[i + 1].strip()) == {"-"}
        if underline:
            literal = line.strip() in C2_SECTIONS
            if literal:
                out += [line, lines[i + 1], ".. code-block:: none", ""]
                continue
        if literal and set(line.strip()) == {"-"}:
            continue
        out.append("    " + line if literal else line)
    lines[:] = out


def setup(app):
    app.connect("source-read", add_historical_note)
    app.connect("autodoc-process-docstring", c2_literal_sections, priority=400)
