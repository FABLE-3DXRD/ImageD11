# Notes for agents working on ImageD11

This file is the single canonical agent guidance. It consolidates the working
notes (merged from agents_v2.md) so there is one place to read before editing.

## Ground rules - never fudge the physics

- Never change the **physics / crystallography** in order to please, make a test
  pass, or match an expectation.
- Never **invent figures, numbers, or merit** that are not derived from the data
  or the actual computation.
- Never **invent "workarounds"** (e.g. adding tolerated slack, scaling metrics,
  or special-casing inputs) to make a result look better. If the geometry or
  matching is wrong, fix the real cause.

If a result "looks bad", investigate the actual cause (geometry, calibration,
indexing, data quality) rather than relaxing criteria or adding ad-hoc
corrections.

## Read before changing

Read the files you will touch, end to end. Fragments and summaries give
confident wrong answers.

State what the current code does before proposing what it should do.

Say which commit you read. A defect only counts if it is on the branch the work
starts from; code on a shelved branch nobody is carrying over is not a defect.
Attribute every claim to the branch and file it lives on.

## Guess or ask

Guess when being wrong costs one edit. Ask when being wrong costs a migration.

Yours to decide: names of variables, functions and tests; where a helper lives;
loop structure; which of two equivalent numpy idioms; formatting.

Ask first, every time:

- the name, dtype or meaning of anything written to a file, or of a columnfile
  column
- the signature or return shape of a function that already has callers
- a default value, where changing it changes results
- which policy to apply to real data: what to drop, merge, mask or pad
- scope: whether a neighbouring problem is in or out

Unsure which list? Ask what undoing it would cost. "Reprocess the data" or "edit
other people's notebooks" puts it on the second.

A question built on an unchecked premise is worse than a guess: it spends
someone's attention and takes them down your wrong path with you. Check the
premise, then ask. A partial answer that leaves the design undetermined is not
an answer - ask again. Two rounds of questions cost less than one wrong
implementation.

## Do not guess, and never write a measurement that is not one

If you do not know, say so. "I have not measured that" is a complete answer and
a better one than a plausible number. Never write an inference in the form of a
measurement.

Before claiming a change is worth making, ask what would have to be true for it
to be worthless, and go and check that first. Implementing before checking
wastes the reviewer's time as well as yours.

Every number carries its conditions. A timing taken with the file in page cache
says nothing about reading over a network or a disk. Say what was warm, what was
synthetic, what was extrapolated, and what was not tested at all.

Do not test a hypothesis on data you generated to fit it. Separate the cost of
the algorithm from the cost of building the test case.

If you keep a measurement in a design, keep what it implies: quoting the
evidence while dropping the conclusion is how a fixed bug comes back.

## A saving in one step is not a saving in the job

A change that removes work from one place usually adds it somewhere else:
trading compute for I/O, or memory for passes. Measure both ends or claim
neither. A saving in one step is often an overhead in the next.

## Do not assume a regular scan

omega and dty are not on a regular grid. Rotation and translation move at
constant speed and the sinogram fills in as the scan runs. In use: rows built up
one ystep per turn; zig-zags with dty constant along a row; interlaced omega,
where a detector with readout takes every other frame and fills the gaps on a
later pass; multi-turn scans whose step may be irrational on purpose so frames
never repeat. See guess_omega_step in sinograms/dataset.py.

DataSet.shape is a guess. guess_shape reshapes omega and dty to
(nrotations, nframes_per_rotation) and warns "irregular scan" when the product
does not match. Do not rely on that shape as an ordering. Frames may arrive in
any order, and a grid column index is not an omega bin.

ostep and ystep are derived by guessbins, not recorded by the beamline.

Omega is circular. Compare with ((a - b + 180) % 360) - 180, and never average
linearly. guessbins applies % 360 only when the span exceeds 360, so a 270-390
scan stays contiguous. Keep that.

## Scale

Peak arrays are large: pk_props is 17 GB and rc at least 9.6 GB at 4.24e8 peaks;
120 GB and 68 GB at 3e9. Do not add a per-peak array for something reachable by
one indirection through a per-frame or per-cell table. Say what a new array
costs per peak before proposing it.

## Backwards compatibility

Data from the last 25 years must keep working. Some is on tape and will come
back.

A contract is anything another person's code or files already depend on:

- names, dtypes and meanings of datasets, attributes and columnfile columns
- the shape and type a public function returns
- the set of datasets a writer writes, and the set a reader requires
- default arguments, where they change results

Adding is safe; renaming, removing, retyping and redefining are not. New meaning
gets a new name, so old and new files can be told apart and migrated. Read paths
must tolerate missing attributes, because old files lack whatever was added
later.

Before changing any of the above, list the callers and say what happens to each.
Changing a default is a results change: measure before and after, report both.

## Python 2 and old numpy

The codebase must stay importable on Python 2.7. When editing, do not use
py3-only syntax:

- No `@` **matrix multiplication** operator - use `numpy.dot(a, b)` instead.
  (`@` decorators on functions/classes are fine; the infix operator is not.)
- No **f-strings** - use `"{}".format(...)` or `%` formatting.
- No bare `print(...)` as a bare call - either avoid printing, or add
  `from __future__ import print_function` to the file.
- No other py3-only features (`yield from`, keyword-only arguments, `nonlocal`,
  type annotations on 2.7 builtins) without guarding by version.

Support numpy below 2 as well as current. ndarray.reshape accepts copy= only
from numpy 2.1; ImageD11.nputils.reshape_no_copy handles the difference. np.float
and np.int were removed in numpy 2.

Check that changed files parse under both Pythons, and test against an old and a
current numpy where behaviour could differ.

## ASCII-only source

Do not introduce non-ASCII characters into `.py` files. They never render
properly across terminals/editors and force a `# coding:` declaration for Python
2.7. Use ASCII equivalents instead: `-`/`--` for dashes, `->` for arrows, `...`
for ellipses, `deg` for degrees, `Angstrom` for Angstrom, `^2`/`1` for
superscripts, and the spelled form for Greek letters (`omega`, `sigma`, `delta`,
`theta`, `epsilon`, `Sigma`). No `# coding: utf-8`-style declarations are needed
because the sources are ASCII. (`.ipynb` files may keep unicode; the browser
usually renders it.)

## Writing

Be terse. Long verbose text does not get read.

Strip working notes. The reader wants the conclusion, not a record of what was
tried and abandoned. Keep "this other approach did not work" out of comments,
commit messages and pull requests.

## Commits

Stage files by name, one at a time. Never `git add -A` or `git add .`. Working
trees here hold large untracked data and test output sitting next to the code,
and sweeping those into a commit is easy to do and tedious to undo.

## Pull requests

FABLE-3DXRD is protected. Push to a fork, your own or someone else's with their
agreement, and leave the upstream pull request to a person.

## People

Everyone here is intelligent and trying to be constructive, including when a
change arrives in a shape that does not fit. Separate the idea from its
implementation: say what is worth keeping, and be specific about what does not
work.

Disagree with evidence. A case where the behaviour differs beats an opinion
about structure. If you cannot produce the case, you may be wrong.

Resolve issues in the design rather than recording them for later.
