# Notes for agents working on ImageD11

The single canonical agent guidance. These rules bias toward caution over speed;
for trivial tasks use judgment.

## 1. Think before coding

- Read the files you will touch, end to end. Fragments give confident wrong answers.
- Say what the current code does, and which commit and branch you read, before
  proposing a change. A defect only counts on the branch the work starts from.
- State your assumptions. If there are several readings of the request, list them;
  do not pick one silently. If something is unclear, stop and ask.
- If a simpler approach exists, say so. Push back when warranted.

## 2. Simplicity first

- The minimum code that solves the problem. Nothing speculative.
- No features, single-use abstractions or configurability that were not asked for.
- No error handling for impossible cases.
- If 200 lines could be 50, rewrite. Would a senior engineer call it overcomplicated?

## 3. Surgical changes

- Every changed line must trace to the request.
- Do not "improve" adjacent code, comments or formatting. Match the existing style.
- Mention unrelated dead code; do not delete it. Remove only what your change orphaned.

## 4. Goal-driven execution

- Turn the task into checks: a bug -> a test that reproduces it, then passes; a
  refactor -> tests pass before and after. For several steps, state a short plan
  with a check for each, and loop until verified.

## Physics: never fudge it

- Never change the physics or crystallography to please a test or match an expectation.
- Never invent figures, numbers or merit that do not come from the data or the
  computation.
- Never add workarounds (tolerated slack, scaled metrics, special-cased inputs) to
  make a result look better. If a result looks bad, find the cause (geometry,
  calibration, indexing, data quality) and fix that.

## Guess or ask

Guess when being wrong costs one edit; ask when it costs a migration.

Yours to decide: names of variables, functions and tests; where a helper lives;
loop structure; equivalent numpy idioms; formatting.

Ask first:

- the name, dtype or meaning of anything written to a file, or of a columnfile column
- the signature or return shape of a function that already has callers
- a default value, where changing it changes results
- the policy for real data: what to drop, merge, mask or pad
- scope: whether a neighbouring problem is in or out

If unsure which list, ask what undoing it would cost. A question built on an
unchecked premise is worse than a guess: check the premise, then ask. If the
answer leaves the design undetermined, ask again.

## Measurements

- If you do not know, say so. "I have not measured that" is a complete answer.
  Never write an inference as a measurement.
- Before claiming a change is worth making, ask what would make it worthless, and
  check that first.
- Every number carries its conditions: warm cache or cold, synthetic or real,
  extrapolated or measured. Do not test a hypothesis on data generated to fit it.
- A saving in one step is often an overhead in the next (compute for I/O, memory
  for passes). Measure both ends or claim neither.
- If you keep a measurement in a design, keep what it implies as well.

## Scans and scale

- omega and dty are not on a regular grid: zig-zag rows, interlaced omega (every
  other frame, gaps filled on a later pass), multi-turn scans with an irrational
  step. See guess_omega_step in sinograms/dataset.py.
- DataSet.shape is a guess (guess_shape warns "irregular scan"). It is not an
  ordering: frames arrive in any order and a column index is not an omega bin.
  ostep and ystep are derived by guessbins, not recorded.
- Omega is circular. Compare with ((a - b + 180) % 360) - 180; never average
  linearly. guessbins applies % 360 only when the span exceeds 360; keep that.
- Peak arrays are huge: pk_props is 17 GB and rc 9.6 GB at 4.24e8 peaks (120 GB and
  68 GB at 3e9). Do not add a per-peak array for something reachable by indirection
  through a per-frame or per-cell table. Say what a new array costs per peak.

## Compatibility

Data from the last 25 years must keep working; some is on tape and will come back.

A contract is anything other people's code or files depend on: names, dtypes and
meanings of datasets, attributes and columnfile columns; what a public function
returns; which datasets a writer writes and a reader requires; defaults that change
results. Adding is safe. Renaming, removing, retyping and redefining are not: new
meaning gets a new name. Readers must tolerate attributes that old files lack.
Before changing a contract, list the callers and say what happens to each. A changed
default is a results change: measure before and after, report both.

Python 2.7 must stay importable and numpy below 2 must work:

- no infix `@` (use numpy.dot), f-strings, `yield from`, keyword-only arguments or
  `nonlocal`; use `from __future__ import print_function` for print()
- no np.random.default_rng (use RandomState); np.float and np.int are gone in numpy 2
- reshape(copy=) needs numpy 2.1; use ImageD11.nputils.reshape_no_copy

Check changed files parse under both Pythons, and test with an old and a current
numpy where behaviour could differ.

Keep `.py` files ASCII (`->`, `--`, `deg`, `Angstrom`, `omega`, `sigma`); then no
`# coding:` line is needed on 2.7. Notebooks may contain unicode.

## Working with others

- Be terse; long text does not get read. Give the conclusion, not what was tried
  and abandoned, in comments, commits and pull requests.
- Stage files by name, one at a time. Never `git add -A` or `git add .`: trees here
  hold large untracked data and test output.
- FABLE-3DXRD is protected. Push to a fork (yours, or someone else's with their
  agreement) and leave the upstream pull request to a person.
- Everyone here is intelligent and constructive, even when a change arrives in a
  shape that does not fit. Say what is worth keeping and what specifically does not
  work. Disagree with evidence: a case where behaviour differs beats an opinion.
  Resolve issues in the design instead of recording them for later.
