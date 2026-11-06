# Working with AI assistants on Piernik

A practical guide for developers who have not used AI coding assistants (Copilot, Claude Code, Cursor, ...) yet.
Assistants are good at reading unfamiliar Fortran, writing boilerplate, scripts and tests, and explaining legacy code; they are
unreliable about numerical correctness and physics. **You stay responsible for every line you commit.**

## 1. How the repository helps the assistant
The repo contains hierarchical `AGENTS.md` files (see [AGENTS.md](../AGENTS.md) for the scheme). Assistants that understand them
(Copilot, Claude Code, Codex, Cursor, ...) load the one closest to the files you work on, so you do not have to explain the project each time.
If your tool does not read `AGENTS.md` automatically, start a session with: *"Read AGENTS.md and the AGENTS.md of the directory I work in."*
When you change structure (move/rename modules, add a directory, change `setup`) update the matching `AGENTS.md` in the same commit.

## 2. Good first uses (low risk)
- "Explain what module `X` does and who calls it" / "trace what happens between `init_piernik` and the first timestep".
- "Which files do I touch to add a new namelist parameter?" – then check the answer in the code.
- Writing a new `problems/<name>/` by adapting a similar one; `problem.par` variants.
- Post-processing scripts in `visual/` or `python/`; shell helpers; documentation and Doxygen comments.
- Reading CI/compiler error logs and compiler-specific portability problems (ifx vs gfortran, macOS).
- Code review of your own diff *before* the pull request ("look for bugs, MPI deadlocks, missing `only:`, QA violations").

## 3. Riskier uses (verify carefully)
- Edits in `src/grid`, `src/multigrid`, `src/scheme`, `src/IO`: MPI communication, AMR, ghost cells and I/O formats are easy to break in ways that appear only on many ranks or with refinement.
- Anything that changes numerical results: an assistant may "simplify" a formula and silently change the physics. Require a justification and run the gold tests.
- Large refactors across many files in one go – split them into steps you can review (and `make QA` after each).

## 4. A workflow that works
1. **Ask for a plan first**, with the files it intends to touch. Correct it before any code is written.
2. **Small tasks**, one purpose per pull request; commit often so you can `git diff`/revert.
3. **Make the assistant verify**: tell it the commands – `./setup sedov -n`, `./setup <problem>` + short run, `make QA`, `make artifacts -j 3`, `make gold -j 3`. Ask it to run them and report the failures instead of guessing.
4. **Read the diff yourself**, including deleted lines and changes to `problem.par`/`piernik.def`. Check that it did not touch unrelated files or reformat code.
5. Before merging to `master`: `make CI -j`.

## 5. Writing good prompts for Piernik
- Name the problem/module/file: *"in `src/fluids/` ..., test with `problems/sedov`"*.
- State constraints: *"keep results bit-identical"*, *"no new dependencies"*, *"gfortran and ifx must both compile"*, *"keep the module name unique"*.
- Give the failure verbatim (compiler message, `./setup` output, the line from the log) rather than a paraphrase.
- Ask for alternatives and trade-offs on design questions; ask "what could break?" afterwards.
- Do not paste secrets, unpublished data or code from other (non-GPL-compatible) projects into assistant sessions.

## 6. Making the best use of your quota
- Keep sessions focused; start a new one for an unrelated task (long histories cost more and dilute context).
- Point at files and line ranges instead of letting the tool scan `problems/` (≈400 files) or `src/` – use the map in `src/AGENTS.md`.
- Use the cheaper/faster model for routine edits and scripts, the strongest one for design, debugging MPI/AMR issues and reviews.
- Let the tool run the build and tests itself, but avoid `make CI` in a loop: it is long. Use `./setup <problem> -n` or `make QA` for fast feedback.
- Improve the `AGENTS.md` files when you notice the assistant repeatedly misunderstands something – this benefits everybody.

## 7. Known traps in this codebase (assistants fall into these)
- Assuming `src/` subdirectories are namespaces: they are not, module names must be unique globally.
- Message prefix is the *file* name, not module name (`[divB:...]` in `divB.F90`); enforced by `bin/checkmessages.sh`.
- Use the `piernik_MPI_*` wrappers where available. `bin/qa.py` checks `implicit none`, `private`/`public`, `only:`, and related Fortran style issues, and reports magic integers; it does not enforce MPI-wrapper use.
- Python 2 scripts in `python/`, `bin/` (see the lists in AGENTS.md) – do not try to "fix" them unless asked.
- `piernik.def` symbols (compile time) vs `problem.par` parameters (run time) are different mechanisms – don't add `#if` where a namelist parameter would do.
- Changing initial conditions or defaults silently changes gold-test results.

## 8. Etiquette
Mention significant AI assistance in the pull request description (helps reviewers decide how closely to look), keep the licence header and
authorship conventions, and never commit generated artefacts (`obj*/`, `runs/`, images, tmp files).
