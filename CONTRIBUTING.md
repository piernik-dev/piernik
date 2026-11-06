# Contributing to Piernik

Short version – details in [doc/GETTING_STARTED.md](doc/GETTING_STARTED.md).

1. Fork the repository on GitHub and work on a branch; do not commit directly to `master`.
2. Keep one purpose per pull request. Easy starter tasks are listed in [doc/BACKLOG.md](doc/BACKLOG.md) (section A).
3. Before opening the PR run `make QA`; for anything that can change code or results also `make CI -j 3` (see `INSTALL.md`).
4. Follow the conventions enforced by QA (licence header, `[file:routine]` message prefix, `implicit none`, `private`/`only:`, MPI wrappers) –
   see [AGENTS.md](AGENTS.md#conventions) and [doc/HOWTO.md](doc/HOWTO.md) for recipes.
5. If results of standard tests change on purpose, say so in the PR (gold tests will differ).
6. If you changed structure (moved/renamed files, new directory, new setup behaviour), update the matching `AGENTS.md` in the same PR.
7. If an AI assistant wrote a significant part of the change, mention it in the PR; you are responsible for the code. See [doc/AI_DEVELOPMENT.md](doc/AI_DEVELOPMENT.md).

PR checklist:
- [ ] `make QA` passes
- [ ] `./setup <affected problem>` builds and a short run works
- [ ] `make CI -j` (or the relevant group) passes
- [ ] documentation / `AGENTS.md` updated if needed
