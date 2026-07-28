## 1. Archive existing docs

- [x] 1.1 Rename `docs/` → `_docs/` (git mv, preserve history)
- [x] 1.2 Update any references to `docs/` in `.gitlab-ci.yml`, `_quarto.yml`, `_quarto-public.yml`, `_quarto-private.yml`, `README.md`, `index.qmd` to point at `_docs/` (or remove if the Quarto build is being retired)
- [x] 1.3 Confirm `_freeze/docs` (Quarto cache) either moves alongside or is safely regenerable; relocate if needed

## 2. Run BestieTemplate scaffold

- [x] 2.1 Install/confirm access to `BestieTemplate.jl` (Copier-based generator)
- [x] 2.2 Run `BestieTemplate.new_pkg_quick` against the repo root with `package_owner = "thealanjason"`, `authors = "Tanish Jain, Alan Correa"`, `JuliaMinVersion = "1.12"`, copyright year 2026, strategy level 4 ("Robust"), non-interactive (`StrategyConfirmIncluded`/`StrategyReviewExcluded` = false)
- [x] 2.3 Review the full diff of generated/modified files before staging anything
- [x] 2.4 Reconcile `Project.toml` (merge generated fields with existing `[deps]`/`[compat]`/`version`, add both authors, bump `julia` compat from `"1.10"` to `"1.12"`)
- [x] 2.4a Patch generated `Project.toml`'s `uuid` back to the existing `d94f865e-9be7-4415-9b4a-0e34ca8c7922` (do not keep a freshly generated UUID)
- [x] 2.5 Reconcile `.gitignore`, `README.md`, `LICENSE` (keep existing content where still accurate; adopt generated badges/structure where useful)
- [x] 2.6 Verify `.copier-answers.yml` records the correct owner, authors, and strategy level

## 3. Documentation skeleton

- [x] 3.1 Create `docs/src/index.md` landing page
- [x] 3.2 Create `docs/src/10-tutorials/` and add `getting-started.md`
- [x] 3.3 Create `docs/src/20-how-to/` directory
- [x] 3.4 Create `docs/src/30-explanation/` directory
- [x] 3.5 Confirm `docs/src/95-reference.md` (or `90-reference/`) `@autodocs`/`@index`/`@contents` blocks are present per BestieTemplate's generated `make.jl`
- [x] 3.6 Update `docs/make.jl` `titles` dict and `pages` list to reflect the new section folders

## 4. Migrate tutorial content

- [x] 4.1 Port `_docs/report/tutorial/tutorial.qmd` (270 lines) into `docs/src/10-tutorials/getting-started.md`, converting Quarto-specific syntax (callouts, code-block execution directives, cross-refs) to plain Documenter Markdown

## 5. Migrate explanation content

- [x] 5.1 Port `_docs/report/00_intro.qmd` into `docs/src/30-explanation/background.md`
- [x] 5.2 Fold the single FAQ entry ("Why Julia?") from `_docs/report/faq.qmd` into `background.md` as a short aside
- [x] 5.3 Port `_docs/report/implementation/implementation.qmd` (397 lines) into `docs/src/30-explanation/numerics.md`
- [x] 5.4 Port `_docs/report/20_interpolators.qmd` into `docs/src/30-explanation/interpolators.md`
- [x] 5.5 Port `_docs/report/30_appendix.qmd` (Jacobian/DCM derivation) into `docs/src/30-explanation/jacobians.md`
- [x] 5.6 Confirm `_docs/report/20_conclusion.qmd` is intentionally not migrated (no body content to lose)
- [x] 5.7 Convert all math (LaTeX) and any Quarto-specific formatting in the above to standard Documenter-compatible Markdown/LaTeX

## 6. Map and write how-to guides

- [x] 6.1 Enumerate current public/exported functions across `src/module/{parser,mesh,geometry,init,initialConditions,interpolators,precomputations,quality,reordering,solver,cache,utils,visualization}`
- [x] 6.2 Cross-reference against `_docs/report/API/API.qmd`'s "for external use" section to identify which capabilities need a how-to page
- [x] 6.3 Write one short how-to page per notable capability (e.g., parsing a mesh, running a solve, exporting a solution) under `docs/src/20-how-to/`

## 7. Docstrings and doctests

- [x] 7.1 For each function referenced in `_docs/report/API/API.qmd`, write a fresh docstring on its current (post-rename) definition, repurposing the old prose description
- [x] 7.2 Add docstrings to any additional small/pure public functions not previously covered by API.qmd
- [x] 7.3 Add `jldoctest` blocks with runnable examples to small/pure functions where feasible
- [x] 7.4 Run `Documenter.doctest` locally against the package and fix failures
- [x] 7.5 Confirm `@autodocs` in the reference page picks up all intended docstrings and none reference stale `FAS.*` naming

## 8. CI wiring

- [x] 8.1 Run `JuliaFormatter` and fix formatting violations across `src/` and `test/`
- [x] 8.2 Run markdownlint/yamllint/lychee locally and fix violations
- [ ] 8.3 Enable and verify `Test.yml` passes on a draft PR
- [x] 8.4 Configure GitHub Pages source and `DOCUMENTER_KEY`/deploy secret for `Docs.yml`
- [ ] 8.5 Verify `Docs.yml` builds and deploys successfully on a draft PR / push to default branch
- [ ] 8.6 Verify `Lint.yml`, `CompatHelper`, `TagBot`, and the Copier update-check workflow all run without error
- [x] 8.7 Set up Codecov integration (`codecov.yml`) and confirm coverage upload works

## 9. Final review

- [x] 9.1 Confirm every substantive point from the original `_docs/report/*` content has a home in the new `docs/src/` structure (cross-check against design.md's migration mapping)
- [x] 9.2 Confirm `_docs/` is retained in full (report/, presentation/, orga/, wiki/, images/, references.bib)
- [x] 9.3 Confirm no stale `FAS.*` naming remains anywhere in the new docs
- [x] 9.4 Grep `docs/src/` for any reference to `_docs` and remove — new docs must not link to the archived material
- [ ] 9.5 Merge and confirm the live GitHub Pages site renders correctly end-to-end
