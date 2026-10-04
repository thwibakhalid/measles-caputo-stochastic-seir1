# How to publish this repository and obtain the DOI

The repository is complete and tested. What remains requires your GitHub
account, so it has to be done from your own machine.

---

## Step 1 — create the repository on GitHub

1. Go to <https://github.com/new>.
2. Repository name: `measles-caputo-stochastic-seir`
3. Description: *Code and data for a Caputo fractional-order stochastic SEIR
   model of measles with piecewise constant vaccination.*
4. Set it to **Public** — a private repository cannot be archived by Zenodo
   and the editor requires a working link.
5. Do **not** tick "Add a README", "Add .gitignore" or "Choose a license":
   all three are already in the bundle and GitHub would create a conflict.
6. Click **Create repository**.

---

## Step 2 — push the files

Unzip the bundle, open a terminal inside the unzipped folder, and run the
following. Replace `YOUR-USERNAME` with your GitHub username.

```bash
git init
git add .
git commit -m "Code and data for the Caputo fractional stochastic SEIR measles model"
git branch -M main
git remote add origin https://github.com/YOUR-USERNAME/measles-caputo-stochastic-seir.git
git push -u origin main
```

If `git push` asks for a password, GitHub no longer accepts account
passwords. Create a token at **Settings → Developer settings → Personal
access tokens → Tokens (classic)**, give it the `repo` scope, and paste the
token when prompted for the password.

---

## Step 3 — connect Zenodo and mint the DOI

A GitHub link alone is not a citable archive: anyone can force-push over it,
so it is not a stable record. Zenodo gives you a permanent DOI, which is what
the journal wants in the Data Availability section.

1. Go to <https://zenodo.org> and sign in **with your GitHub account**.
2. Open <https://zenodo.org/account/settings/github/>.
3. Find `measles-caputo-stochastic-seir` in the list and switch the toggle
   **On**.
4. Back on GitHub, go to your repository → **Releases** → **Create a new
   release**.
   * Tag: `v1.0.0`
   * Title: `Version accompanying the JNSPS submission`
   * Click **Publish release**.
5. Within a few minutes Zenodo archives the release and issues a DOI. It
   appears on the Zenodo page for the repository, in the form
   `10.5281/zenodo.XXXXXXX`.

---

## Step 4 — put the DOI into the manuscript

In `_body.tex`, the Data Availability section currently reads

```
\url{https://doi.org/10.5281/zenodo.XXXXXXX}
(archived release; the DOI will be activated on acceptance).
```

Replace `XXXXXXX` with the number Zenodo gave you and delete the
parenthetical remark, so that it reads

```
\url{https://doi.org/10.5281/zenodo.1234567}.
```

Then recompile the manuscript and regenerate both PDFs.

**Do not submit while the placeholder is still there.** The editor has
stated that manuscripts containing references whose links do not resolve to
the cited item will be withdrawn from consideration, and an unresolved DOI
in the Data Availability section falls squarely under that warning.

---

## Checklist before you submit

- [ ] Repository is public and the URL opens in a private browsing window
- [ ] `python src/figure_sudan.py` runs on a clean clone
- [ ] Zenodo release created and DOI issued
- [ ] DOI substituted in the manuscript and the manuscript recompiled
- [ ] Both PDFs (clean and marked) regenerated after the substitution
