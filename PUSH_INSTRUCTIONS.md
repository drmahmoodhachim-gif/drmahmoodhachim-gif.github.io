# How to Push to GitHub (One-Time Setup)

The repository **drmahmoodhachim-gif.github.io** does not exist yet on GitHub. Follow these steps:

---

## Step 1: Create the Repository on GitHub

1. Go to **https://github.com/new**
2. **Repository name:** `drmahmoodhachim-gif.github.io` *(must be exact for GitHub Pages)*
3. **Description:** RNA-seq tutorial for PhD students
4. Choose **Public**
5. **Do NOT** add README, .gitignore, or license (we already have them)
6. Click **Create repository**

---

## Step 2: Push Your Local Files

Open **Command Prompt** or **PowerShell** and run:

```bash
cd C:\Users\chime\drmahmoodhachim-gif.github.io
git push -u origin main
```

If prompted for login, use your GitHub username and a **Personal Access Token** (not your password).

---

## Step 3: Enable GitHub Pages (if needed)

1. On your repo page: **Settings** → **Pages**
2. **Source:** Deploy from a branch
3. **Branch:** main → / (root)
4. Click **Save**

---

## Your Site Will Be Live At:

**https://drmahmoodhachim-gif.github.io**

Students can visit this URL to download the tutorial and read the instructions.
