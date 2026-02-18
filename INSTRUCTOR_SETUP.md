# Instructor Setup: Result Submission Link

To allow students to upload their results, create a Google Form and add the link to the website.

---

## Step 1: Create a Google Form

1. Go to **https://forms.google.com**
2. Click **Blank** (or use a template)
3. Add these questions:

| Question | Type | Settings |
|----------|------|----------|
| Name | Short answer | Required |
| Project (e.g. Epilepsy, Leukemia) | Short answer | Required |
| Number of significant DEGs found | Short answer | Optional |
| Brief description (2–3 sentences about your findings) | Paragraph | Required |
| Upload: volcano_plot.png | File upload | Required, 10 MB max |
| Upload: pathway_barplot.png | File upload | Required, 10 MB max |
| Upload: significant_deg.csv | File upload | Required, 10 MB max |

4. Click **Send** → copy the form link (e.g. `https://forms.gle/xxxxx`)

---

## Step 2: Add the Link to the Website

1. Open **submit.html** in the repository
2. Find the line: `href="#"` in the "Submit Results Here" button
3. Replace `#` with your actual Google Form link (e.g. `https://forms.gle/xxxxx`)
4. **Email fallback:** If students report they can't upload (e.g. form doesn't support files), find the "Submit via Email" section and replace `YOUR_EMAIL@university.edu` in the mailto link with your email
5. Save and push to GitHub

---

## Alternative: Microsoft Forms

If you use Microsoft 365:

1. Go to **https://forms.office.com**
2. Create a form with the same questions as above
3. Add **File upload** questions (if available in your plan)
4. Copy the form link and paste it into submit.html

---

## Alternative: Email Submission

If you prefer email, change the Submit button in submit.html to:

```html
<a href="mailto:YOUR_EMAIL@university.edu?subject=RNA-seq%20Tutorial%20Results&body=Name:%20%0AProject:%20%0ADescription:%20%0A%0APlease%20attach:%20volcano_plot.png,%20pathway_barplot.png,%20significant_deg.csv" class="btn">Submit via Email</a>
```

Replace `YOUR_EMAIL@university.edu` with your email.
