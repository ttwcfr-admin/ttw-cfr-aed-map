# TTW CFR defibrillator map

A public web map of defibrillator coverage across TN1–TN4 and TN9–TN12. It shows where new 24/7 units, or opening existing units 24/7, would bring the most residents within reach. It rebuilds itself every time you upload a new register file.

## One-off setup (about 10 minutes, free)

1. Create a free account at https://github.com/signup. A group email such as `ttwcfr@…` is best, so the account isn't tied to one person.
2. Create a new **public** repository called `ttw-cfr-aed-map` (top-right **+** → **New repository**). Don't add a README.
3. Upload the files: on the empty repository page click **uploading an existing file**. Drag in **everything inside** the `ttw-cfr-aed-map` folder (`data`, `docs`, `scripts`, `.github`, `README.md`, `requirements.txt`), then click **Commit changes**.
   *On a Mac, press Cmd+Shift+. in Finder so the hidden `.github` folder is visible before you drag.*
4. Turn on the website: **Settings → Pages → Build and deployment → Source: Deploy from a branch → Branch: `main`, folder `/docs` → Save**.
5. Let the rebuild job write to the repo: **Settings → Actions → General → Workflow permissions → Read and write permissions → Save**.

After a minute or two the map is live at `https://<your-account>.github.io/ttw-cfr-aed-map/`. Bookmark it and share it.

## Every month

1. Download the latest national defibrillator extract (the `.xlsx` or `.csv` file).
2. In the repository, open the `data/defibrillators` folder → **Add file → Upload files** → drop the file in → **Commit changes**.
   * Put the extract date in the file name if you can (for example `defibrillator_data_2026-11-01.xlsx`). The map header and month-by-month trend use it.
   * GitHub's web upload is limited to 25 MB per file. If the extract is bigger, zip it first and upload the `.zip`.
3. Wait about 3 minutes. The **Actions** tab shows a green tick when the rebuild has finished. Refresh the map.

The rebuild uses the newest file, deletes the older ones (they stay in the repository history), and adds a row to `docs/data/history.csv`. The "Month by month" panel tracks coverage from that file.

If the Actions tab shows a red cross, open the failed run. The last lines say what went wrong, usually a missing column.

## Using the data in ArcGIS Online or Google My Maps

Every rebuild publishes these files, which always hold the latest version:

| File | Contents |
|---|---|
| `…/data/defibrillators.geojson` / `.csv` | Every registered defibrillator in and around the area, with access type |
| `…/data/priority_sites.geojson` / `.csv` | Recommended new 24/7 locations, ranked |
| `…/data/quick_wins.geojson` / `.csv` | Existing units to open 24/7, ranked |
| `…/data/postcode_districts.geojson` | TN1–TN4, TN9–TN12 outlines (approximate) |
| `…/data/output_areas.geojson` | Census 2021 Output Areas with residents and density |
| `…/data/need.png` | The "where a new AED would help most" heat map image |

`…` is your site address, for example `https://ttwcfr.github.io/ttw-cfr-aed-map`.

* **ArcGIS Online (free public account):** Map Viewer → **Add → Add layer from URL** → paste a `.geojson` link. The layer reads the live file, so it updates every month without you doing anything. A free public account can't host data itself, so the GitHub copy does that job.
* **Google My Maps:** **Import** the `.csv` files. My Maps won't refresh on its own, so re-import them each month.

## What's in the repository

* `data/static/` holds population by postcode (Census 2021), postcode district outlines and base-map data. These don't change month to month.
* `data/defibrillators/` holds the register extract you upload.
* `data/brand/` holds the TTW CFR symbol and favicon, taken from the brand guidelines.
* `scripts/build.py` contains the model, plus `map_template.html` for the page.
* `docs/` is the published website, written by the build. Don't edit it by hand.

The model's settings sit at the top of `scripts/build.py`: retrieval distance decay, how much limited-hours units count, the age and deprivation weights, and the number of recommendations.

## Sources

Defibrillators: The Circuit national defibrillator network. Population: ONS Census 2021 (TS001, TS007A) via Nomis. Postcodes: ONS Postcode Directory. Output Area boundaries: ONS. Roads: OS Open Roads. Contains OS data © Crown copyright and database right 2026. Place names and rail: © OpenStreetMap contributors. Map tiles: © OpenStreetMap contributors. Retrieval distance: Deakin et al., *Heart* 2018;104:1339–43.
