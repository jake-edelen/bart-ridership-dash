# Finished BART Ridership Dashboard Plan

## Summary

Build toward a final portfolio-grade app with a **processed-first architecture**: raw BART files are large and externally distributed through **GitHub Releases**, preprocessing creates compact app-ready files, and Dash reads processed outputs only. The finished app will support route-year-aware maps, monthly ridership selection, monthly timelapse playback, and station-level ridership forecasting.

## Target Architecture

- Make

   

  scripts/prepare_data.py

   

  the source of truth for data processing.

  - Raw hourly OD CSVs and raw Excel workbooks stay in data/raw/ and remain ignored by Git.
  - Processed app-ready files are written to data/processed/ and selected outputs are committed.
  - Dash should not parse raw Excel, raw KML/GDB, or giant hourly CSVs during normal startup.

- Use hourly OD data as the canonical ridership source where available.

  - Monthly station ridership = hourly OD grouped by Origin, summed by year/month/station.
  - January 2018 uses the hourly OD source, not the workbook average weekday fallback.
  - Excel Total Trips / Total Trips OD workbooks become validation/fallback sources, not the preferred source.

- Add a small raw-data manifest for reproducibility.

  - Include expected filenames, years, release URLs, file sizes/checksums if practical.
  - Add a download/setup path using GitHub Releases, but keep the app runnable from committed processed data.

## Implementation Milestones

1. **Data Pipeline Hardening**

   - Refactor load_app_data() so app startup reads processed station, route, and ridership files only.
   - Generate station_ridership_monthly_summary from hourly OD CSVs grouped by Origin.
   - Keep columns explicit: year, month, period, station_code, station_name, ridership, source_type.
   - Add workbook-vs-hourly validation reports for monthly totals and station totals.
   - Update .gitignore so raw data remains ignored but approved processed outputs can be committed.

2. **Dashboard Polish**

   - Finalize top-row layout so the black map, legend, and controls do not visually fight.
   - Rename labels clearly:
     - Route Service Year
     - Ridership Period
     - Selected Route
     - Selected Stations
   - Replace the bulky description box with a compact note explaining:
     - route geometry comes from service year
     - ridership comes from selected ridership period
   - Keep current visual roles:
     - top map = station ridership over route context
     - colored map = route reference
     - bar chart = top 10 or selected station comparison

3. **Monthly Ridership + Timelapse**

   - Keep route service year and ridership period separate by default.
   - Add monthly timelapse controls:
     - period slider
     - play/pause
     - speed control if simple
   - Timelapse changes ridership period month-by-month.
   - Route service year defaults to matching the current frame’s year during playback, unless user explicitly locks route year.
   - Preserve selected route/station filters during playback.

4. **Station-Level Forecasting**

   - Add offline model training, not training inside Dash.

   - First forecast target: monthly station ridership.

   - Build

      

     scripts/train_forecast.py

      

     to create:

     - station-level forecast file
     - model metrics file
     - feature importance or coefficient summary if model supports it

   - Use a defensible model stack:

     - seasonal naive baseline
     - scikit-learn tabular model with lag features, rolling means, month, year, station encoding

   - App shows actual vs forecast for selected station(s), with model metric summary included in the UI or docs.

5. **Documentation + Portfolio Readiness**

   - Update README.md with purpose, screenshots, run steps, app controls, architecture, and limitations.
   - Update DATA_SETUP.md with GitHub Release raw-data download instructions.
   - Add a data methodology section:
     - OD matrix meaning
     - entry-station metric
     - hourly OD aggregation by Origin
     - workbook validation role
   - Keep notebooks/data_audit.ipynb as an explainable audit notebook with monthly totals, source checks, and visual sanity checks.

## Test Plan

- Data pipeline tests:
  - app can start with data/raw/ absent if processed files exist.
  - processed monthly ridership contains no Entries, Exits, Grand Total, or null station rows.
  - January 2018 monthly total from hourly OD equals expected grouped total.
  - workbook total-trip validation is within expected tolerance for overlapping months.
- Dashboard tests:
  - route dropdown includes Oakland Airport connector.
  - ridership year/month controls update figure titles and marker values.
  - timelapse callback updates ridership period without changing selected station filters.
  - marker sizes remain capped for high-ridership months.
- Route tests:
  - service-year route options still match known BART eras.
  - colored map offsets remain display-only.
  - route geometry still passes station-sequence validation.
- Forecast tests:
  - feature table has expected lag/rolling columns.
  - train script writes forecasts and metrics.
  - model beats or is clearly compared against seasonal naive baseline.
  - app reads forecast outputs without retraining.

## Assumptions And Defaults

- Finished app means dashboard polish, processed monthly ridership, monthly timelapse, and station-level forecasting.
- Raw hourly OD CSVs are distributed through GitHub Releases, not committed to Git.
- App runtime is processed-first and should not require raw files for normal use.
- Canonical station ridership metric is **entries by origin station**.
- Monthly timelapse is the final polished animation scope; daily/hourly playback remains future work.
- Forecasting is station-level monthly ridership, trained offline and displayed from processed prediction files.



Yes, that is exactly where OD summaries become useful.

A flow feature would answer:

```text
Where are riders coming from?
Where are they going?
How does that movement change over time?
```

**Possible Flow Visuals**
1. **Selected origin -> all destinations**
   User selects `Embarcadero` as origin.
   Map shows lines from Embarcadero to destination stations, with line width/color based on trip volume.

2. **All origins -> selected destination**
   User selects `SFO`.
   Map shows where airport-bound riders are coming from.

3. **Selected origin + destination**
   User selects `19th St/Oakland -> Embarcadero`.
   Chart/timelapse shows flow volume over:
   - hours in a day
   - days in a week
   - months in a year
   - years over time

4. **Top OD pairs**
   App shows top 10 station pairs for selected period.

5. **Commuting rhythm view**
   Compare AM vs PM directionality:
   - morning: East Bay -> SF
   - evening: SF -> East Bay

**Timelapse Version**
A user could select:

```text
Origin: 19th St/Oakland
Destination: Embarcadero
Range: January 2018
Granularity: Hour
```

Then the app plays through hours/days and updates:

```text
line thickness = trips for that hour
marker size = station entries/exits
bar chart = top destinations from selected origin
```

Or for one day:

```text
6 AM -> low flow
8 AM -> strong 19th -> Embarcadero commute flow
5 PM -> reverse or different flow pattern
```

**Modeling Flows**
OD data also enables better modeling than station totals alone.

Station model:

```text
Predict total entries at Embarcadero next month.
```

OD model:

```text
Predict trips from 19th St/Oakland to Embarcadero next hour/day/month.
```

Potential features:
- origin station
- destination station
- hour
- day of week
- month
- route/service era
- lagged OD ridership
- rolling average OD ridership
- holiday/event flags later

Models could be:
- baseline seasonal average
- gradient boosting / XGBoost
- temporal neural net later if justified

**What Data Table You’d Need**
For flow timelapse, you’d want Parquet like:

```text
od_ridership_hourly.parquet
```

Columns:

```text
Date
Year
Month
DayOfWeek
Hour
Origin
Origin Name
Destination
Destination Name
Ridership
```

Then you can filter quickly:

```text
Origin == selected_origin
Destination == selected_destination
Date between start and end
Hour == selected_hour
```

**Scope Advice**
For now, do not build full OD flow animation yet.

Best staged path:

```text
1. Build station_hourly summaries.
2. Build OD_monthly summaries.
3. Add a static OD analysis view: selected origin -> top destinations.
4. Then add time controls.
5. Then consider animated flow.
```

That gives you value before committing to a complex animation system.





Yes. OD should be a separate tab.

If you put station ridership, route reference, time lapse, and OD flows all on one screen, the dashboard will become confusing fast. These are related, but they answer different questions.

**Recommended App Structure**
Use top-level tabs:

```text
1. System Overview
2. Station Trends
3. OD Flows
4. Data Quality / Methodology
```

**Tab 1: System Overview**
Purpose:

```text
What does BART ridership look like across the network for a selected period?
```

This is closest to your current dashboard.

UI:
- Ridership year
- Ridership month
- Route selector
- Station comparison selector

Views:
- station ridership map
- colored route reference map
- top 10 / selected station comparison bar chart

Keep this clean and polished.

**Tab 2: Station Trends**
Purpose:

```text
How does ridership at one or more stations change over time?
```

UI:
- station multi-select
- date range
- granularity: monthly / daily / hourly
- optional route/service era overlay

Views:
- time series line chart
- weekday/weekend comparison
- hour-of-day heatmap
- station ranking table

This is where your station-hour Parquet summaries pay off.

**Tab 3: OD Flows**
Purpose:

```text
Where are riders traveling from and to?
```

UI:
- mode selector:
  - origin -> destinations
  - origins -> destination
  - origin + destination pair
  - top OD pairs
- origin dropdown
- destination dropdown
- date range
- granularity: monthly / daily / hourly
- top N slider

Views:
- flow map with lines
- top destination/origin bar chart
- OD pair time series
- optional flow table

For first implementation, avoid animation. Start with static filtered OD flows.

**Tab 4: Data Quality / Methodology**
Purpose:

```text
Can this project explain and defend its data choices?
```

UI:
- source comparison selector
- threshold slider for differences
- year/month selector

Views:
- hourly OD vs workbook validation table
- station-level difference chart
- source-of-truth note
- known limitations

This is underrated for portfolio value. It shows you audited your own pipeline.

**Where Time Lapse Fits**
I would not make time lapse its own first-class tab yet. Add it later inside:

```text
Station Trends
OD Flows
```

Example:
- Station Trends: play through monthly/daily/hourly station map.
- OD Flows: play through selected OD flows over time.

**Current Dashboard Integration**
Your current app should become `System Overview`.

The top-level UX should look like:

```text
[System Overview] [Station Trends] [OD Flows] [Data Quality]
```

This lets each feature have a clear job.

**Best Build Order**
1. Polish current dashboard as `System Overview`.
2. Add `Data Quality / Methodology` tab using the validation CSV.
3. Add `Station Trends` using station monthly/daily/hourly summaries.
4. Add `OD Flows` with static OD monthly first.
5. Add time lapse controls only after static views are stable.

That sequence keeps the app coherent and prevents the UI from becoming a pile of controls.







ML should be a later tab, not mixed into the core dashboard.

Recommended structure eventually:

```text
1. System Overview
2. Station Trends
3. OD Flows
4. Forecasting
5. Data Quality / Methodology
```

**Where ML Fits**
The ML/algorithms piece belongs in `Forecasting`.

Purpose:

```text
Predict future ridership and explain what historical patterns drive the forecast.
```

It should not be the first thing users see. First, the app should establish:

```text
What happened?
Where did it happen?
How did flows change?
Can the data be trusted?
```

Then ML answers:

```text
What might happen next?
```

**What To Model First**
Start simple and defensible.

Best first target:

```text
Monthly station ridership forecast
```

Example:

```text
Predict next 3-12 months of ridership for selected station(s).
```

Why:
- simpler than hourly OD
- easier to validate
- visually clear
- good portfolio value

Avoid starting with LSTM. A strong sequence would be:

```text
1. Seasonal naive baseline
2. Linear / Ridge regression with calendar features
3. Random Forest or XGBoost
4. Optional LSTM only after baselines
```

**Forecasting Tab UI**
Controls:
- target type:
  - systemwide ridership
  - station ridership
  - OD pair ridership later
- station selector
- forecast horizon
- model selector:
  - seasonal baseline
  - regression
  - XGBoost
  - LSTM later
- train/test split display

Views:
- actual vs predicted line chart
- forecast interval if available
- error metrics:
  - MAE
  - RMSE
  - MAPE or sMAPE
- feature importance for tree/regression models
- residual plot

**How It Connects To Existing Work**
Your preprocessing layer feeds the ML tab:

```text
station_ridership_monthly.parquet
station_ridership_daily.parquet
od_ridership_monthly.parquet later
```

Your validation work supports the model credibility:

```text
source of truth is defined
station mappings are cleaned
monthly totals are audited
```

**OD Modeling Later**
After station forecasting works, you can model OD flows:

```text
Predict trips from Origin -> Destination over time.
```

That is more advanced and more impressive, but much heavier.

Good OD modeling targets:
- top 50 OD pairs only
- selected route corridor
- airport-bound trips
- commute-direction flows

**Portfolio Advice**
A clean forecasting tab with a baseline + XGBoost, tested against held-out months, is stronger than an LSTM demo with unclear assumptions.

A good resume bullet later:

> Built a forecasting module for station-level BART ridership using lag, rolling-average, and calendar features, benchmarking baseline seasonal models against gradient-boosted regressors with MAE/RMSE evaluation.

That sounds much more credible than:

> Used LSTM to predict ridership.

Because it shows process, not just algorithm name.