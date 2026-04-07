# EpiSheet — UI design specification

> This document supersedes the Excel colour conventions described in CLAUDE.md.
> The UI is a modern card-based web interface, not a visual replica of the spreadsheet.
> All design decisions here were finalised through iterative mockup review.

---

## Colour palette

| Name | Hex | Usage |
|------|-----|-------|
| Petrol | `#005f6a` | Primary brand colour, page titles, active states, primary result card |
| Blue | `#165c7d` | Top navigation bar, secondary actions, result labels |
| Petrol light | `#e8f5f7` | Input panel header background |
| Blue light | `#edf2f8` | Results panel header background, secondary metric cards |
| Teal light | `#e6f4f5` | Tertiary metric cards (e.g. χ² test) |
| Purple light | `#f0ecfb` | Quaternary metric cards (e.g. homogeneity) |
| Green light | `#f5f9f0` | Interpretation block background |
| Green border | `#b8d9a0` | Interpretation block border |
| Cell border | `#d0dde3` | 2×2 table cell borders |
| Cell background | `#eef3f5` | 2×2 table header/total cell backgrounds |

---

## Typography

- Font: system default (inherit from browser/Shiny)
- Page title: 18px, font-weight 500, colour `#005f6a`
- Page subtitle: 12px, muted, line-height 1.5
- Panel headers: 11px, font-weight 500, uppercase, letter-spacing 0.06em
- Result labels: 10px, font-weight 500
- Result values: 22px, font-weight 500
- Result sub-text (CI, unit): 10px, muted
- Table headers: 11px, font-weight 500, colour `#4a6070`
- Breadcrumb: 11px, muted

---

## App structure

### Two-screen flow

**Screen 1 — Landing page (calculator selection)**
- Full-width card grid, no sidebar
- Cards grouped by category (Stratified analysis / Sample size & power / Tools & advanced)
- Each category has a colour family: teal, blue, purple
- Each card: icon (30×30px, rounded 7px, coloured background) + calculator name (12px, 500) + short description (10px, muted)
- Clicking a card navigates to Screen 2

**Screen 2 — Calculator**
- "← All calculators" back link at top
- Page title (petrol) + subtitle (muted)
- Two-column layout: input panel (left) + results panel (right)
- Interpretation block full-width below both panels

### Top navigation bar
- Background: `#165c7d`
- Left: brand icon + "EpiSheet" (15px, 500) + "Based on Episheet by K. Rothman" (10px, muted, pushed to right via margin-left: auto)
- Right of brand ref: "Glossary" button (white text, semi-transparent background, book icon)

---

## Components

### Input panel

```css
.panel {
  border: 0.5px solid #d0dde3;
  border-radius: 12px;
  overflow: hidden;
}
.panel-head {
  padding: 9px 14px;
  border-bottom: 0.5px solid #d0dde3;
  font-size: 11px;
  font-weight: 500;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  background: #e8f5f7;   /* teal for input panels */
  color: #005f6a;
}
.panel-body {
  padding: 14px;
}
```

### 2×2 table

```css
.tbl {
  border-collapse: separate;
  border-spacing: 3px;
  width: 100%;
  margin-bottom: 12px;
}
.tbl th {
  font-size: 11px;
  font-weight: 500;
  color: #4a6070;
  padding: 6px 10px;
  text-align: center;
  background: #eef3f5;
  border: 1px solid #d0dde3;
  border-radius: 6px;
}
.tbl th.rh { text-align: left; }
.tbl td {
  border: 1px solid #d0dde3;
  border-radius: 6px;
  padding: 0;
  text-align: center;
  vertical-align: middle;
}
.tbl td.cl {
  /* row label cells */
  font-size: 11px;
  color: #4a6070;
  background: #eef3f5;
  padding: 6px 10px;
  text-align: left;
}
.tbl td.tot {
  /* total cells */
  background: #eef3f5;
  font-size: 12px;
  font-weight: 500;
  color: #4a6070;
  padding: 6px 8px;
  border: 1px solid #d0dde3;
  border-radius: 6px;
}
/* cell letter annotations (a, b, c, d) */
.cell-letter {
  position: absolute;
  top: 3px;
  left: 5px;
  font-size: 9px;
  color: #005f6a;
  opacity: 0.5;
  font-weight: 500;
}
/* input fields inside cells */
.tbl input {
  width: 58px;
  font-size: 13px;
  border: none;
  background: #ffffff;
  color: #1a1a1a;
  text-align: center;
  outline: none;
  padding: 4px 6px;
}
.tbl input:focus {
  background: #eef3f5;
  border-radius: 3px;
}
```

Key point: use `border-collapse: separate` + `border-spacing: 3px` to allow `border-radius` on individual cells. This creates the rounded, card-like cell appearance.

### Calculate button

```css
.calc-btn {
  width: 100%;
  padding: 10px;
  background: #ffffff;
  color: #165c7d;
  border: 1.5px solid #165c7d;
  border-radius: 8px;
  font-size: 13px;
  font-weight: 500;
  cursor: pointer;
  letter-spacing: 0.02em;
}
.calc-btn:hover {
  background: #edf2f8;
}
```

The button is outline style (white background, blue border, blue text) — not solid filled. This distinguishes it from the results and gives it visual weight without overwhelming the form.

### Results panel

```css
.panel-head.ph-blue {
  background: #edf2f8;
  color: #165c7d;
}
```

Results use a 2-column metric grid:

```css
.metrics {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 8px;
}
.metric {
  padding: 11px 13px;
  border-radius: 8px;
}
```

**Result card layout (rate data example):**

| Position | Metric | CSS class | Spans |
|----------|--------|-----------|-------|
| Row 1, full width | Rate ratio (MH) | `m-primary` | `grid-column: span 2` |
| Row 2, left | Rate difference | `m-blue` | 1 column |
| Row 2, right | χ² test (Woolf) | `m-teal` | 1 column |
| Row 3, left | Exposed rate | `m-blue` | 1 column |
| Row 3, right | Homogeneity χ² | `m-purple` | 1 column |

**Do not use `wide` / `grid-column: span 2` on homogeneity** — it should sit beside Exposed rate.

Metric card colours:

```css
.m-primary { background: #005f6a; }   /* petrol — main result */
.m-blue    { background: #e8f0f7; }   /* blue tint */
.m-teal    { background: #e6f4f5; }   /* teal tint */
.m-purple  { background: #f0ecfb; }   /* purple tint */
```

Text colours per card type:

| Class | Label colour | Value colour | Sub-text colour |
|-------|-------------|--------------|-----------------|
| `m-primary` | `rgba(255,255,255,0.7)` | `#fff` | `rgba(255,255,255,0.55)` |
| `m-blue` | `#165c7d` | `#165c7d` | `#7a9bbf` |
| `m-teal` | `#005f6a` | `#005f6a` | `#4a9aa5` |
| `m-purple` | `#6b3fa0` | `#6b3fa0` | `#9b6fc8` |

### Info markers (ⓘ)

Inline after each result label. Very small, same colour as the label but at 55% opacity.

```css
.info-mark {
  font-size: 9px;
  cursor: pointer;
  opacity: 0.55;
  background: none;
  border: none;
  padding: 0 0 0 2px;
  line-height: 1;
}
.info-mark:hover { opacity: 1; }
```

Tooltip appears below the marker, 190px wide, plain white background with a subtle border. Only one tooltip visible at a time; click anywhere outside to dismiss.

### Interpretation block

Full-width, below the two-column layout. Green tint to signal "conclusion".

```css
.conclusion {
  margin-top: 16px;
  background: #f5f9f0;
  border: 0.5px solid #b8d9a0;
  border-radius: 12px;
  padding: 14px 16px;
}
.conc-head {
  font-size: 11px;
  font-weight: 500;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  color: #3b6d11;
  margin-bottom: 8px;
}
.conc-text {
  font-size: 12px;
  line-height: 1.6;
}
```

The interpretation text is **generated dynamically** from the calculated values. It reads the numbers and writes a plain-English paragraph explaining what they mean. Key numbers are bolded inline. This is a core feature for educational use — students can check their interpretation against the app's.

### Glossary panel

Triggered by the "Glossary" button in the top bar. Slides in from the right as a fixed overlay panel (300px wide). Contains plain-language definitions of every result term, including general concepts like "person-years" and "95% CI". Each entry: term (12px, petrol, 500) + definition (11px, muted, line-height 1.6).

---

## Layout constraints

- `max-width: 900px` on the app container, centred
- Two-column main layout: `grid-template-columns: 1fr 1fr; gap: 16px`
- Do not use a sidebar for navigation — all navigation is via the landing page card grid and the top bar

---

## Landing page card grid

```css
.card-grid {
  display: grid;
  grid-template-columns: repeat(4, minmax(0, 1fr));
  gap: 8px;
  margin-bottom: 20px;
}
.lcard {
  background: white;
  border: 0.5px solid #d0dde3;
  border-radius: 12px;
  padding: 12px 13px;
  cursor: pointer;
}
.lcard:hover { border-color: #a0b8c8; }
.lcard-icon {
  width: 30px;
  height: 30px;
  border-radius: 7px;
  margin-bottom: 8px;
}
/* Icon background colours by category */
.ic-teal   { background: #cce9ec; }  /* Stratified analysis */
.ic-blue   { background: #d0e6f5; }  /* Sample size & power */
.ic-purple { background: #e8e4f8; }  /* Tools & advanced */
.ic-gray   { background: #f0f0ec; }  /* Coming soon / misc */

.lcard-name { font-size: 12px; font-weight: 500; margin-bottom: 4px; }
.lcard-desc { font-size: 10px; color: #6a7f8a; line-height: 1.4; }
```

Group labels above each row of cards: 10px, uppercase, letter-spacing 0.08em, muted.

---

## What NOT to do

- Do not use yellow/orange/red cell colours from the Excel original — those conventions are replaced by this design system
- Do not use a sidebar for calculator navigation
- Do not show "Sheet 1", "Sheet 2" etc. — use the calculator name only
- Do not use `border-collapse: collapse` on the 2×2 table — rounded corners require `border-collapse: separate`
- Do not make the homogeneity χ² card span full width — it should sit beside "Exposed rate"
- Do not use a solid filled Calculate button — use the outline style
- Do not hardcode CSS variable names like `--color-border-tertiary` in files meant to run in a plain browser — use concrete hex values
- Do not use a `navbarPage()` dropdown as the primary navigation — use the landing page card grid instead

---

## Reference file

The HTML prototype (`episheet_calculator_v8.html`) is the design reference.
All Shiny modules should match its visual output.
