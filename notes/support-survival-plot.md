# Correct ALT-support survival figures (#613)

Show a visible point when a category contains one supported variant. Aggregate
equal depths and plot the number of calls at or above each unique observed depth;
the curve must drop immediately above a threshold, not at the next observed one.
Keep zero-support calls in the category totals without putting zero on a log axis.
Show an explicit empty state when no positive support exists. Verify singleton,
tied and zero-only cases numerically and inspect the rendered PDF figures.
