# MIT Maker Portfolio demo bundle

## Included files
- `kp_convergence_temporal_demo.ipynb`: browser-runnable notebook.
- `toy_genomes_synthetic.csv`: 16 synthetic teaching rows.
- `README_DEMO.md`: these instructions.

## Run in Google Colab (no local installation)
1. Open https://colab.research.google.com/
2. Select **File → Upload notebook** and choose `kp_convergence_temporal_demo.ipynb`.
3. Run cells in order. If the CSV is not present, the first cell prompts you to upload `toy_genomes_synthetic.csv`.
4. Select **Runtime → Run all** to execute the demo.
5. Save a copy to Drive if you want to preserve the executed notebook.

Colab usually includes pandas, NumPy, matplotlib, and scikit-learn. If an import fails, run `%pip install pandas numpy matplotlib scikit-learn` in a separate cell and restart the runtime.

## Add it to your GitHub repository
Upload both files under `demo/`. Link the notebook from the main README. A Colab URL can be built as:
`https://colab.research.google.com/github/Sourcify-y/kp-convergence-prediction/blob/main/demo/kp_convergence_temporal_demo.ipynb`

GitHub can preview the notebook; Colab is the browser execution environment.

## Scientific-integrity statement
The dataset and label are synthetic and illustrative. The notebook does not use real *Klebsiella pneumoniae* data, run HyperTraPS-CT, perform genomic annotation, or reproduce the paper's results. Never quote its metrics as evidence about resistance–virulence convergence. A genuine research run needs the actual frozen inputs, cohort/split manifest, documented preprocessing, model settings, dependency versions, and evaluation protocol.
