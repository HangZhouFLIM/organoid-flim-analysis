# Reproduction environment

Use Python 3.10; the verified interpreter is Python 3.10.9 on Windows. `requirements-lock.txt` pins all 20 installed dependencies, including transitive dependencies. `tested_environment.json` records the checked versions and the lock-file checksum. Install with `python -m pip install -r environments/requirements-lock.txt` and check consistency with `python -m pip check`.

NumPy, pandas and SciPy implement numerical summaries and statistics; Matplotlib generates PNG/SVG displays; Pillow and tifffile read image assets; scikit-image supports the supplied segmentation functions. These versions describe the portable reader environment. They do not establish the historical LAS X, GSLab, Cellpose or model-weight versions used to create the original exports.

The reader uses deposited accepted masks. No additional segmentation-model download is needed to reproduce the configured measurements. Missing historical vendor settings and model provenance are recorded as gaps rather than inferred from the current Python installation.
