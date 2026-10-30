# Tested batch package, version 2.0

- Provided 3000 × 4096 real pair: B 6545 / R&B 6544; separate raw initial 0001. Completed in about 23 seconds here, including the illustrated report; machine speed varies.
- Fixed pore denominator: 4,298,995 pixels. Nominal liquid 67.1709%; air 32.8291%; review 35.1960%. Parameter sweeps disabled in the batch default. Review fraction therefore differs from the older fully swept single-pair example; nominal labels agree.
- Twelve frame-index planning checks; thirty reference/mask checks; eleven batch workflow checks, including stopped failures, KeyboardInterrupt, continue policy, missing frames and resume/source mismatches.
- Ten actual synthetic pair analyses (118 assertions): raw binarization, exact external binary masks, reversed B parity and single-pair compatibility. Opposite parity produced identical phase labels for the same physical B/R&B images.
- Saturation arithmetic, denominator conservation, review ranges, clipping scope, plotting and atomic CSV failure checks passed.
- Additional extraction/relocation check is recorded in PACKAGE_VERIFICATION.json.

These tests establish software behavior and consistency, not scientific accuracy. Only one acquired experimental pair was provided; no acquired long series has been tested. Synthetic images are test fixtures, not experimental phase ground truth. The method retains ambiguity/review labels and uses 2D projected areas.

[Open the tested report](example_batch_results/run_20261029_124806_581510/pairs/pair_000000_B6545_RB6544/GUIDE.html) · [Open saturation CSV](example_batch_results/run_20261029_124806_581510/saturation.csv)
