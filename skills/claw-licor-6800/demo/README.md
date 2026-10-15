# Demo data

`--demo` reads two LI-6800 Excel logs that ship with the PhotoGEA R package
(MIT licence), loaded in place with `system.file("extdata", package = "PhotoGEA")`
rather than copied here:

- `c3_aci_1.xlsx`: three 16-point A-Ci curves (soybean and tobacco) logged on
  instrument `ripe1`; the file header carries console and head serial numbers.
- `c3_aci_2.xlsx`: three curves logged on instrument `ripe2`; PhotoGEA's authors
  deleted its Remarks sheet (so there are no serial numbers) and changed the
  oxygen value after observation 10 as a reader test.

PhotoGEA's documentation states that noise was added to these unpublished
measurements, so no biological conclusions should be drawn from them. The demo
uses `--curve-cols instrument,species,plot --design-cols species --instrument-col instrument`.

The tests also read PhotoGEA's `plaintext_licor_file` and
`plaintext_licor_file_v2` (a log closed and reopened mid-measurement), and build
an "unevaluated Excel" copy of `c3_aci_1.xlsx` in a temporary folder.
