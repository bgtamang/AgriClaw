# Demo data

The three JSON files are unmodified responses from the USDA-NRCS Soil Data Access REST service (`https://sdmdataaccess.sc.egov.usda.gov/tabular/post.rest`, format `JSON+COLUMNNAME`), saved on 2026-10-15. The SQL is the same text that `scripts/soil.py` sends (`SQL_POINT`, `SQL_MAPUNIT`, `SQL_HORIZONS`):

- `sda_point_40.08_-88.20.json`: the map unit key at 40.08°N, 88.20°W, near Champaign, Illinois (mukey 242963).
- `sda_mapunit_242963.json`: map unit, survey area and `muaggatt` row for Drummer silty clay loam, 0 to 2 percent slopes (Champaign County, IL019, survey version 2026-08-23).
- `sda_horizons_242963.json`: all components (Drummer 94%, Peotone 3%, Harpster 3%) with horizons and the representative texture class.

`tests/fixtures/` holds the same files plus responses for mukey 426520 (Leola loamy sand, Portage County, Wisconsin; includes organic horizons) and an empty response (`{}`) for a point in Lake Michigan.

License: SSURGO is a work of the US Government (USDA-NRCS) and is in the public domain. Citation: Soil Survey Staff, Natural Resources Conservation Service, United States Department of Agriculture. Soil Survey Geographic (SSURGO) Database. Available online. Accessed 2026-10-15.
