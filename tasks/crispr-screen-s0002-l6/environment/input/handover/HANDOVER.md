# Hand-over, VAL-2 drug screen (K.B.)

Sorry for the state of this, my last week was short.

- `count_guides.py` - my counting script. `counts_draft.tsv` is its output for the eight libraries.
- `analyse_screen.py` - gene table and hit list from a count table and the sample sheet.
  `gene_draft.tsv` and `hits_mailed.tsv` are its output for `counts_draft.tsv`.
- `hits_mailed.tsv` went to the PI by mail on Friday as "top 40 modifier candidates". Nobody has
  ordered reagents on it yet as far as I know.

Reads assigned per library in `counts_draft.tsv`:

- LIB6: 111279
- LIB7: 115554
- LIB8: 101881
- LIB9: 106875
- LIB10: 108057
- LIB11: 107956
- LIB12: 0
- LIB13: 113328

Not done: the reference genes, anything on the re-sequenced library beyond counting it
the way the script does. I used the dropout workflow we always use and adapted it for two arms.
